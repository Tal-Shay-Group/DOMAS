import logging
import os
import queue
import random
import shutil
import tempfile
import threading
import time
import warnings
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

import utils
from utils import (
    FEATURE_JUNCTION, FEATURE_RETAINED_INTRON, FEATURE_TYPE_COLUMN,
    find_matching_junction_indices,
)
from generate_gene_pdf import GeneVisualization, prepare_gene_data_bulk

# Suppress FutureWarning about DataFrame concatenation behavior
warnings.filterwarnings('ignore', category=FutureWarning, message='.*DataFrame concatenation.*')


logger = logging.getLogger(__name__)

DOMAIN_NAME_COLUMNS = ['interpro', 'pfam', 'cdd', 'smart', 'tigr', 'CDD_id']
DOMAIN_NAME_PREFIX_PRIORITY = ['IPR', 'pfam', 'cd', 'smart', 'tigr', 'CDD']


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Phase 1: junction <-> exon matching
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Phase 1.5: "most like canonical" transcript selection, an alternative to
# picking the longest-CDS transcript
# ---------------------------------------------------------------------------

def _exon_coord_set(df_exons):
    return set(zip(df_exons['genomic_start_tx'], df_exons['genomic_end_tx']))


def _outside_range_exon_set(df_exons, min_bp, max_bp):
    """Exons of df_exons that fall entirely outside [min_bp, max_bp], as a set of
    (genomic_start_tx, genomic_end_tx) coordinate pairs."""
    outside = df_exons[(df_exons['genomic_end_tx'] < min_bp) | (df_exons['genomic_start_tx'] > max_bp)]
    return _exon_coord_set(outside)


def select_longest_cds(candidate_ids, cds_length_by_transcript):
    """
    The longest-CDS transcript among `candidate_ids`. Ties fall to the lowest id,
    so the choice does not depend on Python's per-process string hash seed.

    @param candidate_ids: transcript ids to choose between (must be non-empty)
    @param cds_length_by_transcript: dict {transcript_id: coding length in bases}
    @return: the chosen transcript id
    """
    return max(candidate_ids, key=lambda tid: (cds_length_by_transcript.get(tid, -1), tid))


def select_most_like_canonical(comparable_transcript_ids, canonical_transcript_id, transcript_exons, junctions,
                                cds_length_by_transcript):
    """
    Pick the comparable transcript that is "most like canonical", as an alternative to
    always taking the longest-CDS one when several transcripts have unique junctions vs
    canonical.

    A candidate qualifies when ALL of its exons lying *outside* the cluster's junction
    range (the span from the earliest to the latest coordinate among the cluster's own
    junctions) exactly match - by genomic start/end, as a set - canonical's exons
    outside that range. Such a transcript differs from canonical only within the spliced
    region, so any domain difference it shows is attributable to the event itself.

    Qualifying is a hard filter, not a preference: a single qualifying candidate is
    taken even if a disqualified one has a longer CDS. Among several qualifying
    candidates the longest-CDS one is taken, ties broken by transcript id for
    determinism. `comparable_transcript_ids` is expected to have been reduced to the
    protein-coding candidates already (step 1 of the priority) - this function
    implements steps 2 and 3 only.

    Returns None when no candidate qualifies - the flag then goes unset rather than
    falling back to the longest-CDS candidate, which is_longest_cds already marks in
    its own right. Callers wanting a single transcript per cluster fall back to that
    flag themselves (see results_stats.select_representative_transcript()).
    """
    # min/max over BOTH coordinates of each pair: normalize_junctions_frame()
    # orders them, but this is reachable with a hand-built list too.
    min_bp = min(min(pair) for pair in junctions)
    max_bp = max(max(pair) for pair in junctions)

    c_exons = transcript_exons[canonical_transcript_id]
    c_outside_set = _outside_range_exon_set(c_exons, min_bp, max_bp)

    qualifying = [
        transcript_id for transcript_id in comparable_transcript_ids
        if _outside_range_exon_set(transcript_exons[transcript_id], min_bp, max_bp) == c_outside_set
    ]

    if not qualifying:
        return None

    return max(
        qualifying,
        key=lambda tid: (cds_length_by_transcript.get(tid, -1), tid),
    )


# ---------------------------------------------------------------------------
# Phase 2: domain boundary determination via coordinate matching
# ---------------------------------------------------------------------------

def find_boundary_exons(df_exons, min_bp, max_bp):
    """
    Return (first_exon, last_exon): the exon whose genomic end is closest to
    (but not before) min_bp, and the exon whose genomic start is closest to
    (but not after) max_bp, allowing a 1bp tolerance on both ends.

    This is called several times per compared transcript (Phase 2 window
    refinement), on the small per-transcript exon slice - implemented with
    numpy array ops + a single positional .iloc[] instead of pandas boolean
    filtering + idxmin/idxmax + .loc[], since pandas' per-call overhead
    dominates at this scale and this runs millions of times over a large
    junctions file.
    """
    ends = df_exons['genomic_end_tx'].to_numpy()
    starts = df_exons['genomic_start_tx'].to_numpy()

    first_mask = ends >= min_bp - 1
    if not first_mask.any():
        raise ValueError("attempt to get argmin of an empty sequence")
    first_pos = np.where(first_mask, ends, np.inf).argmin()

    last_mask = starts <= max_bp + 1
    if not last_mask.any():
        raise ValueError("attempt to get argmax of an empty sequence")
    last_pos = np.where(last_mask, starts, -np.inf).argmax()

    first_exon = df_exons.iloc[first_pos]
    last_exon = df_exons.iloc[last_pos]
    return first_exon, last_exon


def _bp_to_cds(exon, bp):
    """Map a genomic bp position within `exon` to its CDS-relative bp offset, clamped to the exon's CDS span."""
    cds = exon['abs_start_CDS'] + (bp - exon['genomic_start_tx'])
    return min(max(cds, exon['abs_start_CDS']), exon['abs_end_CDS'])


def get_aa_range(first_exon, last_exon, min_bp=None, max_bp=None):
    """
    Amino-acid range (inclusive) spanned by the first/last exons of a window.

    If min_bp/max_bp are given, the range is narrowed to the AA positions
    corresponding to those genomic bp coordinates (clamped within the
    first/last exon's CDS span) instead of the exons' full CDS span.

    first_exon/last_exon are chosen by genomic position (the exon nearest the
    window's lower/upper genomic bound, see find_boundary_exons). On a minus
    strand transcript that does NOT mean first_exon precedes last_exon in CDS
    order - the exon nearest the genomically-lower bound can be the LATER
    exon in the transcript. When min_bp/max_bp aren't given we therefore
    can't just use first_exon's start and last_exon's end (those could be two
    CDS-adjacent values straddling a single codon instead of spanning the
    full region); we take the min/max across all four boundary CDS values
    instead, which is correct regardless of strand.
    """
    if min_bp is None and max_bp is None:
        cds_bounds = [
            first_exon['abs_start_CDS'], first_exon['abs_end_CDS'],
            last_exon['abs_start_CDS'], last_exon['abs_end_CDS'],
        ]
        start_cds = min(cds_bounds)
        end_cds = max(cds_bounds)
    else:
        start_cds = first_exon['abs_start_CDS'] if min_bp is None else _bp_to_cds(first_exon, min_bp)
        end_cds = last_exon['abs_end_CDS'] if max_bp is None else _bp_to_cds(last_exon, max_bp)
    min_aa = min(start_cds, end_cds) // 3
    max_aa = max(start_cds, end_cds) // 3
    return min_aa, max_aa


def _cds_to_bp(exon, cds_bp):
    """Map a CDS-relative bp position to its genomic position within `exon`'s CDS span."""
    cds_bp = min(max(cds_bp, exon['abs_start_CDS']), exon['abs_end_CDS'])
    return exon['genomic_start_tx'] + (cds_bp - exon['abs_start_CDS'])


def find_bp_range_for_domains(df_exons, domains_in_region):
    """
    Genomic (start, end) bp positions - start <= end - spanning the start of
    the nearest-starting domain and the end of the furthest-reaching domain in
    `domains_in_region`, or (None, None) if there are no domains with a
    defined AA range.
    """
    # numpy-array filtering instead of pandas boolean-mask DataFrames - this
    # runs (conditionally, in Phase 2's window-refinement round) for every
    # compared transcript, and df_exons/domains_in_region are always tiny
    # per-transcript slices, so pandas' per-call overhead otherwise dominates.
    aa_start = domains_in_region['AA_start'].to_numpy()
    valid = aa_start != 0
    if not valid.any():
        return None, None
    aa_end = domains_in_region['AA_end'].to_numpy()
    min_domain_bp = aa_start[valid].min() * 3
    max_domain_bp = aa_end[valid].max() * 3

    starts = df_exons['abs_start_CDS'].to_numpy()
    ends = df_exons['abs_end_CDS'].to_numpy()
    first_mask = (starts <= min_domain_bp) & (ends >= min_domain_bp)
    last_mask = (starts <= max_domain_bp) & (ends >= max_domain_bp)
    if not first_mask.any() or not last_mask.any():
        return None, None

    # argmax on a bool array returns the position of the first True value,
    # matching the original's `.iloc[0]` of the filtered rows (first match
    # in df_exons' original order).
    first_exon = df_exons.iloc[np.argmax(first_mask)]
    last_exon = df_exons.iloc[np.argmax(last_mask)]

    # On a minus strand the domain's lower CDS bound maps to a higher genomic
    # position, so the pair can come back reversed. Callers pool it into a
    # min/max across transcripts, so return it in genomic (low, high) order.
    bp_a = _cds_to_bp(first_exon, min_domain_bp)
    bp_b = _cds_to_bp(last_exon, max_domain_bp)
    return min(bp_a, bp_b), max(bp_a, bp_b)


def build_filtered_domain_lookup(domain_lookup):
    """Wrap a domain lookup so filter_representative_domains() runs once per
    transcript and its result is kept.

    Returns (lookup, kept), where `lookup(transcript_id)` yields the transcript's
    domains after the ladder and `kept` accumulates {transcript_id: frame} for
    every transcript asked about. The comparison and the drawing then share one
    decision about which entries are domains instead of each reaching it
    separately - the PDF reads `kept` rather than re-running the filter over its
    own copy of the rows.
    """
    kept = {}

    def lookup(transcript_id):
        if transcript_id not in kept:
            kept[transcript_id] = filter_representative_domains(domain_lookup(transcript_id))
        return kept[transcript_id]

    return lookup, kept


def build_exon_lookup(df_exons):
    """
    Precompute, once per analyze_junctions() run, a transcript_id -> exons
    lookup so per-cluster/per-transcript filtering of the full `df_exons`
    DataFrame (the dominant cost of analyze()) is replaced by O(1)
    dict lookups.
    """
    by_ensembl = {tid: g for tid, g in df_exons.groupby('transcript_ensembl_id')}
    by_refseq = {tid: g for tid, g in df_exons.groupby('transcript_refseq_id')}
    empty = df_exons.iloc[0:0]
    # A transcript id that's ambiguous between the ensembl/refseq groupings is
    # looked up once per cluster it appears in - across a gene with many
    # clusters that's the same pd.concat + drop_duplicates() re-run
    # identically every time. Memoize the merged result per transcript id
    # instead of recomputing it on every lookup() call.
    merged_cache = {}

    def lookup(transcript_id):
        a = by_ensembl.get(transcript_id)
        b = by_refseq.get(transcript_id)
        if a is not None and b is not None:
            merged = merged_cache.get(transcript_id)
            if merged is None:
                merged = pd.concat([a, b]).drop_duplicates()
                merged_cache[transcript_id] = merged
            return merged
        return a if a is not None else (b if b is not None else empty)

    return lookup


def build_domain_lookup(df_domains):
    """
    Precompute, once per analyze_junctions() run, a transcript_id -> domains
    lookup so per-cluster filtering of the full `df_domains` DataFrame is
    replaced by O(1) dict lookups.
    """
    by_transcript = {
        tid: g.rename(columns={'transcript_ensembl_id_version': 'transcript_ensembl_id'})
        for tid, g in df_domains.groupby('transcript_ensembl_id_version')
    }
    empty = df_domains.rename(columns={'transcript_ensembl_id_version': 'transcript_ensembl_id'}).iloc[0:0]

    def lookup(transcript_id):
        return by_transcript.get(transcript_id, empty)

    return lookup


def _domains_in_aa_range(df_domains, min_aa, max_aa):
    # Boolean-mask indexing already returns a new, independent DataFrame (not
    # a view), so the extra .copy() here was a redundant allocation on a
    # function called 4x per compared transcript.
    return df_domains[(df_domains['AA_end'] >= min_aa) & (df_domains['AA_start'] <= max_aa)]




# The only InterPro entry types DOMAS considers, from RepresentativeDomains.type
# (sourced from interpro.xml.gz). Domain and Repeat are the curated structural-
# functional units - the things a splicing event can remove.
#
# Everything else is dropped: Family and Homologous_superfamily name what a
# protein IS rather than delimiting a unit within it; Active_site / Binding_site
# / Conserved_site / PTM are residue positions rather than regions; and a
# member-database signature (Pfam, CDD, G3DSA, PANTHER, SUPERFAMILY) is not a
# curated unit at all - InterPro types many of them as homologous superfamilies
# in their own right, so keeping them while dropping the InterPro entries that
# say the same thing would filter by who issued the accession rather than by
# what the entry is.
_PRIMARY_ENTRY_TYPES = frozenset({'Domain', 'Repeat'})

# Two entries with the same accession are one physical domain, and the shorter
# dropped, only when they overlap by at least this fraction of the shorter one.
# A majority, not any overlap: two tandem copies of a repeat sharing a boundary
# residue are two domains.
_SAME_ID_OVERLAP = 0.5


# One admitted class, and two labels for what is not admitted, kept apart because
# they are not-a-domain for different reasons: TIER_MEMBER is an accession
# InterPro did not issue, TIER_IGNORED an InterPro entry of a type that is not a
# structural unit. Both are dropped by filter_representative_domains(); the
# labels exist so a caller looking at an unfiltered frame - the PDF - can say why
# an entry is not there.
TIER_PRIMARY, TIER_MEMBER, TIER_IGNORED = '1', '2', '-'


def domain_entry_tiers(df_domains):
    """Each row's tier on the ladder filter_representative_domains() applies, as a
    Series of TIER_* labels indexed like `df_domains`.

      TIER_PRIMARY : an InterPro Domain or Repeat entry - the only kind kept
      TIER_MEMBER  : a member-database hit (Pfam, CDD, G3DSA, PANTHER, SUPERFAMILY,
                     ...), i.e. any accession InterPro did not issue - dropped
      TIER_IGNORED : every other InterPro entry - Family, Homologous_superfamily,
                     the residue features, and an IPR of unknown type - dropped

    None when the frame carries no domain_id/type columns at all - the
    DomainEvent/DomainType tables, which have no entry type to rank by - the same
    condition under which the filter returns its input untouched. A
    RepresentativeDomains frame whose rows are ALL untyped is not that case: every
    one of them is a member-DB signature, and they are ranked (and dropped) as
    such rather than waved through.

    Shared with the PDF, so a domain is labelled with the tier the analysis judged
    it on rather than one re-derived alongside it.
    """
    if df_domains is None or len(df_domains) == 0:
        return None
    if 'domain_id' not in df_domains.columns or 'type' not in df_domains.columns:
        return None

    is_ipr = df_domains['domain_id'].astype(str).str.startswith('IPR')
    etype = df_domains['type']
    tiers = pd.Series(TIER_MEMBER, index=df_domains.index)   # non-IPR by default
    tiers[is_ipr] = TIER_IGNORED                             # IPR, incl. unknown type
    tiers[is_ipr & etype.isin(_PRIMARY_ENTRY_TYPES)] = TIER_PRIMARY
    return tiers


def _aa_overlap(s1, e1, s2, e2):
    """True if [s1,e1] and [s2,e2] overlap by at least one residue."""
    return s1 <= e2 and s2 <= e1


def _aa_overlap_fraction(s1, e1, s2, e2):
    """Residues shared by [s1,e1] and [s2,e2] as a fraction of the SHORTER of the
    two - so a short entry sitting inside a long one scores 1.0, not the small
    fraction of the long one it happens to cover. Coordinates are inclusive.
    0.0 when they don't overlap."""
    if not _aa_overlap(s1, e1, s2, e2):
        return 0.0
    overlap = min(e1, e2) - max(s1, s2) + 1
    shorter_length = min(e1 - s1 + 1, e2 - s2 + 1)
    if shorter_length <= 0:
        return 0.0
    return overlap / shorter_length


def filter_representative_domains(df_domains):
    """
    Reduce a single transcript's representative domains to a clean domain set,
    ranked by the curated InterPro entry `type`.

    Only the InterPro Domain and Repeat entries survive - the curated structural
    units. Everything else is dropped outright, and a protein annotated with
    nothing else has no domains, which is the honest answer rather than a gap
    papered over with a weaker assignment:

      - Family and Homologous_superfamily say what the protein IS rather than
        delimiting a unit within it.
      - Active_site / Binding_site / Conserved_site / PTM are residue positions
        rather than regions.
      - Member-database hits (Pfam, CDD, G3DSA, PANTHER, SUPERFAMILY, ...) are
        signatures, not curated units, and InterPro types many of them as
        homologous superfamilies itself - G3DSA:3.30.160.60 and
        G3DSA:1.20.140.150 both come back `homologous_superfamily` with
        `integrated: null` from its API. Keeping them while dropping the InterPro
        entries that say the same thing filtered by who issued the accession
        rather than by what the entry is.

    Dropping them costs coverage: they were ~a fifth of the reported rows, and
    proteins whose only annotation is a member-DB hit now have no domains at all.
    That is the intended trade - a reported domain change now always rests on a
    curated InterPro Domain or Repeat.

    Then collapse genuine duplicates: two kept rows with the SAME domain_id whose
    overlap covers at least _SAME_ID_OVERLAP of the shorter one -> keep the longer.
    Same accession at disjoint or barely-touching positions (e.g. two tandem RRM
    instances) is kept as two domains. Return sorted by AA_start.

    Deliberately NOT handled here: cross-transcript identity when canonical and
    compared annotate one physical domain under different accessions, and
    surfacing an event region that carries no InterPro Domain entry at all. Both
    need a sequence-level model rather than the accessions DoChaP stores.

    Requires 'domain_id' and 'type' columns; a frame without them (the
    DomainEvent/DomainType tables) is returned unchanged. A frame whose `type` is
    entirely NULL is NOT waved through - those rows are all member-DB signatures,
    and a protein annotated with nothing else has no domains. A frame of
    one row is NOT short-circuited: whether that row is a domain does not depend on
    what surrounds it, so a lone Family or site entry is dropped like any other.
    """
    if 'domain_id' not in df_domains.columns or 'type' not in df_domains.columns:
        return df_domains

    df = df_domains
    dom_id = df['domain_id'].astype(str)
    starts = df['AA_start']
    ends = df['AA_end']

    # One definition of the tiers, shared with the PDF - see domain_entry_tiers().
    # Only TIER_PRIMARY is a domain; TIER_MEMBER and TIER_IGNORED rows are not
    # ranked against it, they are simply dropped.
    tiers = domain_entry_tiers(df)
    if tiers is None:
        # An empty frame: the columns are there but there is nothing to rank, so
        # `tiers` is None and `None == TIER_PRIMARY` is a plain False, which is
        # not a valid index. Previously masked by the `type.isna().all()` guard -
        # vacuously true on an empty frame - which this rule no longer wants.
        return df_domains
    keep = df.index[tiers == TIER_PRIMARY].tolist()   # InterPro Domain/Repeat

    # collapse genuine duplicates (same accession, overlapping by a majority of the
    # shorter entry) -> keep the longer
    keep.sort(key=lambda i: (starts[i], ends[i]))
    dropped = set()
    for a in range(len(keep)):
        ia = keep[a]
        if ia in dropped:
            continue
        for b in range(a + 1, len(keep)):
            ib = keep[b]
            if ib in dropped or dom_id[ia] != dom_id[ib]:
                continue
            if _aa_overlap_fraction(starts[ia], ends[ia], starts[ib], ends[ib]) >= _SAME_ID_OVERLAP:
                len_a = ends[ia] - starts[ia]
                len_b = ends[ib] - starts[ib]
                dropped.add(ib if len_a >= len_b else ia)
                if ia in dropped:
                    break

    keep = [i for i in keep if i not in dropped]
    return df.loc[keep].sort_values('AA_start')


def _min_skip_none(values):
    present = [v for v in values if v is not None]
    return min(present) if present else None


def _max_skip_none(values):
    present = [v for v in values if v is not None]
    return max(present) if present else None


def find_relevant_domain_windows(transcript_exons, domain_lookup, canonical_transcript_id, transcript_id,
                                  canonical_junctions, transcript_junctions, junctions):
    """
    Determine the genomic window around the differing junctions and return the
    domains of the canonical and compared transcript that fall within it.

    Two rounds, not an iteration to convergence: round 1 collects the domains
    overlapping the boundary exons of the event span - the exons of BOTH
    transcripts pooled into one genomic span, so neither is windowed over a
    stretch the other is not; round 2 widens that window to the union of the
    event span and those domains' own genomic span - again pooled across both
    transcripts, so each is windowed identically - and re-collects.
    Round 2 is skipped when round 1 found no domains in either transcript.
    Both rounds select from the already-reduced representative domain set (see
    filter_representative_domains()).

    Returns (t_domains_in_region, c_domains_in_region), each with a 'length'
    column added (AA_end - AA_start + 1).
    """
    junction_idxs = canonical_junctions | transcript_junctions
    # min/max over BOTH coordinates, as above: taking min-of-starts and
    # max-of-ends built a truncated - sometimes inverted - window from a pair
    # written end-first.
    min_bp = min(min(junctions[idx]) for idx in junction_idxs)
    max_bp = max(max(junctions[idx]) for idx in junction_idxs)

    t_exons = transcript_exons[transcript_id]
    c_exons = transcript_exons[canonical_transcript_id]

    # Round 1: window spanning the boundary exons of the differing junctions.
    # The two transcripts bound the event with exons of their own, which need not
    # line up - one may splice where the other reads through. Taking each
    # transcript's own pair would window them over different stretches of the
    # gene and compare domain sets drawn from different regions, so the bounding
    # exons of both are pooled into one genomic span and each transcript is
    # windowed by that. The span is then re-snapped to whole exons per transcript,
    # keeping round 1's defining property: it never cuts an exon in half.
    t_first_exon, t_last_exon = find_boundary_exons(t_exons, min_bp, max_bp)
    c_first_exon, c_last_exon = find_boundary_exons(c_exons, min_bp, max_bp)

    # min/max over all four bounding exons, not over a first/last pair: strand is
    # irrelevant here because find_boundary_exons() works in genomic coordinates
    # (its "first" is the genomically-leftmost exon on either strand), but an
    # event span lying wholly inside one transcript's intron makes that
    # transcript's "first" exon fall to the RIGHT of its "last" one, and pairing
    # them would invert the span.
    bounding = (t_first_exon, t_last_exon, c_first_exon, c_last_exon)
    common_first_bp = min(e['genomic_start_tx'] for e in bounding)
    common_last_bp = max(e['genomic_end_tx'] for e in bounding)

    t_min_aa, t_max_aa = get_aa_range(*find_boundary_exons(t_exons, common_first_bp, common_last_bp))
    c_min_aa, c_max_aa = get_aa_range(*find_boundary_exons(c_exons, common_first_bp, common_last_bp))

    # The domain set is reduced by curated InterPro entry type, not by the
    # geometry of the hits.
    # Already reduced by the ladder - see build_filtered_domain_lookup(). The
    # domain set is chosen by curated InterPro entry type, not by the geometry of
    # the hits, so it does not depend on the window computed above.
    df_t_domains = domain_lookup(transcript_id)
    df_c_domains = domain_lookup(canonical_transcript_id)

    t_domains_round1 = _domains_in_aa_range(df_t_domains, t_min_aa, t_max_aa)
    c_domains_round1 = _domains_in_aa_range(df_c_domains, c_min_aa, c_max_aa)

    # Round 2: widen the window to cover the domains found in round 1.
    # Skip refinement if no domains found in round 1
    if t_domains_round1.empty and c_domains_round1.empty:
        t_domains_round2, c_domains_round2 = t_domains_round1, c_domains_round1
    else:
        t_min_bp, t_max_bp = find_bp_range_for_domains(t_exons, t_domains_round1)
        c_min_bp, c_max_bp = find_bp_range_for_domains(c_exons, c_domains_round1)
        common_min_bp = _min_skip_none([min_bp, t_min_bp, c_min_bp])
        common_max_bp = _max_skip_none([max_bp, t_max_bp, c_max_bp])

        if common_min_bp is None or common_max_bp is None:
            t_domains_round2, c_domains_round2 = t_domains_round1, c_domains_round1
        else:
            t_first_exon2, t_last_exon2 = find_boundary_exons(t_exons, common_min_bp, common_max_bp)
            c_first_exon2, c_last_exon2 = find_boundary_exons(c_exons, common_min_bp, common_max_bp)

            t_min_aa2, t_max_aa2 = get_aa_range(t_first_exon2, t_last_exon2, common_min_bp, common_max_bp)
            c_min_aa2, c_max_aa2 = get_aa_range(c_first_exon2, c_last_exon2, common_min_bp, common_max_bp)

            t_domains_round2 = _domains_in_aa_range(df_t_domains, t_min_aa2, t_max_aa2)
            c_domains_round2 = _domains_in_aa_range(df_c_domains, c_min_aa2, c_max_aa2)

    # _domains_in_aa_range() returns a boolean-mask slice: already independent,
    # but carrying pandas' "copy of a slice" marker, which trips
    # SettingWithCopyWarning on the ['length'] assignment below. Only the round-2
    # frames are mutated, so copy here rather than in a helper called 4x per pair.
    t_domains_round2 = t_domains_round2.copy()
    c_domains_round2 = c_domains_round2.copy()
    t_domains_round2['length'] = t_domains_round2['AA_end'] - t_domains_round2['AA_start'] + 1
    c_domains_round2['length'] = c_domains_round2['AA_end'] - c_domains_round2['AA_start'] + 1
    return t_domains_round2, c_domains_round2


# ---------------------------------------------------------------------------
# Phase 3: domain comparison and classification
# ---------------------------------------------------------------------------



def _domain_name_sets(df_domains, name_cols=DOMAIN_NAME_COLUMNS):
    """
    Vectorized replacement for {idx: domain_name_set(row) for idx, row in
    df_domains.iterrows()}. .iterrows() builds a full (mixed-dtype) Series
    per row - one of the most expensive ways to iterate a DataFrame. Pulling
    each column to a numpy array once and indexing by position instead
    avoids that Series construction; called for every compared transcript in
    Phase 3, so this matters at scale.
    """
    if df_domains.empty:
        return {}
    columns = [df_domains[col].to_numpy(dtype=object) for col in name_cols]
    result = {}
    for pos, idx in enumerate(df_domains.index):
        names = set()
        for col_values in columns:
            val = col_values[pos]
            if val is None or pd.isna(val):
                continue
            for name in str(val).split(';'):
                name = name.strip()
                if name and name not in ('None', 'nan'):
                    names.add(name)
        result[idx] = names
    return result


def group_by_shared_names(items_with_names):
    """
    items_with_names: iterable of (key, set_of_names).
    Group keys together if their name-sets overlap (transitively).
    Returns: list of lists of keys.
    """
    groups = []  # list of [keys, names]
    for key, names in items_with_names:
        matching = [g for g in groups if g[1] & names]
        if not matching:
            groups.append([[key], set(names)])
            continue
        merged = matching[0]
        for other in matching[1:]:
            merged[0].extend(other[0])
            merged[1] |= other[1]
            groups.remove(other)
        merged[0].append(key)
        merged[1] |= names
    return [keys for keys, _ in groups]


def total_covered_length(df, idxs, start_col='AA_start', end_col='AA_end'):
    """
    Total length covered by the union of [start, end] intervals (inclusive),
    merging overlapping intervals so overlaps aren't double-counted.
    """
    if not idxs:
        return None
    # .at[] is a much cheaper scalar accessor than .loc[] (skips the general
    # multi-axis alignment machinery) - this runs once per matched domain
    # group in Phase 3, so the per-call savings add up over a large file.
    intervals = sorted((df.at[i, start_col], df.at[i, end_col]) for i in idxs)
    total = 0
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= cur_end:  # overlapping interval -> merge
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start + 1
            cur_start, cur_end = start, end
    return total + cur_end - cur_start + 1


def classify_length_pair(t_length, c_length):
    if t_length == c_length:
        return 'unchanged'
    return 'longer' if t_length > c_length else 'shorter'


def classify_domain_change(c_count, t_count, c_length, t_length):
    """
    Classify a group of matched domains (c_count in canonical, t_count in the
    compared transcript) into one of the Phase 3 result categories.
    """
    if c_count == 0:
        return 'added_domain'
    if t_count == 0:
        return 'dropped_domain'
    if c_count == 1 and t_count == 1:
        return classify_length_pair(t_length, c_length)
    if c_count == 1:
        return 'split_domain'
    if t_count == 1:
        return 'merged_domain'
    if c_count == t_count:
        return classify_length_pair(t_length, c_length) + '_domains'
    return 'increased_domain_number' if c_count < t_count else 'reduced_domain_number'


def choose_domain_display_name(names, prefixes=DOMAIN_NAME_PREFIX_PRIORITY):
    # Sort for a deterministic choice: iteration order over a `set` of names
    # depends on Python's per-process string hash seed, which would otherwise
    # make the chosen name vary between runs.
    sorted_names = sorted(names)
    for prefix in prefixes:
        for name in sorted_names:
            if name.lower().startswith(prefix.lower()):
                return name
    return sorted_names[0] if sorted_names else None


def _group_text(c_domains, c_idxs, t_domains, t_idxs, column):
    """One identity group's text from `column`, taken from whichever source supplied
    the domains: RepresentativeDomains under representative domains, DomainType
    otherwise - both reach here under the same column names, `short_description` for
    the entry's name and `description` for its prose.

    A group can hold several entries (a repeat present twice, a domain split in
    two), and a dropped or new domain has entries on one side only. Canonical is
    read first so the text describes the reference where there is one; distinct
    values are joined rather than picked between.
    """
    values = []
    for df, idxs in ((c_domains, c_idxs), (t_domains, t_idxs)):
        if column not in df.columns:
            continue
        for i in idxs:
            value = df.at[i, column]
            if pd.isna(value):
                continue
            text = str(value).strip()
            if text and text.lower() not in ('nan', 'none') and text not in values:
                values.append(text)
    return '; '.join(values) if values else None


def compare_domains(domain_lookup, transcript_exons, canonical_transcript_id, transcript_id,
                     canonical_junctions, transcript_junctions, junctions):
    """
    Compare the domains of `transcript_id` against `canonical_transcript_id`
    within the Phase 2 window, and classify each group of matched/unmatched
    domains.

    Domains are grouped into "identity groups": two domains belong to the same
    group if they share at least one non-empty name (pfam/interpro/smart/etc),
    transitively.

    Yields one event dict per group, classified as:
    - C=0, T>0          -> 'added_domain'
    - C>0, T=0          -> 'dropped_domain'
    - C=1, T=1          -> 'unchanged' / 'longer' / 'shorter' (by length)
    - C=1, T>1          -> 'split_domain'
    - C>1, T=1          -> 'merged_domain'
    - C>1, T>1, C==T    -> 'unchanged_domains' / 'longer_domains' / 'shorter_domains' (by total length)
    - C>1, T>1, C<T     -> 'increased_domain_number'
    - C>1, T>1, C>T     -> 'reduced_domain_number'
    """
    t_domains, c_domains = find_relevant_domain_windows(
        transcript_exons, domain_lookup, canonical_transcript_id, transcript_id,
        canonical_junctions, transcript_junctions, junctions,
    )

    canonical_names = _domain_name_sets(c_domains)
    transcript_names = _domain_name_sets(t_domains)

    tagged_items = (
        [(('C', idx), names) for idx, names in canonical_names.items()]
        + [(('T', idx), names) for idx, names in transcript_names.items()]
    )

    for members in group_by_shared_names(tagged_items):
        c_idxs = sorted((idx for kind, idx in members if kind == 'C'), key=lambda i: c_domains.loc[i, 'AA_start'])
        t_idxs = sorted((idx for kind, idx in members if kind == 'T'), key=lambda i: t_domains.loc[i, 'AA_start'])

        c_count, t_count = len(c_idxs), len(t_idxs)
        c_length = total_covered_length(c_domains, c_idxs)
        t_length = total_covered_length(t_domains, t_idxs)

        names = set()
        for i in c_idxs:
            names |= canonical_names[i]
        for i in t_idxs:
            names |= transcript_names[i]

        yield {
            'event': classify_domain_change(c_count, t_count, c_length, t_length),
            'alternative_transcript_id': transcript_id,
            'domain_id': choose_domain_display_name(names),
            'domain_name': _group_text(c_domains, c_idxs, t_domains, t_idxs, 'short_description'),
            'domain_description': _group_text(c_domains, c_idxs, t_domains, t_idxs, 'description'),
            'canonical_domain_length': c_length,
            'alternative_domain_length': t_length,
            'canonical_domains_number': c_count,
            'alternative_domains_number': t_count,
        }


def _assert_specie_matches_database(df_junctions, gene_specie):
    """Abort when the species carried on the junctions contradicts DoChaP.

    Catches a wrong -species for any gene the database holds - unlike the Ensembl
    prefix check, which is blind to GeneID-keyed genes. Genes absent from the
    database say nothing and are left to the gene_not_in_db path.
    """
    if not gene_specie or 'specie' not in df_junctions.columns:
        return

    expected = df_junctions['specie'].map(
        lambda s: utils.SPECIE_DB_NAME.get(s) if isinstance(s, str) else None)
    actual = df_junctions['gene_ensembl_id'].map(gene_specie)
    mismatched = df_junctions[expected.notna() & actual.notna() & (expected != actual)]
    if mismatched.empty:
        return

    found = sorted({utils.SPECIE_FROM_DB_NAME.get(s, s)
                    for s in actual[mismatched.index].dropna().unique()})
    stated = sorted(mismatched['specie'].dropna().unique())
    examples = ', '.join(str(g) for g in mismatched['gene_ensembl_id'].unique()[:3])
    raise ValueError(
        f"Species mismatch: {len(mismatched)} of {len(df_junctions)} rows are stated as "
        f"{'/'.join(stated)} but their genes are {'/'.join(found)} in the database "
        f"(e.g. {examples}). Re-run with the species the data actually came from."
    )


# Stands in for one side of a rank label when that boundary of the feature falls
# on no exon edge of the reference transcript - an alternative splice site, which
# is exactly what an event's junction often is.
UNKNOWN_EXON = '*'

# Marks the reference transcript's final exon, as the internal format's own rank
# column does ('E13Last'): a junction reaching the end of the transcript rather
# than an interior exon that happens to be numbered 13.
LAST_EXON_SUFFIX = 'Last'


def exon_pair_label(df_exons, feature):
    """Name the exons of `df_exons` a feature joins: 'E2_E3', or 'E2_E4' where it
    skips one, or '*_E5' where its first boundary lands on no exon edge. The
    reference's final exon carries a 'Last' suffix - 'E11_E13Last'.

    The two sides are named in the order the feature states them - the exon at
    start_position first, then the one at end_position - which is what the
    internal format's own rank_h/rank_m do. That order is not always the genomic
    one: the internal file writes a minus-strand junction high coordinate first,
    and its label follows suit.

    Each boundary is matched against the exon edge facing the intron, with the
    same 1bp tolerance find_matching_junction_indices() allows: the lower
    coordinate meets an exon's genomic_end_tx, the higher one an exon's
    genomic_start_tx. Strand does not enter into it. A retained intron labels the
    same way - against a transcript that splices it out, its two ends are that
    intron's flanking exon edges.

    The reference is a transcript, not the feature's own: the label says where in
    that transcript's exon numbering the event sits, so a junction absent from it
    still gets named as long as its ends land on its exon edges.

    None when the reference has no exons at all, and '*_*' when neither end lands
    on one.
    """
    if df_exons is None or df_exons.empty:
        return None

    exon_orders = df_exons['order_in_transcript'].to_numpy()
    last_order = exon_orders.max()
    # The exon below the intron ends where the intron starts; the one above it
    # starts where the intron ends.
    lower_edges = df_exons['genomic_end_tx'].to_numpy()
    upper_edges = df_exons['genomic_start_tx'].to_numpy()

    start_position, end_position = feature
    low = min(start_position, end_position)

    def _side(bp):
        edges = lower_edges if bp == low else upper_edges
        orders = exon_orders[np.abs(edges - bp) <= 1]
        # Exactly one exon, as in the matcher: a boundary two exons share names
        # neither of them.
        if len(orders) != 1:
            return UNKNOWN_EXON
        suffix = LAST_EXON_SUFFIX if orders[0] == last_order else ''
        return f'E{int(orders[0])}{suffix}'

    return f'{_side(start_position)}_{_side(end_position)}'


# Values of the canonical_junction_in_cds / alternative_junction_in_cds columns - where a group's
# junctions sit relative to a transcript's coding sequence.
CDS_IN = 'yes'            # every junction lies wholly inside the CDS
CDS_PARTIAL = 'partial'   # a junction straddles a CDS boundary, or the group is mixed
CDS_OUT = 'no'            # no junction touches the CDS - all of them in a UTR
CDS_NONE = 'no_cds'       # the transcript has no annotated protein, so no CDS at all


def _cds_spans_by_transcript(df_gene_transcripts, transcript_ids, coding_by_transcript):
    """{transcript_id: (low_bp, high_bp)} - each transcript's coding sequence in
    genomic coordinates. None when the frame carries no CDS columns at all, which
    a hand-built one need not.

    A transcript with no annotated protein is left out rather than mapped to a
    span: DoChaP fills cds_start/cds_end with the transcript's own bounds for
    those, so keeping them would read every non-coding transcript as coding end
    to end.
    """
    if 'cds_start' not in df_gene_transcripts.columns or 'cds_end' not in df_gene_transcripts.columns:
        return None

    starts = pd.to_numeric(df_gene_transcripts['cds_start'], errors='coerce')
    ends = pd.to_numeric(df_gene_transcripts['cds_end'], errors='coerce')

    spans = {}
    for transcript_id, start, end in zip(transcript_ids, starts, ends):
        if not coding_by_transcript.get(transcript_id, True):
            continue
        if pd.isna(start) or pd.isna(end):
            continue
        spans[transcript_id] = (min(int(start), int(end)), max(int(start), int(end)))
    return spans


def _is_missing_gene_id(gene_id):
    """True when an event names no gene at all. A reader may leave this as None,
    NaN, or a blank/placeholder string, so all three are treated alike."""
    if gene_id is None:
        return True
    if not isinstance(gene_id, str) and pd.isna(gene_id):
        return True
    return str(gene_id).strip().lower() in ('', 'nan', 'none', 'na', '.')


class ClusterAnalysisResult:
    def __init__(self, cluster_name, gene_ensembl_id, gene_symbol, chromosome=None, as_event_type=None, specie=None, strand=None):
        self.cluster_name = cluster_name
        self.gene_ensembl_id = gene_ensembl_id
        self.gene_symbol = gene_symbol
        self.chromosome = chromosome
        self.as_event_type = as_event_type
        self.specie = specie
        self.strand = strand
        self.canonical_transcript_id = None
        self.junctions = []
        # Per-feature type, parallel to self.junctions (see FEATURE_JUNCTION /
        # FEATURE_RETAINED_INTRON). None means "every feature is a junction",
        # which is what a junctions frame without the column yields.
        self.feature_types = None
        # Whether to work out the three optional columns - rank,
        # canonical_junction_in_cds and alternative_junction_in_cds. Off by default: they are
        # computed per group and per compared transcript, and the writer leaves
        # them out of the CSV entirely unless the run asked for them.
        self.extra_columns = False
        # {transcript_id: (low_bp, high_bp)} genomic CDS spans, filled from the
        # gene's Transcripts rows - see _cds_spans_by_transcript(). None until
        # then, and afterwards too when the frame names no CDS columns.
        self.cds_spans = None
        # The canonical transcript's exons, once it is resolved. The reference
        # every rank label and every logged junction is named against.
        self.canonical_exons = None
        # What the analysis worked out and the drawing should not work out again:
        # {transcript_id: [junction indices it carries]} and {transcript_id:
        # domains left by the ladder}. Filled during analyze(); read by the PDF
        # (see JunctionsAnalysis._generate_pdfs).
        self.matched_features = {}
        self.kept_domains = {}
        self.events = []

    def add_event(self, event, alternative_transcript_id=None, domain_id=None, domain_name=None, domain_description=None,
                  canonical_domain_length=None, alternative_domain_length=None,
                  canonical_domains_number=None, alternative_domains_number=None, is_longest_cds=None,
                  is_most_like_canonical=None, group=None, rank=None,
                  canonical_junction_in_cds=None, alternative_junction_in_cds=None):
        self.events.append((event, alternative_transcript_id, group, rank, domain_id, domain_name, domain_description, canonical_domain_length,
                            alternative_domain_length, canonical_domains_number, alternative_domains_number,
                            canonical_junction_in_cds, alternative_junction_in_cds,
                            is_longest_cds, is_most_like_canonical))

    def analyze(self, df_gene_transcripts, canonical_transcript_ids, exon_lookup, domain_lookup,
                canonical_rank=None, write_all_comparable=False):
        """
        Run the DOMAS algorithm for this cluster:

        Phase 1 - find the canonical transcript and, for every other transcript
        of the gene, the junctions (if any) that match its exon structure and
        the subset of those that are unique compared to the canonical transcript.
        Where DoChaP flags no canonical transcript for the gene, the longest-CDS
        transcript stands in for it (protein-coding candidates first), so the
        cluster is still analyzed rather than dropped.

        Phase 1.5 - split the comparable transcripts into the cluster's distinct
        events: transcripts adding the same set of features to the canonical form one
        group, and a group whose set is a subset of another's is dropped as the
        lesser account of the same region (see _group_by_unique_features()). Each
        surviving group then gets its own representative, tagged with whether the
        longest-CDS rule and/or the most-like-canonical rule (see
        select_most_like_canonical()) picked it. Both rules run over the group's
        protein-coding candidates where there are any (step 1 of the priority),
        falling through to all of them where there are none. Exactly one transcript
        per group is tagged is_longest_cds; is_most_like_canonical is left unset
        across a group when none of its transcripts qualifies.

        "Longest CDS" means the coding length in bases, not the genomic span from
        cds_start to cds_end - see _load_exons_and_cds_lengths().

        Phase 2/3 - for each transcript with a unique junction, determine the
        relevant genomic window and compare its domains against the canonical
        transcript's, recording one event per domain group.

        Each step below records its own outcome event and reports whether the
        cluster can go on; the steps run in order and any of the first three can
        end the analysis.
        """
        resolved = self._resolve_gene_transcripts(df_gene_transcripts)
        if resolved is None:
            return
        gene_transcript_ids, coding_by_transcript = resolved

        transcript_exons, cds_length_by_transcript = self._load_exons_and_cds_lengths(
            gene_transcript_ids, exon_lookup)

        if not self._resolve_canonical(gene_transcript_ids, canonical_transcript_ids, canonical_rank,
                                       coding_by_transcript, cds_length_by_transcript):
            return

        self.canonical_exons = transcript_exons.get(self.canonical_transcript_id)

        transcript_junctions, canonical_junctions = self._match_features_to_transcripts(transcript_exons)
        # Recorded whatever happens next: which junctions a transcript carries is
        # worth drawing even for one that never gets compared.
        self.matched_features = {
            tid: [self.junctions[i] for i in sorted(idxs) if i < len(self.junctions)]
            for tid, idxs in transcript_junctions.items()
        }
        if canonical_junctions is None:
            return

        # One pass of the ladder per transcript, shared with the drawing.
        domain_lookup, self.kept_domains = build_filtered_domain_lookup(domain_lookup)

        unique_by_transcript = self._find_comparable_transcripts(
            transcript_junctions, canonical_junctions)
        if not unique_by_transcript:
            # An event-level outcome, alongside the per-transcript reasons
            # _find_comparable_transcripts() has already recorded: every transcript
            # of the gene either carries no feature of the event or carries only
            # features the canonical one has too, so there is nothing to compare
            # the canonical against and the cluster yields no group at all.
            self.add_event('no_unique_transcript')
            logger.debug(
                f"No transcript with a unique junction for cluster {self.cluster_name}, "
                f"specie {self.specie}. Skipping analysis."
            )
            return

        # One comparison per distinct event in the cluster, not one per cluster.
        for group_index, group_features, group_transcript_ids in self._group_by_unique_features(unique_by_transcript):
            longest_cds_transcript_id, most_like_canonical_transcript_id = self._select_representatives(
                group_transcript_ids, coding_by_transcript, cds_length_by_transcript, transcript_exons)

            # Same priority as selected_comparable_rows() and
            # results_stats.select_representative_transcript(): most-like-canonical
            # where one qualifies, else longest-CDS. Applied here rather than at write
            # time so the domains of the transcripts that would be discarded are never
            # fetched or compared.
            compared = group_transcript_ids
            if not write_all_comparable:
                selected = most_like_canonical_transcript_id or longest_cds_transcript_id
                if selected is not None:
                    compared = [selected]

            rank_label = self._rank_label(group_features)
            self._log_group(group_index, compared)
            self._record_not_chosen(group_index, group_features, group_transcript_ids,
                                    compared, rank_label,
                                    most_like_canonical_transcript_id is not None)

            self._compare_transcripts(
                compared, transcript_junctions, canonical_junctions,
                transcript_exons, domain_lookup,
                longest_cds_transcript_id, most_like_canonical_transcript_id,
                group_index, rank_label, group_features)

    def _resolve_gene_transcripts(self, df_gene_transcripts):
        """The gene's usable transcript ids and which of them are protein-coding,
        or None when the cluster cannot be analysed at all (the reason is recorded
        as the cluster's event)."""
        # No gene named at all, as opposed to one named and not found: LeafCutter
        # clusters are built annotation-free, so overlapping nothing annotated is
        # an expected outcome, not a failed lookup. A missing id alone does not
        # say which it is - a named gene whose symbol did not resolve also arrives
        # without one - so the symbol decides.
        if _is_missing_gene_id(self.gene_ensembl_id) and _is_missing_gene_id(self.gene_symbol):
            self.add_event('no_gene_specified')
            logger.debug(f"No gene named for cluster {self.cluster_name}, specie {self.specie}. Skipping analysis.")
            return None

        # Check if gene exists in the database at all. Done before any column is
        # read, so an absent gene can be signalled with a plain empty frame (or
        # None) rather than one carrying the DB's columns.
        if df_gene_transcripts is None or df_gene_transcripts.empty:
            self.add_event('gene_not_in_db')
            logger.debug(f"Gene {self.gene_ensembl_id} ({self.gene_symbol}) not found in database for cluster {self.cluster_name}, specie {self.specie}. Skipping analysis.")
            return None

        # Use an order-preserving dedup (not `set`) so the order in which
        # transcripts are processed - and therefore the order of the output
        # rows - doesn't depend on Python's per-process string hash seed.
        # Invalid placeholder ids (e.g. NaN for transcripts with neither an
        # ensembl nor a refseq id) are dropped so they can't spuriously match
        # another gene's similarly-invalid "canonical" id.
        invalid_ids = {'', 'nan', 'None'}
        combined_ids = df_gene_transcripts.transcript_ensembl_id.fillna(df_gene_transcripts.transcript_refseq_id)
        gene_transcript_ids = [
            tid for tid in dict.fromkeys(combined_ids)
            if tid is not None and not pd.isna(tid) and tid not in invalid_ids
        ]

        # Transcripts carrying an annotated protein. The v1 priority puts
        # protein-coding first, ahead of most-like-canonical and longest CDS: a
        # transcript with no protein has no domains, so comparing it to the
        # canonical one trivially "drops" every domain. DoChaP populates
        # cds_start/cds_end for non-coding transcripts too, so CDS length cannot
        # stand in - protein-id presence is the signal. Absent columns mean
        # "unknown", not "non-coding", so a hand-built frame keeps its candidates.
        protein_columns = [c for c in ('protein_ensembl_id', 'protein_refseq_id')
                           if c in df_gene_transcripts.columns]
        if protein_columns:
            has_protein = pd.Series(False, index=df_gene_transcripts.index)
            for column in protein_columns:
                values = df_gene_transcripts[column]
                has_protein |= values.notna() & ~values.astype(str).str.strip().isin(
                    ['', 'nan', 'None'])
            coding_by_transcript = dict(zip(combined_ids, has_protein))
        else:
            coding_by_transcript = {tid: True for tid in combined_ids}

        # Where each transcript's coding sequence starts and ends, for the
        # canonical_junction_in_cds / alternative_junction_in_cds columns. Read here because this is
        # where the frame's protein columns have already been resolved.
        self.cds_spans = _cds_spans_by_transcript(
            df_gene_transcripts, combined_ids, coding_by_transcript)

        if len(gene_transcript_ids) == 1:
            self.add_event('only_one_transcript')
            logger.debug(f"Only one transcript found for cluster {self.cluster_name}, specie {self.specie}.")
            return None

        return gene_transcript_ids, coding_by_transcript

    @staticmethod
    def _load_exons_and_cds_lengths(gene_transcript_ids, exon_lookup):
        """Each transcript's exons, and its coding length in bases.

        The length is the largest CDS-relative exon offset. Interior coding exons
        are wholly coding, so that equals the summed exonic CDS (verified on 4,000
        coding transcripts). NOT cds_end - cds_start, a genomic span that counts
        the introns between the first and last coding exon - a median ~10x the
        coding length, ranking transcripts partly by intron content.
        """
        transcript_exons = {
            transcript_id: exon_lookup(transcript_id)
            for transcript_id in gene_transcript_ids
        }

        cds_length_by_transcript = {}
        for transcript_id, exons in transcript_exons.items():
            length = -1
            if len(exons) and 'abs_end_CDS' in exons.columns:
                largest_offset = pd.to_numeric(exons['abs_end_CDS'], errors='coerce').max()
                if pd.notna(largest_offset):
                    length = largest_offset
            cds_length_by_transcript[transcript_id] = length

        return transcript_exons, cds_length_by_transcript

    def _resolve_canonical(self, gene_transcript_ids, canonical_transcript_ids, canonical_rank,
                           coding_by_transcript, cds_length_by_transcript):
        """Set self.canonical_transcript_id. False when the gene has none and none
        can stand in, the cluster being recorded as no_canonical_transcript."""
        gene_canonical_ids = canonical_transcript_ids.intersection(gene_transcript_ids)
        if not gene_canonical_ids:
            # No transcript flagged canonical - common for genes annotated by
            # RefSeq alone, where neither MANE Select nor RefSeq Select names a
            # representative. The longest-CDS transcript stands in, by the rule
            # used for the comparable transcripts: protein-coding candidates
            # first, then longest coding sequence, then lowest id. A substitute,
            # not an annotation: "canonical" here means DOMAS's choice.
            coding_ids = [tid for tid in gene_transcript_ids if coding_by_transcript.get(tid, True)]
            fallback_candidates = coding_ids or gene_transcript_ids
            if not fallback_candidates:
                self.add_event('no_canonical_transcript')
                logger.debug(f"No canonical transcript found for cluster {self.cluster_name}, specie {self.specie}. Skipping analysis.")
                return False
            self.canonical_transcript_id = select_longest_cds(fallback_candidates, cds_length_by_transcript)
            logger.warning(
                f"No canonical transcript for cluster {self.cluster_name}, specie {self.specie}. "
                f"Using the longest-CDS transcript {self.canonical_transcript_id} "
                f"(CDS {cds_length_by_transcript.get(self.canonical_transcript_id, -1)} bases) instead."
            )
            return True

        # A gene can carry more than one canonical transcript, common outside
        # human: ~4,800 mouse and ~7,100 rat genes have one flagged by RefSeq
        # (canonical=1) and another by Ensembl (2); MANE makes the two agree in
        # human, merging them into one canonical=3 row. Prefer 3, then 2, then 1 -
        # CanonicalEnum ranks in that order - with ties on the lowest id.
        ranked_ids = sorted(gene_canonical_ids)
        if canonical_rank:
            self.canonical_transcript_id = max(
                ranked_ids, key=lambda tid: canonical_rank.get(tid, 0))
        else:
            self.canonical_transcript_id = ranked_ids[0]
        if len(gene_canonical_ids) > 1:
            logger.warning(
                f"Multiple canonical transcripts found for cluster {self.cluster_name}, "
                f"specie {self.specie}: "
                f"{ {tid: canonical_rank.get(tid) for tid in ranked_ids} if canonical_rank else ranked_ids}. "
                f"Using {self.canonical_transcript_id} (highest canonical flag, then lowest id)."
            )
        return True

    def _match_features_to_transcripts(self, transcript_exons):
        """Which of the event's features each transcript carries, and which of them
        the canonical transcript carries. Features matching no transcript at all are
        recorded as feature_not_mapped. The canonical set is None - ending the
        analysis - when the canonical transcript carries none of them.
        """
        transcript_junctions = {
            transcript_id: find_matching_junction_indices(exons, self.junctions, strand=self.strand or '+',
                                                          feature_types=self.feature_types)
            for transcript_id, exons in transcript_exons.items()
        }

        unmapped = 0
        for idx, junction in enumerate(self.junctions):
            if not any(idx in junction_idxs for junction_idxs in transcript_junctions.values()):
                logger.debug(f"Junction {junction} in cluster {self.cluster_name} does not map to any transcript. ")
                self.add_event('feature_not_mapped', None)
                unmapped += 1

        # One feature of an event failing to map is ordinary - a novel junction is
        # exactly what a splicing tool reports. EVERY feature failing is not: it
        # says the coordinates do not fit this gene at all, which is a wrong
        # build, a wrong orientation or a wrong gene rather than novel biology.
        # Warned rather than logged at debug because that is a property of the
        # input worth interrupting for: the whole cluster is about to be dropped.
        if unmapped and unmapped == len(self.junctions) > 1:
            logger.warning(
                "Cluster %s, specie %s: none of its %d features maps to any of %s's "
                "%d transcripts. Check the coordinates against the gene's build and "
                "orientation.",
                self.cluster_name, self.specie, unmapped,
                self.gene_symbol or self.gene_ensembl_id, len(transcript_junctions),
            )

        canonical_junctions = transcript_junctions.get(self.canonical_transcript_id, set())
        if not canonical_junctions:
            self.add_event('no_canonical_features')
            logger.debug(f"No canonical junctions found for cluster {self.cluster_name}, specie {self.specie}. Skipping analysis.")
            return transcript_junctions, None

        return transcript_junctions, canonical_junctions

    def _find_comparable_transcripts(self, transcript_junctions, canonical_junctions):
        """{transcript_id: the features it carries that the canonical one lacks}, for
        the transcripts that carry at least one. The rest are recorded as carrying no
        feature of the event, or none unique to them.

        The unique set is returned, not just the id: which features a transcript adds
        is what separates one event in the cluster from another (see
        _group_by_unique_features()).
        """
        unique_by_transcript = {}
        for transcript_id, junction_idxs in transcript_junctions.items():
            if transcript_id == self.canonical_transcript_id:
                continue
            if not junction_idxs:
                logger.debug(f"Transcript {transcript_id} in cluster {self.cluster_name}, specie {self.specie} does not have any junctions. ")
                self.add_event('transcript_doesnt_have_features', alternative_transcript_id=transcript_id)
                continue

            unique_junctions = junction_idxs - canonical_junctions
            if not unique_junctions:
                logger.debug(
                    f"Transcript {transcript_id} in cluster {self.cluster_name}, specie {self.specie} does not have any unique junctions "
                    "compared to the canonical transcript. Skipping this transcript for comparison."
                )
                self.add_event('no_unique_features', alternative_transcript_id=transcript_id)
                continue

            unique_by_transcript[transcript_id] = frozenset(unique_junctions)
        return unique_by_transcript

    def _group_by_unique_features(self, unique_by_transcript):
        """The cluster's distinct events, as [(group_index, transcript_ids)].

        A LeafCutter cluster is a set of junctions that share a splice site, and it
        routinely holds more than one event: transcripts that add different features
        to the canonical are describing different things and cannot be represented by
        a single comparison. Transcripts are therefore grouped by the exact set of
        features they add, and each surviving group is compared in its own right.

        A group whose feature set is a proper subset of another group's is dropped:
        both describe the same region, and the larger set is the fuller account of
        it, so the smaller one would only report a partial version of the same
        change. Its transcripts are recorded as subsumed_by_larger_event rather than
        dropped silently. Subset is transitive, so one pass over the pairs is enough.

        Groups are numbered from 1 in order of their features, which keeps the index
        stable between runs rather than dependent on dict iteration order.
        """
        if not unique_by_transcript:
            return []

        by_features = {}
        for transcript_id, features in unique_by_transcript.items():
            by_features.setdefault(features, []).append(transcript_id)

        kept = [features for features in by_features
                if not any(features < other for other in by_features)]

        # Numbered before the subsumption pass, not after, so a dropped group can
        # be reported against the number of the group that displaced it.
        groups = [(index, features, by_features[features])
                  for index, features in enumerate(sorted(kept, key=sorted), start=1)]
        group_number = {features: index for index, features, _ in groups}

        # Logged before the subsumption pass so the log defines a group number
        # before anything refers to it.
        for index, features, transcript_ids in groups:
            logger.info(
                "Cluster %s, specie %s: group %d is defined by %s; carried by %s.",
                self.cluster_name, self.specie, index,
                self._features_text(features), ', '.join(sorted(transcript_ids)),
            )

        for features, transcript_ids in by_features.items():
            if features in group_number:
                continue
            # Every dropped set is a proper subset of at least one kept set -
            # subset is transitive, so following the chain up ends at a maximal
            # one - but it can be a subset of several, and naming them all says
            # more than picking one would.
            subsuming = sorted(group_number[other] for other in kept if features < other)
            for transcript_id in sorted(transcript_ids):
                logger.info(
                    "Cluster %s, specie %s: transcript %s adds %s - a subset of group %s, "
                    "which gives the fuller account of the same region. Not compared "
                    "(subsumed_by_larger_event).",
                    self.cluster_name, self.specie, transcript_id,
                    self._features_text(features),
                    ' and '.join(str(number) for number in subsuming) or '?',
                )
                self.add_event('subsumed_by_larger_event', alternative_transcript_id=transcript_id)

        return groups

    def _select_representatives(self, comparable_transcript_ids, coding_by_transcript,
                                cds_length_by_transcript, transcript_exons):
        """(longest_cds, most_like_canonical) - which transcript each rule picks, or
        None where it picks nothing. Under write_all_comparable these only tag the
        rows; otherwise they also decide which single transcript is compared at all,
        so the work for the rest is never done."""
        if not comparable_transcript_ids:
            return None, None

        # Step 1 of the priority: prefer protein-coding candidates. A priority,
        # not a hard filter - where no candidate is coding they all tie here and
        # selection falls through to the structural and length steps, so a
        # cluster still resolves to one transcript rather than none.
        coding_candidates = [tid for tid in comparable_transcript_ids
                             if coding_by_transcript.get(tid, True)]
        selection_candidates = coding_candidates or comparable_transcript_ids

        longest_cds_transcript_id = select_longest_cds(
            selection_candidates, cds_length_by_transcript)
        most_like_canonical_transcript_id = select_most_like_canonical(
            selection_candidates, self.canonical_transcript_id, transcript_exons,
            self.junctions, cds_length_by_transcript,
        )
        return longest_cds_transcript_id, most_like_canonical_transcript_id

    def _features_text(self, features):
        """The junctions of one group, written out for the log: coordinates, each
        with the canonical transcript's exon pair where its ends land on one.

        Computed from the exon table whatever the run asked for - the rank column
        is optional, the log is where the reader reconstructs how a cluster was
        split and should always say which junctions define which group.
        """
        parts = []
        for index in sorted(features):
            if index >= len(self.junctions):
                continue
            start, end = self.junctions[index]
            label = exon_pair_label(self.canonical_exons, self.junctions[index])
            parts.append(f"{start}-{end}" + (f" [{label}]" if label else ""))
        return ', '.join(parts) or 'no junction'

    def _log_group(self, group_index, compared):
        """Which of a group's transcripts the run actually compares. Separate from
        the group's definition, logged in _group_by_unique_features(): the choice
        is only made once the group's representatives have been selected."""
        logger.info(
            "Cluster %s, specie %s: group %d is represented by %s.",
            self.cluster_name, self.specie, group_index, ', '.join(sorted(compared)),
        )

    def _record_not_chosen(self, group_index, features, group_transcript_ids, compared,
                           rank_label, most_like_canonical_exists):
        """Record every transcript of a group that the selection rule passed over.

        Without this they would leave no trace: only one transcript per group is
        compared by default, and the rest simply never reach the output. The rows
        are non-comparisons - the transcript carries the group's junctions and
        could have been compared, it just was not the one picked.
        """
        not_chosen = [tid for tid in group_transcript_ids if tid not in set(compared)]
        if not not_chosen:
            return

        rule = 'most like the canonical' if most_like_canonical_exists else 'longest CDS'
        for transcript_id in sorted(not_chosen):
            logger.info(
                "Cluster %s, specie %s: transcript %s carries group %d (%s) but %s was "
                "picked to represent it (%s). Not compared (transcript_not_chosen).",
                self.cluster_name, self.specie, transcript_id, group_index,
                self._features_text(features), ', '.join(sorted(compared)), rule,
            )
            self.add_event('transcript_not_chosen', alternative_transcript_id=transcript_id,
                           group=group_index, rank=rank_label,
                           canonical_junction_in_cds=self._junctions_in_cds(
                               self.canonical_transcript_id, features),
                           alternative_junction_in_cds=self._junctions_in_cds(
                               transcript_id, features))

    def _rank_label(self, features):
        """The rank labels of the features that define one group, joined.

        Each names the canonical transcript's exons that one junction joins (E2_E3,
        *_E5 - see exon_pair_label()). A row covers one group, so it carries the
        ranks of the junctions that separate that group from the cluster's others,
        not the cluster's whole set. The canonical transcript is the reference for
        the same reason it is everywhere else here: the event is described as a
        departure from it.
        """
        if not self.extra_columns:
            return None
        labels = []
        for index in sorted(features):
            if index >= len(self.junctions):
                continue
            label = exon_pair_label(self.canonical_exons, self.junctions[index])
            if label and label not in labels:
                labels.append(label)
        return '; '.join(labels) if labels else None

    def _junctions_in_cds(self, transcript_id, features):
        """Where the junctions defining one group sit relative to a transcript's
        coding sequence: CDS_IN when every one of them lies wholly between its
        cds_start and cds_end, CDS_OUT when none does (all in a UTR, or outside the
        transcript altogether), CDS_PARTIAL when the group is mixed or a single
        junction straddles a CDS boundary - the start or stop codon falling between
        its two splice sites. CDS_NONE for a transcript with no annotated protein:
        there is no coding region for the event to fall in.

        The test is positional, so it answers for the canonical transcript as well
        as the compared one even though a group's junctions are by definition
        absent from the canonical: what it asks is whether the region the event
        changes is coding in that transcript.

        None - and the column is left out of the CSV - unless the run asked for
        the extra columns, and where the transcripts frame names no CDS at all.
        """
        if not self.extra_columns or self.cds_spans is None:
            return None
        span = self.cds_spans.get(transcript_id)
        if span is None:
            return CDS_NONE

        cds_start, cds_end = span
        labels = set()
        for index in sorted(features):
            if index >= len(self.junctions):
                continue
            start, end = self.junctions[index]
            ends_inside = ((cds_start <= min(start, end) <= cds_end)
                           + (cds_start <= max(start, end) <= cds_end))
            labels.add({2: CDS_IN, 1: CDS_PARTIAL, 0: CDS_OUT}[ends_inside])

        if not labels:
            return None
        return labels.pop() if len(labels) == 1 else CDS_PARTIAL

    def _compare_transcripts(self, comparable_transcript_ids, transcript_junctions,
                             canonical_junctions, transcript_exons, domain_lookup,
                             longest_cds_transcript_id, most_like_canonical_transcript_id,
                             group_index, rank_label=None, group_features=()):
        """Compare each transcript's domains against the canonical transcript's,
        recording one event per domain group - or no_domains_in_region where the
        comparison happened but neither side has a domain there.

        The tags are per group: with several events in a cluster each one has its own
        longest-CDS and most-like-canonical transcript."""
        # Same junctions for both columns - the group's own - measured against a
        # different transcript's CDS. The canonical one is fixed across the group.
        canonical_in_cds = self._junctions_in_cds(self.canonical_transcript_id, group_features)

        for transcript_id in comparable_transcript_ids:
            junction_idxs = transcript_junctions[transcript_id]
            is_longest_cds = transcript_id == longest_cds_transcript_id
            is_most_like_canonical = transcript_id == most_like_canonical_transcript_id
            transcript_in_cds = self._junctions_in_cds(transcript_id, group_features)
            events = list(compare_domains(
                domain_lookup, transcript_exons, self.canonical_transcript_id, transcript_id,
                canonical_junctions, junction_idxs, self.junctions,
            ))
            if events:
                for event in events:
                    self.add_event(**event, is_longest_cds=is_longest_cds,
                                   is_most_like_canonical=is_most_like_canonical,
                                   group=group_index, rank=rank_label,
                                   canonical_junction_in_cds=canonical_in_cds,
                                   alternative_junction_in_cds=transcript_in_cds)
            else:
                self.add_event('no_domains_in_region', alternative_transcript_id=transcript_id,
                                is_longest_cds=is_longest_cds,
                                is_most_like_canonical=is_most_like_canonical,
                                group=group_index, rank=rank_label,
                                canonical_junction_in_cds=canonical_in_cds,
                                alternative_junction_in_cds=transcript_in_cds)


    def get_results_df(self):
        df = pd.DataFrame(
            self.events,
            columns=['event', 'alternative_transcript_id', 'group', 'rank', 'domain_id', 'domain_name', 'domain_description',
                        'canonical_domain_length', 'alternative_domain_length',
                        'canonical_domains_number', 'alternative_domains_number',
                        'canonical_junction_in_cds', 'alternative_junction_in_cds',
                        'is_longest_cds', 'is_most_like_canonical']
        )
        # Nullable integer, not float: the rows that belong to no group (the
        # cluster-level outcomes) leave it empty, and a plain int column holding
        # NaN would print every group index as "1.0".
        df['group'] = df['group'].astype('Int64')
        return df
        

def _analyze_single_cluster(cluster_tuple, exon_lookup=None, domain_lookup=None, canonical_transcript_ids=None,
                           gene_strand=None, transcripts_by_gene=None, canonical_rank=None,
                           write_all_comparable=False, extra_columns=False):
    """Analyze a single cluster."""
    _, cluster_df = cluster_tuple

    gene_ensembl_id = cluster_df.gene_ensembl_id.iat[0]
    gene_symbol = cluster_df.gene_symbol.iat[0]
    event_type = cluster_df.event_type.iat[0]
    specie = cluster_df.specie.iat[0]
    strand = gene_strand.get(gene_ensembl_id)

    cluster_result = ClusterAnalysisResult(
        cluster_df.cluster_name.iat[0], gene_ensembl_id, gene_symbol,
        as_event_type=event_type, specie=specie, strand=strand
    )
    cluster_result.junctions = list(zip(cluster_df['start_position'], cluster_df['end_position']))
    cluster_result.feature_types = [
        FEATURE_JUNCTION if pd.isna(value) else str(value)
        for value in cluster_df[FEATURE_TYPE_COLUMN]
    ]
    cluster_result.extra_columns = extra_columns

    df_gene_transcripts = transcripts_by_gene.get(gene_ensembl_id)
    cluster_result.analyze(df_gene_transcripts, canonical_transcript_ids, exon_lookup, domain_lookup,
                           canonical_rank=canonical_rank,
                           write_all_comparable=write_all_comparable)

    return cluster_result


# Events recorded for a transcript/cluster that is NOT a real comparison to the
# canonical transcript: cluster-level non-comparisons (gene/canonical/junction
# problems) plus per-transcript skips (no junctions / no unique junction). Rows
# carrying these events are the "non-comparable" / "not chosen" transcripts that
# analyze_junctions(filter_non_comparable=True) drops from the output CSV.
NON_COMPARISON_EVENTS = frozenset({
    'no_gene_specified', 'gene_not_in_db', 'no_canonical_transcript', 'only_one_transcript',
    'no_canonical_features', 'feature_not_mapped', 'no_unique_transcript',
    'transcript_doesnt_have_features', 'no_unique_features', 'subsumed_by_larger_event',
    # Carries its group's junctions and could have been compared, but the
    # selection rule picked another transcript of the same group.
    'transcript_not_chosen',
})


def selected_comparable_rows(df_cluster_results):
    """One cluster's result rows with the COMPARISON rows of each group reduced to
    the single transcript the selection rule picks there: the one tagged
    is_most_like_canonical where any is, otherwise the one tagged is_longest_cds.
    Same priority as results_stats.select_representative_transcript(), applied
    before writing. A cluster holding several distinct events keeps one transcript
    per event, not one overall.

    Non-comparison rows (no junctions, no unique junction, a cluster-level outcome)
    are left alone - which of those to write is filter_non_comparable's decision. So
    is a cluster with no tagged row, i.e. no comparable transcript.
    """
    if df_cluster_results.empty:
        return df_cluster_results
    is_comparison = ~df_cluster_results['event'].isin(NON_COMPARISON_EVENTS)
    if not is_comparison.any():
        return df_cluster_results

    comparisons = df_cluster_results[is_comparison]
    # Per group, not per cluster: a cluster holding several distinct events keeps
    # one transcript for each of them.
    keep = pd.Series(False, index=df_cluster_results.index)
    for group, rows in comparisons.groupby('group', dropna=False):
        for column in ('is_most_like_canonical', 'is_longest_cds'):
            tagged = rows.loc[rows[column] == True, 'alternative_transcript_id'].dropna()  # noqa: E712
            if len(tagged):
                keep |= is_comparison & (df_cluster_results['group'] == group) \
                        & (df_cluster_results['alternative_transcript_id'] == tagged.iat[0])
                break
        else:
            keep |= is_comparison & (df_cluster_results['group'] == group)
    return df_cluster_results[~is_comparison | keep]


def _csv_writer_worker(result_queue, output_path, df_results_columns, logger_instance=None,
                       filter_non_comparable=False, write_all_comparable=False):
    """
    Dedicated writer thread that processes results from a queue and writes to CSV.

    Runs continuously until it receives a None sentinel value.
    Writes results incrementally as they arrive from compute workers.

    filter_non_comparable: if True, rows whose event is in NON_COMPARISON_EVENTS
    (transcripts that were not actually compared to canonical) are dropped.
    write_all_comparable: if False (the default), the comparison rows of each
    cluster are reduced to the selected transcript - see selected_comparable_rows().
    """
    log = logger_instance or logger
    output_dir = tempfile.mkdtemp(prefix='domas_csv_')

    try:
        chunk_num = 0
        total_rows = 0

        while True:
            # Get results from queue (blocks until available)
            chunk_results = result_queue.get()

            # None is sentinel for end-of-stream
            if chunk_results is None:
                break

            # Convert results to DataFrame and write
            if chunk_results:
                result_frames = []
                for cluster_result in chunk_results:
                    df_cluster_results = cluster_result.get_results_df()
                    if df_cluster_results.empty:
                        continue

                    if not write_all_comparable:
                        df_cluster_results = selected_comparable_rows(df_cluster_results)

                    if filter_non_comparable:
                        df_cluster_results = df_cluster_results[
                            ~df_cluster_results['event'].isin(NON_COMPARISON_EVENTS)]
                        if df_cluster_results.empty:
                            continue

                    df_transformed = (
                        df_cluster_results
                        # The per-cluster frame's own 'event' column holds the
                        # outcome label and is renamed first, so the assign below
                        # is free to write the cluster id under that name.
                        .rename(columns={'event': 'event_type'})
                        .assign(
                            event=cluster_result.cluster_name,
                            gene_symbol=cluster_result.gene_symbol,
                            canonical_transcript_id=cluster_result.canonical_transcript_id,
                            specie=cluster_result.specie,
                        )
                    )
                    result_frames.append(df_transformed)

                if result_frames:
                    # Filter out empty DataFrames before concatenating to avoid FutureWarning
                    non_empty_frames = [df for df in result_frames if not df.empty]
                    if non_empty_frames:
                        df_chunk = pd.concat(non_empty_frames, ignore_index=True)

                        # Validate and select columns
                        missing_cols = set(df_results_columns) - set(df_chunk.columns)
                        if missing_cols:
                            raise ValueError(f"Missing columns in result: {missing_cols}")
                        df_chunk = df_chunk[df_results_columns]

                        chunk_path = os.path.join(output_dir, f'chunk_{chunk_num:04d}.csv')
                        df_chunk.to_csv(chunk_path, index=False)

                        chunk_rows = len(df_chunk)
                        total_rows += chunk_rows
                        log.info(f"[Writer] Wrote chunk {chunk_num} ({chunk_rows} rows)")

                        chunk_num += 1

        # Combine all chunks into final CSV. Each chunk was already written by
        # pandas with the same header/column order, so this is a plain
        # file-level concatenation (header from the first chunk, data-only
        # from the rest) rather than reading every chunk back into pandas
        # with read_csv, concatenating in memory, and re-serializing with
        # to_csv - which doubles the I/O and peak memory for no benefit on a
        # large results file.
        chunk_files = sorted([os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith('.csv')])

        if chunk_files:
            log.info(f"[Writer] Combining {len(chunk_files)} chunks...")
            with open(output_path, 'w', newline='') as out_f:
                for i, chunk_path in enumerate(chunk_files):
                    with open(chunk_path, 'r', newline='') as in_f:
                        if i > 0:
                            next(in_f)  # skip this chunk's header line
                        shutil.copyfileobj(in_f, out_f)
            log.info(f"[Writer] Final CSV written: {output_path} ({total_rows} rows)")
        else:
            df_all_results = pd.DataFrame(columns=df_results_columns)
            df_all_results.to_csv(output_path, index=False)
            log.info(f"[Writer] Empty results CSV written: {output_path}")

    finally:
        # Drop the temporary chunk files.
        shutil.rmtree(output_dir, ignore_errors=True)


# Per-worker state, populated once by _init_worker(): df_exons/df_domains can be
# large, and submit() args would re-pickle them for every chunk.
_worker_state = {}


def _init_worker(df_exons, df_domains, canonical_transcript_ids, gene_strand, transcripts_by_gene,
                 canonical_rank=None, write_all_comparable=False, extra_columns=False):
    """ProcessPoolExecutor initializer - runs once when each worker process starts."""
    # A spawned worker starts with no logging configuration, so its records would
    # otherwise reach logging's last-resort handler and the console. domas.py
    # names the run's log file in the environment; append to it.
    log_file = os.environ.get('DOMAS_LOG_FILE')
    if log_file:
        logging.basicConfig(
            filename=log_file, filemode='a', level=logging.INFO,
            format="%(asctime)s %(levelname)s %(name)s %(message)s", force=True,
        )
    _worker_state['exon_lookup'] = build_exon_lookup(df_exons)
    _worker_state['domain_lookup'] = build_domain_lookup(df_domains)
    _worker_state['canonical_transcript_ids'] = canonical_transcript_ids
    _worker_state['canonical_rank'] = canonical_rank
    _worker_state['gene_strand'] = gene_strand
    _worker_state['transcripts_by_gene'] = transcripts_by_gene
    _worker_state['write_all_comparable'] = write_all_comparable
    _worker_state['extra_columns'] = extra_columns


def _process_cluster_chunk(chunk_info):
    """Process a chunk of clusters sequentially (worker function for ProcessPoolExecutor).

    Reads shared, read-only state populated once per worker by _init_worker().

    Args:
        chunk_info: tuple of (worker_id, chunk_index, total_chunks, chunk)
    """
    worker_id, chunk_index, total_chunks, chunk = chunk_info
    chunk_size = len(chunk)

    exon_lookup = _worker_state['exon_lookup']
    domain_lookup = _worker_state['domain_lookup']
    canonical_transcript_ids = _worker_state['canonical_transcript_ids']
    canonical_rank = _worker_state.get('canonical_rank')
    gene_strand = _worker_state['gene_strand']
    transcripts_by_gene = _worker_state['transcripts_by_gene']
    write_all_comparable = _worker_state.get('write_all_comparable', False)
    extra_columns = _worker_state.get('extra_columns', False)

    chunk_results = []
    processed_in_chunk = 0

    for cluster_tuple in chunk:
        try:
            result = _analyze_single_cluster(
                cluster_tuple,
                exon_lookup=exon_lookup,
                domain_lookup=domain_lookup,
                canonical_transcript_ids=canonical_transcript_ids,
                gene_strand=gene_strand,
                transcripts_by_gene=transcripts_by_gene,
                canonical_rank=canonical_rank,
                write_all_comparable=write_all_comparable,
                extra_columns=extra_columns,
            )
            chunk_results.append(result)
        except (KeyError, ValueError, AttributeError, TypeError) as e:
            # Handle expected data issues (malformed input, missing fields, etc.)
            cluster_name = cluster_tuple[1].cluster_name.iat[0] if len(cluster_tuple) > 1 else "unknown"
            logger.error(f"[Worker {worker_id}] Error processing cluster {cluster_name}: {type(e).__name__}: {e}")
            # Continue processing other clusters instead of crashing
            continue

        processed_in_chunk += 1

        # Progress logging every 10000 clusters within the chunk
        if processed_in_chunk % 10000 == 0:
            logger.info(f"[Worker {worker_id}] Chunk {chunk_index + 1}/{total_chunks}: processed {processed_in_chunk}/{chunk_size}")

    return chunk_results


class JunctionsAnalysis:
    def __init__(self, con, logger_instance=None, gene_visualization_cls=GeneVisualization):
        self.con = con
        self.logger = logger_instance or logger
        self.gene_visualization_cls = gene_visualization_cls

    def _load_junctions_data(self, df_junctions, specie=None):
        """Validate the junctions DataFrame. Reading junctions from a file (plain CSV,
        hadas-format Excel, IOE, ...) is alternative_splicing.py's responsibility -
        callers pass an already-loaded DataFrame here."""
        df_junctions = df_junctions.copy()

        if 'cluster_name' not in df_junctions.columns and 'cluster' in df_junctions.columns:
            df_junctions = df_junctions.rename(columns={'cluster': 'cluster_name'})

        # One contract for every reader: required columns validated, optional ones
        # filled with their declared default, specie derived where a reader left it
        # blank. Past this point the columns are simply there, so downstream code
        # does not each carry its own "if 'x' in df.columns" default.
        return utils.normalize_junctions_frame(df_junctions, specie=specie)

    def _filter_junctions_by_transcript_count(self, df_junctions, filter_transcript_count):
        """Filter junctions to genes with exactly filter_transcript_count transcripts."""
        if filter_transcript_count <= 0:
            return df_junctions.copy()

        # Push the per-gene transcript count into SQL instead of pulling the
        # entire transcripts table (`select *`) into pandas just to compute a
        # value_counts() over one column.
        df_gene_counts = pd.read_sql_query(
            'select gene_GeneID_id, count(*) as gene_count from transcripts group by gene_GeneID_id',
            self.con,
        )
        genes_with_count = df_gene_counts.loc[
            df_gene_counts['gene_count'] == filter_transcript_count, 'gene_GeneID_id'
        ].tolist()
        return df_junctions[df_junctions['gene_ensembl_id'].isin(genes_with_count)]

    def _load_database_data(self, gene_ids, use_ensembl_only=False):
        """Load genes, transcripts, domains, and exons from database."""
        clause, params = utils.gene_id_clause(gene_ids)
        df_genes = pd.read_sql_query(
            f"SELECT gene_ensembl_id, gene_GeneID_id, specie, strand FROM Genes WHERE {clause}",
            self.con, params=params
        )
        gene_strand = dict(zip(utils.combined_gene_ids(df_genes), df_genes['strand']))
        gene_specie = dict(zip(utils.combined_gene_ids(df_genes), df_genes['specie']))

        df_transcripts = utils.get_genes_df_transcripts(self.con, gene_ids)

        if use_ensembl_only:
            invalid_ids = {'', 'nan', 'None'}
            has_ensembl_id = df_transcripts.transcript_ensembl_id.notna() & \
                ~df_transcripts.transcript_ensembl_id.isin(invalid_ids)
            df_transcripts = df_transcripts[has_ensembl_id]

        # Combine ensembl and refseq IDs, filtering out invalid entries
        invalid_ids = {'', 'nan', 'None'}
        transcript_ids = set(
            df_transcripts.transcript_ensembl_id.fillna(df_transcripts.transcript_refseq_id)
        ) - {None} - invalid_ids

        df_domains = utils.get_domains_db(
            self.con, transcript_ids
        )

        df_exons = utils.get_exons_for_transcripts(self.con, transcript_ids)

        return df_genes, df_transcripts, df_domains, df_exons, gene_strand, gene_specie

    def _prepare_lookup_structures(self, df_transcripts, df_exons, df_domains):
        """Build lookup structures and transcript groupings."""
        invalid_ids = {'', 'nan', 'None'}
        all_transcript_ids = df_transcripts.transcript_ensembl_id.fillna(df_transcripts.transcript_refseq_id)
        valid_mask = all_transcript_ids.notna() & ~all_transcript_ids.isin(invalid_ids)
        is_canonical = valid_mask & (df_transcripts.canonical != 0)
        canonical_transcript_ids = set(all_transcript_ids[is_canonical])
        # id -> DoChaP canonical flag (1 RefSeq / 2 Ensembl / 3 both), so a gene
        # with several canonical transcripts can prefer the strongest. Built here
        # rather than passed as a frame: it crosses a process boundary per worker.
        canonical_rank = dict(zip(all_transcript_ids[is_canonical],
                                  df_transcripts.canonical[is_canonical].astype(int)))

        # Grouped by the combined gene id, so a gene carrying only a GeneID is
        # reachable under the id the junctions frame holds for it. Grouping on
        # gene_ensembl_id alone dropped those genes: groupby() skips NaN keys, so
        # every one of their transcripts vanished and the cluster reported
        # gene_not_in_db.
        transcripts_by_gene = {
            gid: g for gid, g in df_transcripts.groupby(utils.combined_gene_ids(df_transcripts))
        }

        return canonical_transcript_ids, canonical_rank, transcripts_by_gene

    def _prepare_cluster_groups(self, df_junctions):
        """Group junctions into clusters."""
        # specie is always present after normalize_junctions_frame(), so clusters
        # group per species unconditionally. Grouping on cluster_name alone would
        # merge same-named clusters from different species in a multi-species run.
        #
        # dropna=False because the value may still be None: it is derived from the
        # Ensembl gene id, and an input keyed by GeneID or by a non-Ensembl
        # identifier yields nothing to derive from. groupby() discards NaN keys by
        # default, which would silently drop those clusters instead of analysing
        # them - the same trap that lost every gene without a gene_ensembl_id.
        group_columns = ['specie', 'cluster_name']
        cluster_groups = list(df_junctions.groupby(group_columns, dropna=False))
        return cluster_groups

    def _run_parallel_analysis(self, cluster_groups, df_exons, df_domains, canonical_transcript_ids,
                               gene_strand, transcripts_by_gene, num_workers, output_path,
                               filter_non_comparable=False, canonical_rank=None,
                               write_all_comparable=False, extra_columns=False):
        """Execute cluster analysis in parallel with dedicated writer thread."""
        total = len(cluster_groups)
        actual_workers = min(num_workers, total)

        # Cluster cost varies a lot and is source-correlated (grouped by .ioe file), so a
        # few big contiguous chunks can leave some workers idle for hours. Deterministic
        # shuffling plus many small chunks lets idle workers pick up more work instead.
        shuffled_groups = cluster_groups.copy()
        random.Random(42).shuffle(shuffled_groups)
        chunk_size = max(1, min(200, -(-total // (actual_workers * 20))))  # ceil division

        chunks = [shuffled_groups[i:i + chunk_size] for i in range(0, total, chunk_size)]
        total_chunks = len(chunks)

        # as_completed() yields chunks as they finish, in no relation to
        # cluster_groups' order. PDF names embed a sequential count, so the
        # original order is restored below.
        has_specie_column = 'specie' in cluster_groups[0][1].columns if cluster_groups else False

        def _group_identity(cluster_name, specie):
            return (specie, cluster_name) if has_specie_column else (cluster_name,)

        original_order = {
            _group_identity(group_df['cluster_name'].iat[0], group_df['specie'].iat[0] if has_specie_column else None): idx
            for idx, (_, group_df) in enumerate(cluster_groups)
        }

        self.logger.log(utils.PROGRESS,
                        f"Analyzing {total} clusters with {actual_workers} workers")
        self.logger.info(f"Processing {total_chunks} chunks (~{chunk_size} clusters per chunk)")

        # Prepare column definitions for CSV. The two selection flags are written
        # only under write_all_comparable, where several transcripts share a
        # cluster and the reader needs to know which one the rule picked. With one
        # comparison row per cluster they would be True on every row of it.
        df_results_columns = ['event', 'group', 'gene_symbol', 'specie', 'event_type', 'canonical_transcript_id',
                              'alternative_transcript_id', 'domain_id', 'domain_name',
                              'canonical_domain_length', 'alternative_domain_length', 'canonical_domains_number', 'alternative_domains_number']
        # The optional three, off unless the run asked for them: which exons the
        # group's junctions join, and whether they fall in each side's CDS.
        if extra_columns:
            df_results_columns.insert(2, 'rank')
            df_results_columns += ['canonical_junction_in_cds', 'alternative_junction_in_cds']
        if write_all_comparable:
            df_results_columns += ['is_longest_cds', 'is_most_like_canonical']
        # Last: a paragraph of InterPro prose per domain, which would otherwise push
        # every column a reader scans for off the right of the screen.
        df_results_columns += ['domain_description']

        # Create queue for results (backpressure if writer lags)
        result_queue = queue.Queue(maxsize=actual_workers * 2)

        # Start dedicated writer thread
        writer_thread = threading.Thread(
            target=_csv_writer_worker,
            args=(result_queue, output_path, df_results_columns, self.logger,
                  filter_non_comparable, write_all_comparable),
            daemon=False
        )
        writer_thread.start()

        chunks_with_info = [
            (chunk_idx % actual_workers, chunk_idx, total_chunks, chunk)
            for chunk_idx, chunk in enumerate(chunks)
        ]

        all_results = []  # Collect results for PDF generation
        processed_count = 0
        last_time = time.perf_counter()

        with ProcessPoolExecutor(
            max_workers=actual_workers,
            initializer=_init_worker,
            initargs=(df_exons, df_domains, canonical_transcript_ids, gene_strand, transcripts_by_gene,
                      canonical_rank, write_all_comparable, extra_columns),
        ) as executor:
            # Submit all tasks - df_exons/df_domains are sent once per worker via
            # initargs above, not re-pickled per chunk here.
            futures = [executor.submit(_process_cluster_chunk, chunk_info) for chunk_info in chunks_with_info]

            # Process results as they complete: send to writer AND collect for PDFs
            for future in as_completed(futures):
                chunk_results = future.result()
                result_queue.put(chunk_results)  # Send to writer thread for CSV
                all_results.extend(chunk_results)  # Collect for PDF generation
                processed_count += len(chunk_results)

                if processed_count % 10000 == 0:
                    cur_time = time.perf_counter()
                    self.logger.log(
                        utils.PROGRESS,
                        f"Analyzed {processed_count}/{total} clusters. "
                        f"Last 10000 took {cur_time - last_time:.2f}s. "
                        f"ETA: {(cur_time - last_time) * (total - processed_count) / 10000 / 60:.1f}min"
                    )
                    last_time = cur_time

        # Signal writer to stop
        result_queue.put(None)
        writer_thread.join()  # Wait for writer to finish

        all_results.sort(
            key=lambda r: original_order.get(_group_identity(r.cluster_name, r.specie if has_specie_column else None), total)
        )

        self.logger.log(utils.PROGRESS,
                        f"Analysis complete: {processed_count}/{total} clusters")
        self._log_unmapped_summary(all_results)
        return all_results

    def _log_unmapped_summary(self, results):
        """How much of the input never reached a transcript, per species.

        A per-cluster warning scrolls past in a large run; this is one line, and
        in a two-species comparison the asymmetry between the two is itself the
        diagnosis - one species at 0% against the other at 48.5% is an input
        problem, not biology. That asymmetry existed for 1,684 of 3,484 clusters
        and nothing reported it (see find_matching_junction_indices).
        """
        per_specie = {}
        for result in results:
            features = len(result.junctions)
            if not features:
                continue
            unmapped = sum(1 for event in result.events if event[0] == 'feature_not_mapped')
            counts = per_specie.setdefault(result.specie, [0, 0, 0, 0])
            counts[0] += unmapped
            counts[1] += features
            counts[2] += 1 if unmapped else 0
            counts[3] += 1 if unmapped == features else 0
        if not per_specie:
            return

        parts = []
        for specie in sorted(per_specie, key=lambda s: (s is None, s)):
            unmapped, features, touched, wholly = per_specie[specie]
            parts.append(
                f"{specie or 'unknown'} {unmapped:,} of {features:,} "
                f"({100 * unmapped / features:.1f}%) in {touched:,} clusters, "
                f"{wholly:,} with NO feature mapped"
            )
        self.logger.log(utils.PROGRESS, "Features unmapped: " + "; ".join(parts))

    # Events recorded for a transcript that was NOT actually compared to the
    # canonical transcript (it was skipped for lacking junctions or lacking a
    # unique junction).
    # Transcripts whose domains were never compared to the canonical - not drawn
    # under -show_only_compared. transcript_not_chosen belongs here: it carries
    # its group's junctions and could have been compared, but another transcript
    # of the group represented it, so its domains were never even fetched.
    _SKIPPED_TRANSCRIPT_EVENTS = {
        'transcript_doesnt_have_features', 'no_unique_features', 'subsumed_by_larger_event',
        'transcript_not_chosen',
    }

    def _comparable_transcript_ids(self, cluster_result):
        """Canonical transcript id plus every transcript id that was actually compared to it."""
        comparable_ids = set()
        if cluster_result.canonical_transcript_id:
            comparable_ids.add(cluster_result.canonical_transcript_id)
        df_cluster_results = cluster_result.get_results_df()
        if len(df_cluster_results) > 0:
            compared_mask = ~df_cluster_results['event'].isin(self._SKIPPED_TRANSCRIPT_EVENTS)
            comparable_ids.update(df_cluster_results.loc[compared_mask, 'alternative_transcript_id'].dropna())
        return comparable_ids

    def _generate_pdfs(self, results, print_genes, restrict_to_comparable=False):
        """Generate PDF visualizations for clusters.

        restrict_to_comparable: if True, each PDF draws only the canonical
        transcript and the ones compared to it; transcripts skipped during
        analysis (no junctions, no unique junction) are omitted entirely.
        """
        print_gene_set = {gene.upper() for gene in print_genes} if print_genes is not None else None

        def _gene_symbol_key(value):
            if isinstance(value, str):
                return value.upper()
            return None

        filtered_results = (
            [cluster_result for cluster_result in results if _gene_symbol_key(cluster_result.gene_symbol) in print_gene_set]
            if print_gene_set is not None
            else results
        )
        preloaded_gene_data = prepare_gene_data_bulk(
            self.con, [cluster_result.gene_ensembl_id for cluster_result in filtered_results if _gene_symbol_key(cluster_result.gene_symbol) is not None],
        )

        gene_visualizations = {}
        for count, cluster_result in enumerate(results, start=1):
            try:
                df_cluster_junctions = pd.DataFrame(cluster_result.junctions, columns=['start', 'end'])
                # Carry the feature type into the visualization, so a retained
                # intron is matched there by containment too. Rebuilt from the
                # coordinate pairs alone, the PDF treated every feature as a
                # junction and showed the retaining transcript as not supporting it.
                #
                # Only when the event actually contains one: the column also lands
                # in the first-page junction table, where a "junction" value on
                # every row of every plain event would be noise. Present, it tells
                # the reader which row the dotted in-exon span belongs to.
                if (cluster_result.feature_types is not None
                        and FEATURE_RETAINED_INTRON in cluster_result.feature_types):
                    df_cluster_junctions[FEATURE_TYPE_COLUMN] = cluster_result.feature_types
                cluster_gene_key = _gene_symbol_key(cluster_result.gene_symbol)
                if print_gene_set is not None and cluster_gene_key not in print_gene_set:
                    continue
                viz = gene_visualizations.get(cluster_result.gene_ensembl_id)
                if viz is None:
                    gene_symbol = cluster_result.gene_symbol
                    gene_ensembl_id = cluster_result.gene_ensembl_id
                    preloaded = preloaded_gene_data.get(gene_ensembl_id) if isinstance(gene_ensembl_id, str) else None
                    viz = self.gene_visualization_cls(
                        self.con, gene_symbol, preloaded=preloaded,
                    )
                    gene_visualizations[gene_ensembl_id] = viz

                file_name = f'{cluster_result.gene_symbol}_{count}_junction_comparison.pdf'
                if cluster_result.as_event_type is not None:
                    file_name = f'{cluster_result.as_event_type}_{file_name}'
                if cluster_result.specie is not None:
                    file_name = f'{cluster_result.specie}_{file_name}'

                transcript_ids = self._comparable_transcript_ids(cluster_result) if restrict_to_comparable else None
                no_comparison_note = None
                if restrict_to_comparable and transcript_ids - {cluster_result.canonical_transcript_id} == set():
                    no_comparison_note = "No transcripts available for comparison to the canonical transcript."
                viz.create_pdf(
                    file_name,
                    protein_only=False,
                    domains_only=False,
                    df_junction=df_cluster_junctions,
                    # The description is a CSV column, not a PDF one: it is a
                    # paragraph of InterPro prose per domain, and the PDF draws
                    # this frame as a narrow per-transcript table.
                    df_results=cluster_result.get_results_df().drop(
                        columns=['domain_description', 'domain_name'], errors='ignore'),
                    transcript_ids=transcript_ids,
                    no_comparison_note=no_comparison_note,
                    # Which transcript the comparison used as the reference - not
                    # always the one the DB flags, where a gene carries two.
                    canonical_transcript_id=cluster_result.canonical_transcript_id,
                    # What the analysis decided, so the drawing does not decide it
                    # again: the domains the ladder kept per transcript, and the
                    # junctions each transcript was found to carry.
                    analysis_domains=cluster_result.kept_domains,
                    analysis_features=cluster_result.matched_features,
                )
            except ValueError as e:
                self.logger.warning(f"Warning: Skipping PDF generation for {cluster_result.gene_symbol}, specie {cluster_result.specie}: {e}")

    def analyze_junctions(self, df_junctions, output_path='as_events_junctions_analysis.csv',
                          specie=None, filter_transcript_count=0, create_pdf=True, print_genes=None,
                          num_workers=4, use_ensembl_only=False, restrict_pdf_to_comparable=False,
                          filter_non_comparable=False,
                          write_all_comparable=False, extra_columns=False):
        """
        Analyze junctions and detect domain changes across alternative transcripts.

        Every comparable transcript is compared to canonical (none are skipped); each
        result row is tagged is_longest_cds/is_most_like_canonical for post-hoc filtering
        (see get_results_df()).

        Args:
            df_junctions: DataFrame of junctions. Reading junctions from a file
                (plain CSV, hadas-format Excel, IOE, ...) is alternative_splicing.py's
                responsibility - pass the already-loaded DataFrame here.
            output_path: Path for output CSV file
            filter_transcript_count: If > 0, only analyze genes with exactly this many transcripts
            create_pdf: Whether to generate PDF visualizations
            print_genes: List of gene symbols to generate PDFs for (or all if None)
            num_workers: Number of parallel workers for analysis
            use_ensembl_only: If True, only consider transcripts that have an
                ensembl id - transcripts with no ensembl id (refseq-only) are
                filtered out before any exon/domain lookups are built, so they
                never participate in the analysis at all.
            restrict_pdf_to_comparable: If True, each generated PDF only draws
                the canonical transcript and the transcripts that were actually
                compared to it - every other transcript of the gene is omitted
                from the visualization entirely.
            write_all_comparable: If True, every comparable transcript is compared
                to the canonical one and the CSV keeps a row per transcript, with
                the is_most_like_canonical / is_longest_cds columns naming the one
                the selection rule picks. If False (default), only that transcript
                is compared - so the domains of the others are never fetched - and
                the two columns are omitted, being True on every written row.
            filter_non_comparable: If True, the output CSV contains only rows for
                transcripts that were actually compared to canonical - rows whose
                event is a non-comparison / skip event (see NON_COMPARISON_EVENTS)
                are dropped. The returned ClusterAnalysisResult objects and any
                PDFs are unaffected; only the written CSV is filtered.
            extra_columns: If True, the CSV carries three further columns, for
                every input format: `rank`, naming the canonical transcript's
                exons the group's junctions join (E2_E4, *_E5), and
                canonical_junction_in_cds / alternative_junction_in_cds, saying whether those
                junctions fall inside the canonical and the compared
                transcript's coding sequence. Omitted by default.

        Returns:
            List of ClusterAnalysisResult objects
        """
        # Validate input
        df_junctions = self._load_junctions_data(df_junctions, specie=specie)
        df_junctions = self._filter_junctions_by_transcript_count(df_junctions, filter_transcript_count)

        gene_ids = df_junctions.gene_ensembl_id.unique().tolist()
        self.logger.log(utils.PROGRESS, f"Analyzing {len(gene_ids)} genes")

        # Load data from database
        df_genes, df_transcripts, df_domains, df_exons, gene_strand, gene_specie = self._load_database_data(
            gene_ids, use_ensembl_only=use_ensembl_only
        )

        # The database knows the species of every gene it holds, including those
        # keyed only by GeneID, which the Ensembl-prefix check cannot read. This is
        # therefore the stronger check on the stated species: a wrong -species is
        # caught on any gene present in DoChaP, not only Ensembl-keyed ones.
        _assert_specie_matches_database(df_junctions, gene_specie)

        # Prepare lookup structures
        canonical_transcript_ids, canonical_rank, transcripts_by_gene = \
            self._prepare_lookup_structures(df_transcripts, df_exons, df_domains)

        # Group junctions into clusters
        cluster_groups = self._prepare_cluster_groups(df_junctions)

        # Run parallel analysis (with dedicated writer thread for CSV output)
        results = self._run_parallel_analysis(
            cluster_groups, df_exons, df_domains, canonical_transcript_ids,
            gene_strand, transcripts_by_gene, num_workers,
            output_path, filter_non_comparable=filter_non_comparable,
            canonical_rank=canonical_rank, write_all_comparable=write_all_comparable,
            extra_columns=extra_columns,
        )

        # Generate PDFs if requested
        if create_pdf:
            self._generate_pdfs(
                results, print_genes, restrict_to_comparable=restrict_pdf_to_comparable,
            )

        return results
