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

def find_matching_junction_indices(df_transcript_exons, junctions, strand='+'):
    """
    Return the set of indices (into `junctions`) of junctions that match this
    transcript's exon structure.

    A junction (start_position, end_position) matches the transcript if there
    are two exons, adjacent in transcript order, such that one exon's
    genomic_end_tx is within 1bp of the junction's intron-left boundary and the
    other exon's genomic_start_tx is within 1bp of the junction's intron-right
    boundary.

    On the positive strand the intron-left boundary is start_position and the
    intron-right boundary is end_position.  On the negative strand the
    transcript runs right-to-left in genomic coordinates, so the roles are
    reversed: the intron-left boundary (in genomic terms) is end_position and
    the intron-right boundary is start_position.
    """
    if df_transcript_exons.empty or not junctions:
        return set()

    exon_starts = df_transcript_exons['genomic_start_tx'].to_numpy()
    exon_ends = df_transcript_exons['genomic_end_tx'].to_numpy()
    exon_orders = df_transcript_exons['order_in_transcript'].to_numpy()

    matched = set()
    for idx, (start_position, end_position) in enumerate(junctions):
        if strand == '-':
            # Negative strand: transcript runs from high to low genomic coords.
            # DB always stores genomic_start_tx < genomic_end_tx, so:
            #   - the upstream exon's intron boundary is its genomic_start_tx (lower bound of the higher exon)
            #   - the downstream exon's intron boundary is its genomic_end_tx (upper bound of the lower exon)
            # Junction end_position is near the upstream exon's genomic_start_tx.
            # Junction start_position is near the downstream exon's genomic_end_tx.
            intron_left, intron_right = end_position, start_position
            upstream_orders = exon_orders[np.abs(exon_starts - intron_left) <= 1]
            downstream_orders = exon_orders[np.abs(exon_ends - intron_right) <= 1]
        else:
            intron_left, intron_right = start_position, end_position
            upstream_orders = exon_orders[np.abs(exon_ends - intron_left) <= 1]
            downstream_orders = exon_orders[np.abs(exon_starts - intron_right) <= 1]
        if len(upstream_orders) == 1 and len(downstream_orders) == 1:
            if abs(int(downstream_orders[0]) - int(upstream_orders[0])) == 1:
                matched.add(idx)
    return matched


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
    determinism.

    Returns None when no candidate qualifies - the flag then goes unset rather than
    falling back to the longest-CDS candidate, which is_longest_cds already marks in
    its own right. Callers wanting a single transcript per cluster fall back to that
    flag themselves (see results_stats.select_representative_transcript()).
    """
    min_bp = min(start for start, end in junctions)
    max_bp = max(end for start, end in junctions)

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

    # On a minus strand transcript, the domain's lower CDS bound maps to a
    # HIGHER genomic position (and vice versa), so these two values can come
    # back genomically reversed. Callers pool this pair into a min/max across
    # several transcripts, so it must be returned in genomic (low, high)
    # order regardless of strand - otherwise the reversed pair is silently
    # dropped from the pooled max() and the window never widens to cover it.
    bp_a = _cds_to_bp(first_exon, min_domain_bp)
    bp_b = _cds_to_bp(last_exon, max_domain_bp)
    return min(bp_a, bp_b), max(bp_a, bp_b)


def build_exon_lookup(df_exons):
    """
    Precompute, once per analyze_junctions() run, a transcript_id -> exons
    lookup so per-cluster/per-transcript filtering of the full `df_exons`
    DataFrame (the dominant cost of analyze_junction()) is replaced by O(1)
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


def _merge_domain_names(*values):
    """Merge one or more ';'-separated identifier strings into a deduplicated ';'-separated string."""
    names = []
    for value in values:
        if value is None or pd.isna(value):
            continue
        for name in str(value).split(';'):
            name = name.strip()
            if name and name not in ('None', 'nan') and name not in names:
                names.append(name)
    return '; '.join(names) if names else None


def collapse_contained_domains(df_domains, tolerance=2, overlap_fraction=0.85):
    """
    Collapse domains (of a single transcript) that are the same physical
    domain instance into a single row - the longer domain - merging the
    redundant domains' identifier names into the kept row so
    compare_domains() can match on any of them.

    Two domains are treated as the same instance if EITHER:
    - one contains the other within `tolerance` AA on either side (handles
      near-identical hits whose boundaries differ by only a couple AA), OR
    - their overlap covers at least `overlap_fraction` of the shorter
      domain's length (handles cases like co-occurring cross-database calls,
      e.g. a CDD hit and an InterPro hit for the same fold, whose boundaries
      can differ by more than `tolerance` AA without being a different
      domain - a fixed AA tolerance alone is brittle here since it flips
      on a few AA of alignment drift regardless of domain size).
    """
    if len(df_domains) <= 1:
        return df_domains

    df_domains = df_domains.copy()
    starts = df_domains['AA_start'].to_numpy()
    ends = df_domains['AA_end'].to_numpy()
    index = df_domains.index.to_numpy()
    name_block = df_domains[DOMAIN_NAME_COLUMNS].to_numpy(dtype=object)

    # Process longest domain first, but break ties with a *stable* key so the
    # collapse result never depends on df_domains' own row order. That row order
    # is hash-seed dependent (it flows from a set() of transcript ids into the
    # domain query), and for equal-length overlapping domains the order decides
    # which one is kept vs merged away - so ordering by length alone let a
    # cluster's domain count (and thus its results.csv row count) vary between
    # runs. Ties break by start, then end, then the row's merged identifier set.
    name_key = [
        ";".join(sorted(
            str(v).strip() for v in name_block[pos]
            if not pd.isna(v) and str(v).strip() not in ("", "None", "nan")
        ))
        for pos in range(len(df_domains))
    ]
    order = sorted(
        range(len(df_domains)),
        key=lambda pos: (starts[pos] - ends[pos], starts[pos], ends[pos], name_key[pos]),
    )

    dropped = set()
    merges = {}  # position -> list of positions merged into it
    for oi in order:
        if oi in dropped:
            continue
        for oj in order:
            if oj == oi or oj in dropped:
                continue
            contained = starts[oi] - tolerance <= starts[oj] and ends[oj] <= ends[oi] + tolerance
            if not contained:
                overlap = min(ends[oi], ends[oj]) - max(starts[oi], starts[oj]) + 1
                if overlap > 0:
                    shorter_length = min(ends[oi] - starts[oi] + 1, ends[oj] - starts[oj] + 1)
                    contained = (overlap / shorter_length) >= overlap_fraction
            if contained:
                merges.setdefault(oi, []).append(oj)
                dropped.add(oj)

    if merges:
        for oi, others in merges.items():
            for ci in range(len(DOMAIN_NAME_COLUMNS)):
                values = [name_block[oi, ci]] + [name_block[oj, ci] for oj in others]
                name_block[oi, ci] = _merge_domain_names(*values)
        for ci, col in enumerate(DOMAIN_NAME_COLUMNS):
            df_domains[col] = name_block[:, ci]

    return df_domains.drop(index=index[list(dropped)] if dropped else [])


# InterPro entry types (from RepresentativeDomains.type, sourced from
# interpro.xml.gz). Domain/Repeat are the real structural-functional units;
# Family/Homologous_superfamily are broader groupings; the site/PTM types are
# residue features, not domains.
_PRIMARY_ENTRY_TYPES = frozenset({'Domain', 'Repeat'})
_SITE_ENTRY_TYPES = frozenset({'Active_site', 'Binding_site', 'Conserved_site', 'PTM'})
# A lower-tier entry is treated as redundant (and dropped) only when MORE THAN
# this fraction of its residues is already covered by higher-tier entries. The
# "majority" midpoint keeps a member-DB/family entry that annotates a region no
# domain covers (e.g. KLF1's cd21581 N-terminal domain, which merely abuts the
# zinc fingers), while still dropping a co-hit that sits on top of a domain
# (e.g. TRIB2's two G3DSA kinase lobes over the Pkinase family).
_REDUNDANT_COVER = 0.5


def _aa_overlap(s1, e1, s2, e2):
    """True if [s1,e1] and [s2,e2] overlap by at least one residue."""
    return s1 <= e2 and s2 <= e1


def _covered_fraction(s, e, intervals):
    """Fraction of residues in [s,e] covered by the union of `intervals`."""
    length = e - s + 1
    if length <= 0:
        return 0.0
    segs = sorted((max(s, a), min(e, b)) for a, b in intervals if a <= e and b >= s)
    covered = 0
    cur = s - 1
    for a, b in segs:
        a = max(a, cur + 1)
        if b >= a:
            covered += b - a + 1
            cur = b
    return covered / length


def filter_representative_domains(df_domains):
    """
    Reduce a single transcript's representative domains to a clean domain set
    using the curated InterPro entry `type`, replacing the geometric
    collapse_contained_domains() heuristic (the 2 AA / 85% overlap rule).

    Priority ladder - per region, keep the best available annotation, and drop a
    lower tier only where a higher tier already covers the MAJORITY of it
    (>_REDUNDANT_COVER). This is "demote, don't delete": a member-DB or family
    entry that annotates a region no domain covers is kept, so we never lose the
    only annotation for a region (e.g. KLF1's cd21581, which carries the change).

      Tier 1 PRIMARY : InterPro Domain / Repeat
      Tier 2 PARENT  : InterPro Family / Homologous_superfamily (+ IPR of unknown
                       type) - kept unless the majority is covered by PRIMARY
      Tier 3 MEMBER  : non-InterPro member-DB hits (G3DSA/PTHR/SSF/cd/PF/...) -
                       kept unless the majority is covered by kept PRIMARY+PARENT
      dropped        : InterPro site/PTM types (residue features, not domains)

    Then collapse genuine duplicates: two kept rows with the SAME domain_id that
    overlap -> keep the longer (same accession at DISJOINT positions, e.g. two
    RRM instances, is kept as two domains). Return sorted by AA_start.

    Deliberately NOT handled here (left to the HMM/Pfam enrichment layer): cross-
    transcript identity when canonical and compared annotate one physical domain
    under different accessions (Pfam family from hmmscan is the stable key), and
    surfacing an event region with no InterPro Domain entry (hmmscan recovers it).

    Requires 'domain_id' and 'type' columns. If either is absent, or `type` is
    entirely NULL (DomainEvent/DomainType fallback rows, or a DB built before the
    type column existed), the frame is returned unchanged.
    """
    if len(df_domains) <= 1:
        return df_domains
    if 'domain_id' not in df_domains.columns or 'type' not in df_domains.columns:
        return df_domains
    if df_domains['type'].isna().all():
        return df_domains

    df = df_domains
    dom_id = df['domain_id'].astype(str)
    is_ipr = dom_id.str.startswith('IPR')
    etype = df['type']
    starts = df['AA_start']
    ends = df['AA_end']

    is_primary = is_ipr & etype.isin(_PRIMARY_ENTRY_TYPES)
    is_site = is_ipr & etype.isin(_SITE_ENTRY_TYPES)
    primary_idx = df.index[is_primary].tolist()
    parent_idx = df.index[is_ipr & ~is_primary & ~is_site].tolist()  # Family/HSF/unknown-type IPR
    member_idx = df.index[~is_ipr].tolist()                          # non-IPR member DBs

    keep = list(primary_idx)
    prim_iv = [(starts[i], ends[i]) for i in primary_idx]
    for pi in parent_idx:                                            # tier 2
        if _covered_fraction(starts[pi], ends[pi], prim_iv) <= _REDUNDANT_COVER:
            keep.append(pi)
    higher_iv = [(starts[i], ends[i]) for i in keep]                 # PRIMARY + kept PARENT
    for mi in member_idx:                                            # tier 3
        if _covered_fraction(starts[mi], ends[mi], higher_iv) <= _REDUNDANT_COVER:
            keep.append(mi)

    # collapse genuine duplicates (same accession, overlapping) -> keep the longer
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
            if _aa_overlap(starts[ia], ends[ia], starts[ib], ends[ib]):
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
    Determine the genomic window around the differing junctions (refined over
    three rounds, per the domain-boundary algorithm) and return the domains of
    the canonical and compared transcript that fall within it.

    Returns (t_domains_in_region, c_domains_in_region), each with a 'length'
    column added (AA_end - AA_start + 1).
    """
    junction_idxs = canonical_junctions | transcript_junctions
    min_bp = min(junctions[idx][0] for idx in junction_idxs)
    max_bp = max(junctions[idx][1] for idx in junction_idxs)

    t_exons = transcript_exons[transcript_id]
    c_exons = transcript_exons[canonical_transcript_id]

    # Round 1: window spanning the boundary exons of the differing junctions.
    t_first_exon, t_last_exon = find_boundary_exons(t_exons, min_bp, max_bp)
    c_first_exon, c_last_exon = find_boundary_exons(c_exons, min_bp, max_bp)

    t_min_aa, t_max_aa = get_aa_range(t_first_exon, t_last_exon)
    c_min_aa, c_max_aa = get_aa_range(c_first_exon, c_last_exon)

    # collapse_contained_domains() (the geometric 2 AA / 85% overlap heuristic)
    # is retained above but no longer called; the domain set is now reduced by
    # curated InterPro entry type via filter_representative_domains().
    df_t_domains = filter_representative_domains(domain_lookup(transcript_id))
    df_c_domains = filter_representative_domains(domain_lookup(canonical_transcript_id))

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

    # _domains_in_aa_range() no longer defensively .copy()s its result
    # (removed as a redundant allocation on a function called up to 4x per
    # compared transcript - boolean-mask indexing already returns a new,
    # independent DataFrame). Its return value still carries pandas'
    # internal "copy of a slice" provenance marker though, which trips
    # SettingWithCopyWarning on the ['length'] assignment below regardless
    # of .loc[] usage (the warning is provenance-based, not a real aliasing
    # check). Since only the round2 frames actually get mutated - not every
    # _domains_in_aa_range() call - .copy() them here, once, right before
    # the mutation, instead of inside the 4x-per-pair helper.
    t_domains_round2 = t_domains_round2.copy()
    c_domains_round2 = c_domains_round2.copy()
    t_domains_round2['length'] = t_domains_round2['AA_end'] - t_domains_round2['AA_start'] + 1
    c_domains_round2['length'] = c_domains_round2['AA_end'] - c_domains_round2['AA_start'] + 1
    return t_domains_round2, c_domains_round2


# ---------------------------------------------------------------------------
# Phase 3: domain comparison and classification
# ---------------------------------------------------------------------------

def domain_name_set(row, name_cols=DOMAIN_NAME_COLUMNS):
    """Return the set of non-empty/non-null identifier names for a domain row.

    Cell values may be a single identifier or a ';'-joined list of identifiers
    merged by collapse_contained_domains()/_merge_domain_names() - split them
    back out here so a merged name (e.g. "G3DSA:2.80.10.50; IPR002209") still
    matches a plain "IPR002209" cell elsewhere, as collapse_contained_domains'
    docstring intends.
    """
    names = set()
    for col in name_cols:
        val = row[col]
        if val is None or pd.isna(val):
            continue
        for name in str(val).split(';'):
            name = name.strip()
            if name and name not in ('None', 'nan'):
                names.add(name)
    return names


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
        return 'same'
    return 'longer' if t_length > c_length else 'shorter'


def classify_domain_change(c_count, t_count, c_length, t_length):
    """
    Classify a group of matched domains (c_count in canonical, t_count in the
    compared transcript) into one of the Phase 3 result categories.
    """
    if c_count == 0:
        return 'added_domain'
    if t_count == 0:
        return 'dropped domain'
    if c_count == 1 and t_count == 1:
        return classify_length_pair(t_length, c_length)
    if c_count == 1:
        return 'split domain'
    if t_count == 1:
        return 'merged domain'
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
    - C=0, T>0          -> 'new domain'
    - C>0, T=0          -> 'dropped domain'
    - C=1, T=1          -> 'same' / 'longer' / 'shorter' (by length)
    - C=1, T>1          -> 'split domain'
    - C>1, T=1          -> 'merged domain'
    - C>1, T>1, C==T    -> 'same_domains' / 'longer_domains' / 'shorter_domains' (by total length)
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
            'transcript_id': transcript_id,
            'domain_name': choose_domain_display_name(names),
            'canonical_domain_length': c_length,
            'transcript_domain_length': t_length,
            'canonical_domains_number': c_count,
            'transcript_domains_number': t_count,
        }


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
        self.events = []

    def add_event(self, event, transcript_id=None, domain_name=None, canonical_domain_length=None, transcript_domain_length=None,
                  canonical_domains_number=None, transcript_domains_number=None, is_longest_cds=None, is_most_like_canonical=None):
        self.events.append((event, transcript_id, domain_name, canonical_domain_length,
                            transcript_domain_length, canonical_domains_number, transcript_domains_number,
                            is_longest_cds, is_most_like_canonical))

    def analyze_junction(self, df_gene_transcripts, canonical_transcript_ids, exon_lookup, domain_lookup):
        """
        Run the DOMAS algorithm for this cluster:

        Phase 1 - find the canonical transcript and, for every other transcript
        of the gene, the junctions (if any) that match its exon structure and
        the subset of those that are unique compared to the canonical transcript.

        Phase 1.5 - tag each comparable transcript with whether it's the one the
        longest-CDS rule and/or the most-like-canonical rule (see
        select_most_like_canonical()) would pick - none are skipped. Exactly one
        transcript is tagged is_longest_cds; is_most_like_canonical is left unset
        across the whole cluster when no transcript qualifies.

        Phase 2/3 - for each transcript with a unique junction, determine the
        relevant genomic window and compare its domains against the canonical
        transcript's, recording one event per domain group.
        """
        # Use an order-preserving dedup (not `set`) so the order in which
        # transcripts are processed - and therefore the order of the output
        # rows - doesn't depend on Python's per-process string hash seed.
        # Invalid placeholder ids (e.g. NaN for transcripts with neither an
        # ensembl nor a refseq id) are dropped so they can't spuriously match
        # another gene's similarly-invalid "canonical" id.
        invalid_ids = {'', 'nan', 'None'}
        gene_transcript_ids = [
            tid for tid in dict.fromkeys(
                df_gene_transcripts.transcript_ensembl_id.fillna(df_gene_transcripts.transcript_refseq_id)
            )
            if tid is not None and not pd.isna(tid) and tid not in invalid_ids
        ]

        # CDS length per transcript id, for the longest-CDS tag below. Missing
        # cds_start/cds_end yields NaN here, normalized to -1 so such transcripts
        # are deprioritized rather than breaking the max() comparison.
        combined_ids = df_gene_transcripts.transcript_ensembl_id.fillna(df_gene_transcripts.transcript_refseq_id)
        cds_length_by_transcript = {
            tid: (length if pd.notna(length) else -1)
            for tid, length in zip(combined_ids, df_gene_transcripts.cds_end - df_gene_transcripts.cds_start)
        }

        # Check if gene exists in the database at all
        if df_gene_transcripts.empty:
            self.add_event('gene_not_in_db')
            logger.debug(f"Gene {self.gene_ensembl_id} ({self.gene_symbol}) not found in database for cluster {self.cluster_name}, specie {self.specie}. Skipping analysis.")
            return

        gene_canonical_ids = canonical_transcript_ids.intersection(gene_transcript_ids)
        if not gene_canonical_ids:
            self.add_event('no_canonical_transcript')
            logger.debug(f"No canonical transcript found for cluster {self.cluster_name}, specie {self.specie}. Skipping analysis.")
            return
        if len(gene_transcript_ids) == 1:
            self.add_event('only_one_transcript')
            logger.debug(f"Only one transcript found for cluster {self.cluster_name}, specie {self.specie}.")
            return
        if len(gene_canonical_ids) > 1:
            logger.warning(f"Multiple canonical transcripts found for cluster {self.cluster_name}, specie {self.specie}: "
                f"{sorted(gene_canonical_ids)}. Using the first one (sorted)."
            )
        self.canonical_transcript_id = sorted(gene_canonical_ids)[0]

        transcript_exons = {
            transcript_id: exon_lookup(transcript_id)
            for transcript_id in gene_transcript_ids
        }

        transcript_junctions = {
            transcript_id: find_matching_junction_indices(exons, self.junctions, strand=self.strand or '+')
            for transcript_id, exons in transcript_exons.items()
        }

        for idx, junction in enumerate(self.junctions):
            if not any(idx in junction_idxs for junction_idxs in transcript_junctions.values()):
                logger.debug(f"Junction {junction} in cluster {self.cluster_name} does not map to any transcript. ")
                self.add_event('junction_not_mapped', None)

        canonical_junctions = transcript_junctions.get(self.canonical_transcript_id, set())
        if not canonical_junctions:
            self.add_event('no_canonical_junctions')
            logger.debug(f"No canonical junctions found for cluster {self.cluster_name}, specie {self.specie}. Skipping analysis.")
            return

        comparable_transcript_ids = []
        for transcript_id, junction_idxs in transcript_junctions.items():
            if transcript_id == self.canonical_transcript_id:
                continue
            if not junction_idxs:
                logger.debug(f"Transcript {transcript_id} in cluster {self.cluster_name}, specie {self.specie} does not have any junctions. ")
                self.add_event('transcript_doesnt_have_junctions', transcript_id=transcript_id)
                continue

            unique_junctions = junction_idxs - canonical_junctions
            if not unique_junctions:
                logger.debug(
                    f"Transcript {transcript_id} in cluster {self.cluster_name}, specie {self.specie} does not have any unique junctions "
                    "compared to the canonical transcript. Skipping this transcript for comparison."
                )
                self.add_event('no_unique_junctions', transcript_id=transcript_id)
                continue

            comparable_transcript_ids.append(transcript_id)

        # For tagging only - every comparable transcript is compared below regardless.
        longest_cds_transcript_id = None
        most_like_canonical_transcript_id = None
        if comparable_transcript_ids:
            longest_cds_transcript_id = max(
                comparable_transcript_ids,
                key=lambda tid: (cds_length_by_transcript.get(tid, -1), tid),
            )
            most_like_canonical_transcript_id = select_most_like_canonical(
                comparable_transcript_ids, self.canonical_transcript_id, transcript_exons,
                self.junctions, cds_length_by_transcript,
            )

        for transcript_id in comparable_transcript_ids:
            junction_idxs = transcript_junctions[transcript_id]
            is_longest_cds = transcript_id == longest_cds_transcript_id
            is_most_like_canonical = transcript_id == most_like_canonical_transcript_id
            events = list(compare_domains(
                domain_lookup, transcript_exons, self.canonical_transcript_id, transcript_id,
                canonical_junctions, junction_idxs, self.junctions,
            ))
            if events:
                for event in events:
                    self.add_event(**event, is_longest_cds=is_longest_cds, is_most_like_canonical=is_most_like_canonical)
            else:
                self.add_event('no_domains_in_region', transcript_id=transcript_id,
                                is_longest_cds=is_longest_cds, is_most_like_canonical=is_most_like_canonical)

    def print_results(self, file_name='analysis_results.txt'):
        with open(file_name, 'a') as f:
            f.write(f"Cluster: {self.cluster_name}, Gene: {self.gene_symbol} ({self.gene_ensembl_id}), Chromosome: {self.chromosome}, Specie: {self.specie}\n")
            f.write(f"\tCanonical Transcript ID: {self.canonical_transcript_id}\n")
            f.write(f"\tJunctions: {self.junctions}\n")
            f.write(f"\tEvents:\n")
            for event, transcript_id, domain_name, canonical_domain_length, transcript_domain_length, \
                canonical_domains_number, transcript_domains_number in self.events:
                msg = f"\t\t{event}: Transcript ID={transcript_id}, Domain Name={domain_name}, "
                msg += f"Canonical Length={canonical_domain_length}, Transcript Length={transcript_domain_length}, "
                msg += f"Canonical Domains Number={canonical_domains_number}, Transcript Domains Number={transcript_domains_number}\n"
                f.write(msg)

    def get_results_df(self):
        return pd.DataFrame(
            self.events,
            columns=['event', 'transcript_id', 'domain_name', 'c_domain_length', 't_domain_length',
                        'c_domains_number', 't_domains_number', 'is_longest_cds', 'is_most_like_canonical']
        )
        

def _analyze_single_cluster(cluster_tuple, exon_lookup=None, domain_lookup=None, canonical_transcript_ids=None,
                           gene_strand=None, transcripts_by_gene=None, empty_transcripts=None):
    """Analyze a single cluster."""
    _, cluster_df = cluster_tuple

    gene_ensembl_id = cluster_df.gene_ensembl_id.iat[0]
    gene_symbol = cluster_df.gene_symbol.iat[0]
    event_type = cluster_df.event_type.iat[0] if 'event_type' in cluster_df.columns else None
    specie = cluster_df.specie.iat[0] if 'specie' in cluster_df.columns else None
    strand = gene_strand.get(gene_ensembl_id)

    cluster_result = ClusterAnalysisResult(
        cluster_df.cluster_name.iat[0], gene_ensembl_id, gene_symbol,
        as_event_type=event_type, specie=specie, strand=strand
    )
    cluster_result.junctions = list(zip(cluster_df['start_position'], cluster_df['end_position']))

    df_gene_transcripts = transcripts_by_gene.get(gene_ensembl_id, empty_transcripts)
    cluster_result.analyze_junction(df_gene_transcripts, canonical_transcript_ids, exon_lookup, domain_lookup)

    return cluster_result


# Events recorded for a transcript/cluster that was NOT a real comparison to the
# canonical transcript: cluster-level non-comparisons (gene/canonical/junction
# problems) plus per-transcript skips (no junctions / no unique junction). Rows
# carrying these events are the "non-comparable" / "not chosen" transcripts that
# analyze_junctions(filter_non_comparable=True) drops from the output CSV.
NON_COMPARISON_EVENTS = frozenset({
    'gene_not_in_db', 'no_canonical_transcript', 'only_one_transcript',
    'no_canonical_junctions', 'junction_not_mapped',
    'transcript_doesnt_have_junctions', 'no_unique_junctions',
})


def _csv_writer_worker(result_queue, output_path, df_results_columns, logger_instance=None,
                       filter_non_comparable=False):
    """
    Dedicated writer thread that processes results from a queue and writes to CSV.

    Runs continuously until it receives a None sentinel value.
    Writes results incrementally as they arrive from compute workers.

    filter_non_comparable: if True, rows whose event is in NON_COMPARISON_EVENTS
    (transcripts that were not actually compared to canonical) are dropped.
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

                    if filter_non_comparable:
                        df_cluster_results = df_cluster_results[
                            ~df_cluster_results['event'].isin(NON_COMPARISON_EVENTS)]
                        if df_cluster_results.empty:
                            continue

                    df_transformed = (
                        df_cluster_results
                        .rename(columns={'event': 'event_type'})
                        .assign(
                            cluster=cluster_result.cluster_name,
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
        # Cleanup temporary chunks
        shutil.rmtree(output_dir, ignore_errors=True)


# Per-worker-process state, populated once by _init_worker() rather than
# per chunk - df_exons/df_domains can be large, so re-pickling them for every
# chunk via submit() args would be wasted work.
_worker_state = {}


def _init_worker(df_exons, df_domains, canonical_transcript_ids, gene_strand, transcripts_by_gene, empty_transcripts):
    """ProcessPoolExecutor initializer - runs once when each worker process starts."""
    _worker_state['exon_lookup'] = build_exon_lookup(df_exons)
    _worker_state['domain_lookup'] = build_domain_lookup(df_domains)
    _worker_state['canonical_transcript_ids'] = canonical_transcript_ids
    _worker_state['gene_strand'] = gene_strand
    _worker_state['transcripts_by_gene'] = transcripts_by_gene
    _worker_state['empty_transcripts'] = empty_transcripts


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
    gene_strand = _worker_state['gene_strand']
    transcripts_by_gene = _worker_state['transcripts_by_gene']
    empty_transcripts = _worker_state['empty_transcripts']

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
                empty_transcripts=empty_transcripts,
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

    def _load_junctions_data(self, df_junctions):
        """Validate the junctions DataFrame. Reading junctions from a file (plain CSV,
        hadas-format Excel, IOE, ...) is alternative_splicing.py's responsibility -
        callers pass an already-loaded DataFrame here."""
        df_junctions = df_junctions.copy()

        if 'cluster_name' not in df_junctions.columns and 'cluster' in df_junctions.columns:
            df_junctions = df_junctions.rename(columns={'cluster': 'cluster_name'})

        required_columns = ['gene_ensembl_id', 'start_position', 'end_position', 'cluster_name']
        for col in required_columns:
            if col not in df_junctions.columns:
                raise ValueError(f"Column '{col}' is required in df_junctions but not found.")

        return df_junctions

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

    def _load_database_data(self, gene_ids, use_ensembl_only=False, use_representative_domains=False):
        """Load genes, transcripts, domains, and exons from database."""
        placeholders = ','.join('?' * len(gene_ids))
        df_genes = pd.read_sql_query(
            f"SELECT gene_ensembl_id, strand FROM Genes WHERE gene_ensembl_id IN ({placeholders})",
            self.con, params=gene_ids
        )
        gene_strand = dict(zip(df_genes['gene_ensembl_id'], df_genes['strand']))

        from alternative_splicing import get_genes_df_transcripts, get_domains_db, get_exons_for_transcripts
        df_transcripts = get_genes_df_transcripts(self.con, gene_ids)

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

        df_domains = get_domains_db(
            self.con, transcript_ids, use_representative_domains=use_representative_domains
        )

        df_exons = get_exons_for_transcripts(self.con, transcript_ids)

        return df_genes, df_transcripts, df_domains, df_exons, gene_strand

    def _prepare_lookup_structures(self, df_transcripts, df_exons, df_domains):
        """Build lookup structures and transcript groupings."""
        invalid_ids = {'', 'nan', 'None'}
        all_transcript_ids = df_transcripts.transcript_ensembl_id.fillna(df_transcripts.transcript_refseq_id)
        valid_mask = all_transcript_ids.notna() & ~all_transcript_ids.isin(invalid_ids)
        canonical_transcript_ids = set(all_transcript_ids[valid_mask & (df_transcripts.canonical != 0)])

        transcripts_by_gene = {gid: g for gid, g in df_transcripts.groupby('gene_ensembl_id')}
        empty_transcripts = df_transcripts.iloc[0:0]

        return canonical_transcript_ids, transcripts_by_gene, empty_transcripts

    def _prepare_cluster_groups(self, df_junctions):
        """Group junctions into clusters."""
        group_columns = ['specie', 'cluster_name'] if 'specie' in df_junctions.columns else ['cluster_name']
        cluster_groups = list(df_junctions.groupby(group_columns))
        return cluster_groups

    def _run_parallel_analysis(self, cluster_groups, df_exons, df_domains, canonical_transcript_ids,
                               gene_strand, transcripts_by_gene, empty_transcripts, num_workers, output_path,
                               filter_non_comparable=False):
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

        # as_completed() returns chunks in whatever order they finish, unrelated to
        # cluster_groups' original order. PDF file names embed a sequential count (see
        # _generate_pdfs), so results are restored to that original order below.
        has_specie_column = 'specie' in cluster_groups[0][1].columns if cluster_groups else False

        def _group_identity(cluster_name, specie):
            return (specie, cluster_name) if has_specie_column else (cluster_name,)

        original_order = {
            _group_identity(group_df['cluster_name'].iat[0], group_df['specie'].iat[0] if has_specie_column else None): idx
            for idx, (_, group_df) in enumerate(cluster_groups)
        }

        self.logger.info(f"Starting analysis with {actual_workers} workers + 1 writer thread")
        self.logger.info(f"Total clusters: {total}")
        self.logger.info(f"Processing {total_chunks} chunks (~{chunk_size} clusters per chunk)")

        # Prepare column definitions for CSV
        df_results_columns = ['cluster', 'gene_symbol', 'specie', 'event_type', 'canonical_transcript_id', 'transcript_id', 'domain_name',
                              'c_domain_length', 't_domain_length', 'c_domains_number', 't_domains_number',
                              'is_longest_cds', 'is_most_like_canonical']

        # Create queue for results (backpressure if writer lags)
        result_queue = queue.Queue(maxsize=actual_workers * 2)

        # Start dedicated writer thread
        writer_thread = threading.Thread(
            target=_csv_writer_worker,
            args=(result_queue, output_path, df_results_columns, self.logger, filter_non_comparable),
            daemon=False
        )
        writer_thread.start()

        chunks_with_info = [
            (chunk_idx % actual_workers, chunk_idx, total_chunks, chunk)
            for chunk_idx, chunk in enumerate(chunks)
        ]

        chunk_worker = _process_cluster_chunk

        all_results = []  # Collect results for PDF generation
        processed_count = 0
        last_time = time.perf_counter()

        with ProcessPoolExecutor(
            max_workers=actual_workers,
            initializer=_init_worker,
            initargs=(df_exons, df_domains, canonical_transcript_ids, gene_strand, transcripts_by_gene, empty_transcripts),
        ) as executor:
            # Submit all tasks - df_exons/df_domains are sent once per worker via
            # initargs above, not re-pickled per chunk here.
            futures = [executor.submit(chunk_worker, chunk_info) for chunk_info in chunks_with_info]

            # Process results as they complete: send to writer AND collect for PDFs
            for future in as_completed(futures):
                chunk_results = future.result()
                result_queue.put(chunk_results)  # Send to writer thread for CSV
                all_results.extend(chunk_results)  # Collect for PDF generation
                processed_count += len(chunk_results)

                if processed_count % 10000 == 0:
                    cur_time = time.perf_counter()
                    self.logger.info(
                        f"[Main] Analyzed {processed_count}/{total} clusters. "
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

        self.logger.info(f"[Main] Analysis complete: processed {processed_count}/{total} clusters")
        return all_results

    # Events recorded for a transcript that was NOT actually compared to the
    # canonical transcript (it was skipped for lacking junctions or lacking a
    # unique junction).
    _SKIPPED_TRANSCRIPT_EVENTS = {
        'transcript_doesnt_have_junctions', 'no_unique_junctions',
    }

    def _comparable_transcript_ids(self, cluster_result):
        """Canonical transcript id plus every transcript id that was actually compared to it."""
        comparable_ids = set()
        if cluster_result.canonical_transcript_id:
            comparable_ids.add(cluster_result.canonical_transcript_id)
        df_cluster_results = cluster_result.get_results_df()
        if len(df_cluster_results) > 0:
            compared_mask = ~df_cluster_results['event'].isin(self._SKIPPED_TRANSCRIPT_EVENTS)
            comparable_ids.update(df_cluster_results.loc[compared_mask, 'transcript_id'].dropna())
        return comparable_ids

    def _generate_pdfs(self, results, print_genes, restrict_to_comparable=False, use_representative_domains=False):
        """Generate PDF visualizations for clusters.

        restrict_to_comparable: if True, each PDF only draws the canonical
        transcript and the transcripts that were actually compared to it -
        transcripts skipped during analysis (no junctions, no unique junction,
        or lost the longest-CDS tie-break) are omitted entirely.
        use_representative_domains: if True, domains drawn in the PDFs come from
        the RepresentativeDomains table where available (matching the domain
        source used for analysis when the same flag is passed to
        analyze_junctions()), falling back to DomainEvent/DomainType per protein.
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
            use_representative_domains=use_representative_domains,
        )

        gene_visualizations = {}
        for count, cluster_result in enumerate(results, start=1):
            try:
                df_cluster_junctions = pd.DataFrame(cluster_result.junctions, columns=['start', 'end'])
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
                        use_representative_domains=use_representative_domains,
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
                    df_results=cluster_result.get_results_df(),
                    transcript_ids=transcript_ids,
                    no_comparison_note=no_comparison_note,
                )
            except ValueError as e:
                self.logger.warning(f"Warning: Skipping PDF generation for {cluster_result.gene_symbol}, specie {cluster_result.specie}: {e}")

    def analyze_junctions(self, df_junctions, output_path='as_events_junctions_analysis.csv',
                          filter_transcript_count=0, create_pdf=True, print_genes=None,
                          num_workers=4, use_ensembl_only=False, restrict_pdf_to_comparable=False,
                          use_representative_domains=False, filter_non_comparable=False):
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
            use_representative_domains: If True, domains come from the
                RepresentativeDomains table instead of DomainEvent/DomainType.
                A protein with no RepresentativeDomains entry falls back to its
                DomainEvent/DomainType domains. If False (default), the
                algorithm is unchanged - domains come from DomainEvent/DomainType
                only, exactly as before.
            filter_non_comparable: If True, the output CSV contains only rows for
                transcripts that were actually compared to canonical - rows whose
                event is a non-comparison / skip event (see NON_COMPARISON_EVENTS)
                are dropped. The returned ClusterAnalysisResult objects and any
                PDFs are unaffected; only the written CSV is filtered.

        Returns:
            List of ClusterAnalysisResult objects
        """
        # Validate input
        df_junctions = self._load_junctions_data(df_junctions)
        df_junctions = self._filter_junctions_by_transcript_count(df_junctions, filter_transcript_count)

        gene_ids = df_junctions.gene_ensembl_id.unique().tolist()
        self.logger.info(f"Analyzing {len(gene_ids)} genes")

        # Load data from database
        df_genes, df_transcripts, df_domains, df_exons, gene_strand = self._load_database_data(
            gene_ids, use_ensembl_only=use_ensembl_only, use_representative_domains=use_representative_domains
        )

        # Prepare lookup structures
        canonical_transcript_ids, transcripts_by_gene, empty_transcripts = \
            self._prepare_lookup_structures(df_transcripts, df_exons, df_domains)

        # Group junctions into clusters
        cluster_groups = self._prepare_cluster_groups(df_junctions)

        # Run parallel analysis (with dedicated writer thread for CSV output)
        results = self._run_parallel_analysis(
            cluster_groups, df_exons, df_domains, canonical_transcript_ids,
            gene_strand, transcripts_by_gene, empty_transcripts, num_workers,
            output_path, filter_non_comparable=filter_non_comparable,
        )

        # Generate PDFs if requested
        if create_pdf:
            self._generate_pdfs(
                results, print_genes, restrict_to_comparable=restrict_pdf_to_comparable,
                use_representative_domains=use_representative_domains,
            )

        return results
