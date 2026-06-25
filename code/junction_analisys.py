import logging
import os
import queue
import shutil
import tempfile
import threading
import time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

from generate_gene_pdf import GeneVisualization, prepare_gene_data_bulk
import domas


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
# Phase 2: domain boundary determination via coordinate matching
# ---------------------------------------------------------------------------

def find_boundary_exons(df_exons, min_bp, max_bp):
    """
    Return (first_exon, last_exon): the exon whose genomic end is closest to
    (but not before) min_bp, and the exon whose genomic start is closest to
    (but not after) max_bp, allowing a 1bp tolerance on both ends.
    """
    first_candidates = df_exons[df_exons['genomic_end_tx'] >= min_bp - 1]
    last_candidates = df_exons[df_exons['genomic_start_tx'] <= max_bp + 1]
    first_exon = first_candidates.loc[first_candidates['genomic_end_tx'].idxmin()]
    last_exon = last_candidates.loc[last_candidates['genomic_start_tx'].idxmax()]
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
    """
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
    Genomic (start, end) bp positions corresponding to the start of the
    nearest-starting domain and the end of the furthest-reaching domain in
    `domains_in_region`, or (None, None) if there are no domains with a
    defined AA range.
    """
    domains = domains_in_region[domains_in_region['AA_start'] != 0]
    if domains.empty:
        return None, None

    min_domain_bp = domains['AA_start'].min() * 3
    max_domain_bp = domains['AA_end'].max() * 3

    first_exon = df_exons[(df_exons['abs_start_CDS'] <= min_domain_bp) & (df_exons['abs_end_CDS'] >= min_domain_bp)]
    last_exon = df_exons[(df_exons['abs_start_CDS'] <= max_domain_bp) & (df_exons['abs_end_CDS'] >= max_domain_bp)]
    if first_exon.empty or last_exon.empty:
        return None, None

    min_bp = _cds_to_bp(first_exon.iloc[0], min_domain_bp)
    max_bp = _cds_to_bp(last_exon.iloc[0], max_domain_bp)
    return min_bp, max_bp


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

    def lookup(transcript_id):
        a = by_ensembl.get(transcript_id)
        b = by_refseq.get(transcript_id)
        if a is not None and b is not None:
            return pd.concat([a, b]).drop_duplicates()
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
    return df_domains[(df_domains['AA_end'] >= min_aa) & (df_domains['AA_start'] <= max_aa)].copy()


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


def collapse_contained_domains(df_domains, tolerance=2):
    """
    Collapse domains (of a single transcript) that contain each other within
    `tolerance` AA on either side into a single row - the longer domain -
    merging the redundant domains' identifier names into the kept row so
    compare_domains() can match on any of them.
    """
    if len(df_domains) <= 1:
        return df_domains

    df_domains = df_domains.copy()
    starts = df_domains['AA_start'].to_numpy()
    ends = df_domains['AA_end'].to_numpy()
    index = df_domains.index.to_numpy()
    order_labels = (df_domains['AA_end'] - df_domains['AA_start']).sort_values(ascending=False).index.tolist()
    pos_of_label = {label: pos for pos, label in enumerate(df_domains.index)}
    order = [pos_of_label[label] for label in order_labels]

    dropped = set()
    merges = {}  # position -> list of positions merged into it
    for oi in order:
        if oi in dropped:
            continue
        for oj in order:
            if oj == oi or oj in dropped:
                continue
            if starts[oi] - tolerance <= starts[oj] and ends[oj] <= ends[oi] + tolerance:
                merges.setdefault(oi, []).append(oj)
                dropped.add(oj)

    if merges:
        name_block = df_domains[DOMAIN_NAME_COLUMNS].to_numpy(dtype=object)
        for oi, others in merges.items():
            for ci in range(len(DOMAIN_NAME_COLUMNS)):
                values = [name_block[oi, ci]] + [name_block[oj, ci] for oj in others]
                name_block[oi, ci] = _merge_domain_names(*values)
        for ci, col in enumerate(DOMAIN_NAME_COLUMNS):
            df_domains[col] = name_block[:, ci]

    return df_domains.drop(index=index[list(dropped)] if dropped else [])


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

    df_t_domains = collapse_contained_domains(domain_lookup(transcript_id))
    df_c_domains = collapse_contained_domains(domain_lookup(canonical_transcript_id))

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

    t_domains_round2['length'] = t_domains_round2['AA_end'] - t_domains_round2['AA_start'] + 1
    c_domains_round2['length'] = c_domains_round2['AA_end'] - c_domains_round2['AA_start'] + 1
    return t_domains_round2, c_domains_round2


# ---------------------------------------------------------------------------
# Phase 3: domain comparison and classification
# ---------------------------------------------------------------------------

def domain_name_set(row, name_cols=DOMAIN_NAME_COLUMNS):
    """Return the set of non-empty/non-null identifier names for a domain row."""
    names = set()
    for col in name_cols:
        val = row[col]
        if val is None or pd.isna(val):
            continue
        text = str(val).strip()
        if text and text not in ('None', 'nan'):
            names.add(val)
    return names


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
    intervals = sorted((df.loc[i, start_col], df.loc[i, end_col]) for i in idxs)
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

    canonical_names = {idx: domain_name_set(row) for idx, row in c_domains.iterrows()}
    transcript_names = {idx: domain_name_set(row) for idx, row in t_domains.iterrows()}

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
                  canonical_domains_number=None, transcript_domains_number=None):
        self.events.append((event, transcript_id, domain_name, canonical_domain_length,
                            transcript_domain_length, canonical_domains_number, transcript_domains_number))

    def analyze_junction(self, df_gene_transcripts, canonical_transcript_ids, exon_lookup, domain_lookup):
        """
        Run the DOMAS algorithm for this cluster:

        Phase 1 - find the canonical transcript and, for every other transcript
        of the gene, the junctions (if any) that match its exon structure and
        the subset of those that are unique compared to the canonical transcript.

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
                self.add_event('junction_not_mapped', idx)

        canonical_junctions = transcript_junctions.get(self.canonical_transcript_id, set())
        if not canonical_junctions:
            self.add_event('no_canonical_junctions')
            logger.debug(f"No canonical junctions found for cluster {self.cluster_name}, specie {self.specie}. Skipping analysis.")
            return

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

            events = list(compare_domains(
                domain_lookup, transcript_exons, self.canonical_transcript_id, transcript_id,
                canonical_junctions, junction_idxs, self.junctions,
            ))
            if events:
                for event in events:
                    self.add_event(**event)
            else:
                self.add_event('no_domains_in_region', transcript_id=transcript_id)

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
                        'c_domains_number', 't_domains_number']
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


def _csv_writer_worker(result_queue, output_path, df_results_columns, logger_instance=None):
    """
    Dedicated writer thread that processes results from a queue and writes to CSV.

    Runs continuously until it receives a None sentinel value.
    Writes results incrementally as they arrive from compute workers.
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
                    df_chunk = pd.concat(result_frames, ignore_index=True)

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

        # Combine all chunks into final CSV
        chunk_files = sorted([os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith('.csv')])

        if chunk_files:
            log.info(f"[Writer] Combining {len(chunk_files)} chunks...")
            df_all_results = pd.concat([pd.read_csv(f) for f in chunk_files], ignore_index=True)
            df_all_results.to_csv(output_path, index=False)
            log.info(f"[Writer] Final CSV written: {output_path} ({len(df_all_results)} rows)")
        else:
            df_all_results = pd.DataFrame(columns=df_results_columns)
            df_all_results.to_csv(output_path, index=False)
            log.info(f"[Writer] Empty results CSV written: {output_path}")

    finally:
        # Cleanup temporary chunks
        shutil.rmtree(output_dir, ignore_errors=True)


def _process_cluster_chunk(chunk_info, df_exons=None, df_domains=None, canonical_transcript_ids=None,
                          gene_strand=None, transcripts_by_gene=None, empty_transcripts=None):
    """Process a chunk of clusters sequentially (worker function for ProcessPoolExecutor).

    Builds lookups inside the worker to avoid pickling issues with closures.

    Args:
        chunk_info: tuple of (worker_id, chunk_index, total_chunks, chunk)
    """
    worker_id, chunk_index, total_chunks, chunk = chunk_info
    chunk_size = len(chunk)

    # Build lookups once per worker process
    exon_lookup = build_exon_lookup(df_exons)
    domain_lookup = build_domain_lookup(df_domains)

    logger.info(f"[Worker {worker_id}] Processing chunk {chunk_index + 1}/{total_chunks} ({chunk_size} clusters)")

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
                empty_transcripts=empty_transcripts
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

    logger.info(f"[Worker {worker_id}] Completed chunk {chunk_index + 1}/{total_chunks} ({chunk_size} clusters)")
    return chunk_results


class JunctionsAnalysis:
    def __init__(self, con, logger_instance=None, gene_visualization_cls=GeneVisualization):
        self.con = con
        self.logger = logger_instance or logger
        self.gene_visualization_cls = gene_visualization_cls

    def _load_junctions_data(self, junctions_csv, df_junctions, hadas_format):
        """Load and validate junctions from CSV or DataFrame."""
        if df_junctions is None and junctions_csv is None:
            raise ValueError("Either df_junctions or junctions_csv must be provided.")
        elif df_junctions is not None and junctions_csv is not None:
            raise ValueError("Only one of df_junctions or junctions_csv should be provided.")
        elif junctions_csv is not None:
            if hadas_format:
                df_junctions = domas.read_hadas_input_file(junctions_csv)
            else:
                df_junctions = pd.read_csv(junctions_csv)
        else:
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

        df_transcripts = pd.read_sql_query('select * from transcripts', self.con)
        gene_count = df_transcripts.value_counts('gene_GeneID_id')
        df_transcripts['gene_count'] = df_transcripts['gene_GeneID_id'].map(gene_count)
        genes_with_count = df_transcripts[df_transcripts['gene_count'] == filter_transcript_count].gene_GeneID_id.unique().tolist()
        return df_junctions[df_junctions['gene_ensembl_id'].isin(genes_with_count)]

    def _load_database_data(self, gene_ids):
        """Load genes, transcripts, domains, and exons from database."""
        placeholders = ','.join('?' * len(gene_ids))
        df_genes = pd.read_sql_query(
            f"SELECT gene_ensembl_id, strand FROM Genes WHERE gene_ensembl_id IN ({placeholders})",
            self.con, params=gene_ids
        )
        gene_strand = dict(zip(df_genes['gene_ensembl_id'], df_genes['strand']))

        df_transcripts = domas.get_genes_df_transcripts(self.con, gene_ids)

        # Combine ensembl and refseq IDs, filtering out invalid entries
        invalid_ids = {'', 'nan', 'None'}
        transcript_ids = set(
            df_transcripts.transcript_ensembl_id.fillna(df_transcripts.transcript_refseq_id)
        ) - {None} - invalid_ids

        df_domains = domas.get_transcript_domains_db(self.con, transcript_ids)

        from alternative_splicing import get_exons_for_transcripts
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
                               gene_strand, transcripts_by_gene, empty_transcripts, num_workers, output_path):
        """Execute cluster analysis in parallel with dedicated writer thread."""
        total = len(cluster_groups)

        # Divide clusters into chunks (one per worker, as many as needed)
        actual_workers = min(num_workers, total)
        chunk_size = max(1, (total + actual_workers - 1) // actual_workers)

        # Create chunks with only as many as needed (no empty chunks)
        chunks = []
        for i in range(0, total, chunk_size):
            chunks.append(cluster_groups[i:i+chunk_size])

        total_chunks = len(chunks)

        self.logger.info(f"Starting analysis with {actual_workers} workers + 1 writer thread")
        self.logger.info(f"Total clusters: {total}")
        self.logger.info(f"Processing {total_chunks} chunks (~{chunk_size} clusters per chunk)")

        # Prepare column definitions for CSV
        df_results_columns = ['cluster', 'gene_symbol', 'specie', 'event_type', 'canonical_transcript_id', 'transcript_id', 'domain_name',
                              'c_domain_length', 't_domain_length', 'c_domains_number', 't_domains_number']

        # Create queue for results (backpressure if writer lags)
        result_queue = queue.Queue(maxsize=actual_workers * 2)

        # Start dedicated writer thread
        writer_thread = threading.Thread(
            target=_csv_writer_worker,
            args=(result_queue, output_path, df_results_columns, self.logger),
            daemon=False
        )
        writer_thread.start()

        chunks_with_info = [
            (chunk_idx % total_chunks, chunk_idx, total_chunks, chunk)
            for chunk_idx, chunk in enumerate(chunks)
        ]

        from functools import partial
        chunk_worker = partial(
            _process_cluster_chunk,
            df_exons=df_exons,
            df_domains=df_domains,
            canonical_transcript_ids=canonical_transcript_ids,
            gene_strand=gene_strand,
            transcripts_by_gene=transcripts_by_gene,
            empty_transcripts=empty_transcripts
        )

        all_results = []  # Collect results for PDF generation
        processed_count = 0
        last_time = time.perf_counter()

        with ProcessPoolExecutor(max_workers=actual_workers) as executor:
            # Submit all tasks
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

        self.logger.info(f"[Main] Analysis complete: processed {processed_count}/{total} clusters")
        return all_results

def _generate_pdfs(self, results, print_genes):
        """Generate PDF visualizations for clusters."""
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
            self.con, [cluster_result.gene_ensembl_id for cluster_result in filtered_results if _gene_symbol_key(cluster_result.gene_symbol) is not None]
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
                    viz = self.gene_visualization_cls(self.con, gene_symbol, preloaded=preloaded)
                    gene_visualizations[gene_ensembl_id] = viz

                file_name = f'{cluster_result.gene_symbol}_{count}_junction_comparison.pdf'
                if cluster_result.as_event_type is not None:
                    file_name = f'{cluster_result.as_event_type}_{file_name}'
                if cluster_result.specie is not None:
                    file_name = f'{cluster_result.specie}_{file_name}'

                viz.create_pdf(
                    file_name,
                    protein_only=False,
                    domains_only=False,
                    df_junction=df_cluster_junctions,
                    df_results=cluster_result.get_results_df(),
                )
            except ValueError as e:
                self.logger.warning(f"Warning: Skipping PDF generation for {cluster_result.gene_symbol}, specie {cluster_result.specie}: {e}")

    def analyze_junctions(self, junctions_csv='as_events_junctions.csv', output_path='as_events_junctions_analysis.csv', df_junctions=None,
                          hadas_format=False, filter_transcript_count=0, create_pdf=True, print_genes=None, num_workers=4):
        """
        Analyze junctions and detect domain changes across alternative transcripts.

        Args:
            junctions_csv: Path to junction CSV file (if df_junctions not provided)
            output_path: Path for output CSV file
            df_junctions: DataFrame of junctions (if not reading from CSV)
            hadas_format: Whether to read CSV in HADAS format
            filter_transcript_count: If > 0, only analyze genes with exactly this many transcripts
            create_pdf: Whether to generate PDF visualizations
            print_genes: List of gene symbols to generate PDFs for (or all if None)
            num_workers: Number of parallel workers for analysis

        Returns:
            List of ClusterAnalysisResult objects
        """
        # Load and validate input
        df_junctions = self._load_junctions_data(junctions_csv, df_junctions, hadas_format)
        df_junctions = self._filter_junctions_by_transcript_count(df_junctions, filter_transcript_count)

        gene_ids = df_junctions.gene_ensembl_id.unique().tolist()
        self.logger.info(f"Analyzing {len(gene_ids)} genes")

        # Load data from database
        df_genes, df_transcripts, df_domains, df_exons, gene_strand = self._load_database_data(gene_ids)

        # Prepare lookup structures
        canonical_transcript_ids, transcripts_by_gene, empty_transcripts = \
            self._prepare_lookup_structures(df_transcripts, df_exons, df_domains)

        # Group junctions into clusters
        cluster_groups = self._prepare_cluster_groups(df_junctions)

        # Run parallel analysis (with dedicated writer thread for CSV output)
        results = self._run_parallel_analysis(
            cluster_groups, df_exons, df_domains, canonical_transcript_ids,
            gene_strand, transcripts_by_gene, empty_transcripts, num_workers,
            output_path
        )

        # Generate PDFs if requested
        if create_pdf:
            self._generate_pdfs(results, print_genes)

        return results
