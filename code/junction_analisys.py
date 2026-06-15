import logging
import time
import numpy as np
import pandas as pd

from generate_gene_pdf import GeneVisualization, prepare_gene_data_bulk
import domas


logger = logging.getLogger(__name__)

DOMAIN_NAME_COLUMNS = ['interpro', 'pfam', 'cdd', 'smart', 'tigr', 'CDD_id']
DOMAIN_NAME_PREFIX_PRIORITY = ['IPR', 'pfam', 'cd', 'smart', 'tigr', 'CDD']


# ---------------------------------------------------------------------------
# Phase 1: junction <-> exon matching
# ---------------------------------------------------------------------------

def find_matching_junction_indices(df_transcript_exons, junctions):
    """
    Return the set of indices (into `junctions`) of junctions that match this
    transcript's exon structure.

    A junction (start_position, end_position) matches the transcript if there
    are two exons, adjacent in transcript order, such that one exon's
    genomic_end_tx is within 1bp of the junction's start_position and the
    other exon's genomic_start_tx is within 1bp of the junction's end_position.
    """
    if df_transcript_exons.empty or not junctions:
        return set()

    exon_starts = df_transcript_exons['genomic_start_tx'].to_numpy()
    exon_ends = df_transcript_exons['genomic_end_tx'].to_numpy()
    exon_orders = df_transcript_exons['order_in_transcript'].to_numpy()

    matched = set()
    for idx, (start_position, end_position) in enumerate(junctions):
        upstream_orders = exon_orders[np.abs(exon_ends - start_position) <= 1]
        downstream_orders = exon_orders[np.abs(exon_starts - end_position) <= 1]
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
        return 'new domain'
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
    def __init__(self, cluster_name, gene_ensembl_id, gene_symbol, chromosome=None, as_event_type=None, specie=None):
        self.cluster_name = cluster_name
        self.gene_ensembl_id = gene_ensembl_id
        self.gene_symbol = gene_symbol
        self.chromosome = chromosome
        self.as_event_type = as_event_type
        self.specie = specie
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
            logger.warning(f"No canonical transcript found for cluster {self.cluster_name}, specie {self.specie}. Skipping analysis.")
            return
        if len(gene_transcript_ids) == 1:
            self.add_event('only_one_transcript')
            logger.info(f"Only one transcript found for cluster {self.cluster_name}, specie {self.specie}.")
            return
        if len(gene_canonical_ids) > 1:
            logger.warning(
                f"Multiple canonical transcripts found for cluster {self.cluster_name}, specie {self.specie}: "
                f"{sorted(gene_canonical_ids)}. Using the first one (sorted)."
            )
        self.canonical_transcript_id = sorted(gene_canonical_ids)[0]

        transcript_exons = {
            transcript_id: exon_lookup(transcript_id)
            for transcript_id in gene_transcript_ids
        }

        transcript_junctions = {
            transcript_id: find_matching_junction_indices(exons, self.junctions)
            for transcript_id, exons in transcript_exons.items()
        }

        for idx, junction in enumerate(self.junctions):
            if not any(idx in junction_idxs for junction_idxs in transcript_junctions.values()):
                logger.warning(
                    f"Junction {junction} in cluster {self.cluster_name} does not map to any transcript. "
                    "This may indicate an issue with the junction definition or exon mapping."
                )
                self.add_event('junction_not_mapped', idx)

        canonical_junctions = transcript_junctions.get(self.canonical_transcript_id, set())
        if not canonical_junctions:
            self.add_event('no_canonical_junctions')
            logger.warning(f"No canonical junctions found for cluster {self.cluster_name}, specie {self.specie}. Skipping analysis.")
            return

        for transcript_id, junction_idxs in transcript_junctions.items():
            if transcript_id == self.canonical_transcript_id:
                continue
            if not junction_idxs:
                logger.warning(
                    f"Transcript {transcript_id} in cluster {self.cluster_name}, specie {self.specie} does not have any junctions. "
                    "This may indicate an issue with the exon mapping for the canonical transcript."
                )
                self.add_event('transcript_doesnt_have_junctions', transcript_id=transcript_id)
                continue

            unique_junctions = junction_idxs - canonical_junctions
            if not unique_junctions:
                logger.info(
                    f"Transcript {transcript_id} in cluster {self.cluster_name}, specie {self.specie} does not have any unique junctions "
                    "compared to the canonical transcript. Skipping this transcript for comparison."
                )
                self.add_event('no_unique_junctions', transcript_id=transcript_id)
                continue

            for event in compare_domains(
                domain_lookup, transcript_exons, self.canonical_transcript_id, transcript_id,
                canonical_junctions, junction_idxs, self.junctions,
            ):
                self.add_event(**event)

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
        

class JunctionsAnalysis:
    def __init__(self, con, logger_instance=None, gene_visualization_cls=GeneVisualization):
        self.con = con
        self.logger = logger_instance or logger
        self.gene_visualization_cls = gene_visualization_cls

    def analyze_junctions(self, junctions_csv='as_events_junctions.csv', output_path='as_events_junctions_analysis.csv', df_junctions=None,
                          hadas_format=False, n=0, create_pdf=True):
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

        #df_junctions = df_junctions[df_junctions['gene_ensembl_id'].str.startswith('ENSG', na=False)]
        # check all needed columns are in df_junctions
        required_columns = ['gene_ensembl_id', 'start_position', 'end_position', 'cluster_name']
        for col in required_columns:
            if col not in df_junctions.columns:
                raise ValueError(f"Column '{col}' is required in df_junctions but not found.")

        if n > 0:
            df_t = pd.read_sql_query('select * from transcripts', self.con)
            gene_count = df_t.value_counts('gene_GeneID_id')
            df_t['gene_count'] = df_t['gene_GeneID_id'].map(gene_count)
            ids_n = df_t[df_t['gene_count'] == n].gene_GeneID_id.unique().tolist()
            df_junctions_n = df_junctions[df_junctions['gene_ensembl_id'].isin(ids_n)]
        else:
            df_junctions_n = df_junctions

        gene_ids = df_junctions_n.gene_ensembl_id.unique().tolist()
        df_transcripts = domas.get_genes_df_transcripts(self.con, gene_ids)
        invalid_ids = {'', 'nan', 'None'}
        all_transcript_ids = df_transcripts.transcript_ensembl_id.fillna(df_transcripts.transcript_refseq_id)
        valid_mask = all_transcript_ids.notna() & ~all_transcript_ids.isin(invalid_ids)
        canonical_transcript_ids = set(all_transcript_ids[valid_mask & (df_transcripts.canonical != 0)])
        transcript_ids = set(all_transcript_ids[valid_mask])
        df_domains = domas.get_transcript_domains_db(self.con, transcript_ids)

        # Local import avoids circular imports with aleternative_splicing module.
        from aleternative_splicing import get_exons_for_transcripts

        df_exons = get_exons_for_transcripts(self.con, transcript_ids)
        exon_lookup = build_exon_lookup(df_exons)
        domain_lookup = build_domain_lookup(df_domains)
        transcripts_by_gene = {gid: g for gid, g in df_transcripts.groupby('gene_ensembl_id')}
        empty_transcripts = df_transcripts.iloc[0:0]

        group_columns = ['specie', 'cluster_name'] if 'specie' in df_junctions_n.columns else ['cluster_name']
        cluster_groups = df_junctions_n.groupby(group_columns)
        total = len(cluster_groups)
        print(f"Analyzing {total} clusters...")
        last_time = time.perf_counter()

        results = []
        for count, (_, cluster_df) in enumerate(cluster_groups, start=1):
            cluster_name = cluster_df.cluster_name.iat[0]
            if count % 100 == 0:
                cur_time = time.perf_counter()
                self.logger.info(f"Analyzing cluster {count}/{total}: {cluster_name}. last 100 clusters duration: {cur_time - last_time:.2f}")
                last_time = cur_time

            gene_ensembl_id = cluster_df.gene_ensembl_id.iat[0]
            gene_symbol = cluster_df.gene_symbol.iat[0]
            event_type = cluster_df.event_type.iat[0] if 'event_type' in cluster_df.columns else None
            specie = cluster_df.specie.iat[0] if 'specie' in cluster_df.columns else None
            cluster_result = ClusterAnalysisResult(cluster_name, gene_ensembl_id, gene_symbol, as_event_type=event_type, specie=specie)
            cluster_result.junctions = list(zip(cluster_df['start_position'], cluster_df['end_position']))
            results.append(cluster_result)

            df_gene_transcripts = transcripts_by_gene.get(gene_ensembl_id, empty_transcripts)
            cluster_result.analyze_junction(df_gene_transcripts, canonical_transcript_ids, exon_lookup, domain_lookup)

        df_results_columns = ['cluster', 'gene_symbol', 'specie', 'event_type', 'canonical_transcript_id', 'transcript_id', 'domain_name',
                              'c_domain_length', 't_domain_length', 'c_domains_number', 't_domains_number']
        if results:
            df_all_results = pd.concat(
                (
                    cluster_result.get_results_df()
                    .rename(columns={'event': 'event_type'})
                    .assign(
                        cluster=cluster_result.cluster_name,
                        gene_symbol=cluster_result.gene_symbol,
                        canonical_transcript_id=cluster_result.canonical_transcript_id,
                        specie=cluster_result.specie,
                    )
                    for cluster_result in results
                ),
                ignore_index=True,
            ).reindex(columns=df_results_columns)
        else:
            df_all_results = pd.DataFrame(columns=df_results_columns)
        df_all_results.to_csv(output_path, index=False)

        if create_pdf:
            preloaded_gene_data = prepare_gene_data_bulk(
                self.con, [cluster_result.gene_symbol for cluster_result in results]
            )

            gene_visualizations = {}
            for count, cluster_result in enumerate(results, start=1):
                try:
                    df_cluster_junctions = pd.DataFrame(cluster_result.junctions, columns=['start', 'end'])
                    viz = gene_visualizations.get(cluster_result.gene_symbol)
                    if viz is None:
                        gene_symbol = cluster_result.gene_symbol
                        preloaded = preloaded_gene_data.get(gene_symbol.lower()) if isinstance(gene_symbol, str) else None
                        viz = self.gene_visualization_cls(self.con, gene_symbol, preloaded=preloaded)
                        gene_visualizations[gene_symbol] = viz
                    # name should have event_type and specie if available, to distinguish different clusters of the same gene and across species
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
                    print(f"Warning: Skipping PDF generation for {cluster_result.gene_symbol}, specie {cluster_result.specie}: {e}")

        return results
