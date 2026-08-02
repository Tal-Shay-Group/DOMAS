"""
Tests for junction_analisys.py, particularly error handling.
"""
import os
import sys

import pytest
import pandas as pd
import logging
from unittest.mock import patch

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.normpath(os.path.join(TESTS_DIR, '..', 'code'))
sys.path.insert(0, CODE_DIR)

from junction_analisys import (  # noqa: E402
    _init_worker, _process_cluster_chunk, ClusterAnalysisResult, JunctionsAnalysis,
    find_matching_junction_indices, get_aa_range, find_bp_range_for_domains,
    classify_domain_change, collapse_contained_domains, group_by_shared_names,
    select_most_like_canonical, DOMAIN_NAME_COLUMNS,
)

# Setup logging to capture messages
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def test_worker_handles_cluster_error_gracefully():
    """
    Test that _process_cluster_chunk continues processing even if one cluster fails.

    This verifies that:
    1. An error in one cluster doesn't crash the worker
    2. Other clusters in the chunk are still processed
    3. The error is logged
    4. Successful results are returned
    """
    # Create mock data
    worker_id = 0
    chunk_index = 0
    total_chunks = 1

    # Create 3 fake cluster tuples
    good_cluster_1 = (
        'group_key_1',
        pd.DataFrame({
            'gene_ensembl_id': ['ENSG00001'],
            'cluster_name': ['cluster_A'],
            'start_position': [1000],
            'end_position': [2000]
        })
    )

    bad_cluster = (
        'group_key_2',
        pd.DataFrame({
            'gene_ensembl_id': ['ENSG00002'],
            'cluster_name': ['cluster_B_will_fail'],
            'start_position': [3000],
            'end_position': [4000]
        })
    )

    good_cluster_2 = (
        'group_key_3',
        pd.DataFrame({
            'gene_ensembl_id': ['ENSG00003'],
            'cluster_name': ['cluster_C'],
            'start_position': [5000],
            'end_position': [6000]
        })
    )

    chunk = [good_cluster_1, bad_cluster, good_cluster_2]
    chunk_info = (worker_id, chunk_index, total_chunks, chunk)

    # Mock the supporting data structures
    df_exons = pd.DataFrame({
        'transcript_ensembl_id': ['ENST001', 'ENST002', 'ENST003'],
        'transcript_refseq_id': [None, None, None],
        'genomic_start_tx': [100, 3100, 5100],
        'genomic_end_tx': [1900, 3900, 5900],
        'order_in_transcript': [1, 1, 1],
        'abs_start_CDS': [100, 3100, 5100],
        'abs_end_CDS': [1900, 3900, 5900],
    })

    df_domains = pd.DataFrame({
        'transcript_ensembl_id_version': ['ENST001', 'ENST003'],
        'AA_start': [1, 5],
        'AA_end': [100, 200],
        'interpro': ['IPR001', 'IPR003']
    })

    df_transcripts = pd.DataFrame({
        'transcript_ensembl_id': ['ENST001', 'ENST002', 'ENST003'],
        'gene_ensembl_id': ['ENSG00001', 'ENSG00002', 'ENSG00003'],
        'canonical': [1, 0, 1]
    })

    canonical_transcript_ids = {'ENST001', 'ENST003'}
    transcripts_by_gene = {
        'ENSG00001': df_transcripts[df_transcripts['gene_ensembl_id'] == 'ENSG00001'],
        'ENSG00002': df_transcripts[df_transcripts['gene_ensembl_id'] == 'ENSG00002'],
        'ENSG00003': df_transcripts[df_transcripts['gene_ensembl_id'] == 'ENSG00003'],
    }
    gene_strand = {'ENSG00001': '+', 'ENSG00002': '+', 'ENSG00003': '+'}

    # Mock _analyze_single_cluster to fail on the bad cluster
    def mock_analyze(cluster_tuple, **kwargs):
        cluster_name = cluster_tuple[1].cluster_name.iat[0]
        if 'will_fail' in cluster_name:
            raise ValueError(f"Intentional test error for {cluster_name}")

        # Return a successful result
        result = ClusterAnalysisResult(
            cluster_name=cluster_name,
            gene_ensembl_id=cluster_tuple[1].gene_ensembl_id.iat[0],
            gene_symbol=f"GENE_{cluster_name}",
        )
        result.events = []
        return result

    # Run the worker with error handling
    with patch('junction_analisys._analyze_single_cluster', side_effect=mock_analyze):
        with patch('junction_analisys.build_exon_lookup') as mock_exon_lookup:
            with patch('junction_analisys.build_domain_lookup') as mock_domain_lookup:
                # Setup mocks
                mock_exon_lookup.return_value = lambda tid: df_exons
                mock_domain_lookup.return_value = lambda tid: df_domains

                # Populate the per-worker-process state _process_cluster_chunk reads,
                # the same way ProcessPoolExecutor's initializer does in a real run.
                _init_worker(df_exons, df_domains, canonical_transcript_ids, gene_strand,
                             transcripts_by_gene)

                # Call the worker
                results = _process_cluster_chunk(chunk_info)

    # Verify results
    print(f"\n✓ Worker processed {len(results)} clusters (should be 2, not 3)")
    assert len(results) == 2, f"Expected 2 successful results, got {len(results)}"

    # Verify the successful clusters are present
    cluster_names = [r.cluster_name for r in results]
    print(f"✓ Cluster names: {cluster_names}")
    assert 'cluster_A' in cluster_names, "cluster_A should be in results"
    assert 'cluster_C' in cluster_names, "cluster_C should be in results"
    assert 'cluster_B_will_fail' not in cluster_names, "cluster_B_will_fail should NOT be in results"

    print("✓ Worker continued processing after error")
    print("✓ Only successful clusters were returned")
    print("\n✅ Error handling test PASSED!")


# ---------------------------------------------------------------------------
# Regression tests for the minus-strand AA-window bugs (get_aa_range /
# find_bp_range_for_domains) - these were silently producing wrong/degenerate
# windows for minus-strand genes, causing real domain overlaps to be missed
# (e.g. CTSG, ARX reported "dropped domain"/"no_domains_in_region" instead of
# the correct comparison). See generate_reference_outputs.py / reference
# diffs from that fix for the real-data confirmation; these lock in the fix
# with small synthetic exon data so a future refactor can't reintroduce it.
# ---------------------------------------------------------------------------

def test_get_aa_range_minus_strand_does_not_collapse_to_single_codon():
    """
    On a minus strand transcript, the exon nearest the window's genomically
    LOWER bound is the LATER exon in CDS order (and vice versa). Using that
    exon's own abs_start_CDS as 'start' and the other exon's abs_end_CDS as
    'end' (the pre-fix behavior) picks two CDS-adjacent values that straddle
    a single codon instead of spanning the whole region.
    """
    # Mirrors the real CTSG case: exon order 2 (CDS 56-245) is genomically
    # lower than exon order 1 (CDS 1-55), because the gene is minus strand.
    exon_order_2 = pd.Series({'abs_start_CDS': 56, 'abs_end_CDS': 245})
    exon_order_1 = pd.Series({'abs_start_CDS': 1, 'abs_end_CDS': 55})

    min_aa, max_aa = get_aa_range(first_exon=exon_order_2, last_exon=exon_order_1)

    # The full CDS span is bases 1-245 -> AA 0-81. The pre-fix bug returned
    # (18, 18) - a single codon at the exon boundary.
    assert (min_aa, max_aa) == (0, 81)


def test_get_aa_range_with_explicit_bp_still_clamps_correctly():
    """Sanity check that the (unchanged) min_bp/max_bp clamping path still works."""
    exon = pd.Series({'abs_start_CDS': 100, 'abs_end_CDS': 200, 'genomic_start_tx': 5000})
    # bp 5050 is 50 bases into the exon -> CDS position 150 -> AA 50.
    min_aa, max_aa = get_aa_range(exon, exon, min_bp=5050, max_bp=5050)
    assert (min_aa, max_aa) == (50, 50)


def test_find_bp_range_for_domains_returns_genomically_ordered_pair():
    """
    On a minus strand transcript, the domain's lower CDS bound maps to a
    HIGHER genomic position. find_bp_range_for_domains must still return
    (min_bp, max_bp) with min_bp <= max_bp, since callers pool this pair into
    a min()/max() across transcripts - an unordered pair silently drops out
    of that pooling and the comparison window never widens to cover it.
    """
    df_exons = pd.DataFrame({
        'abs_start_CDS': [56, 637],
        'abs_end_CDS': [245, 810],
        'genomic_start_tx': [24575264, 24573560],
    })
    # A domain spanning AA 35-255 -> CDS bp 105-765, which falls inside the
    # two exons above (genomically reversed relative to CDS order).
    domains_in_region = pd.DataFrame({'AA_start': [35], 'AA_end': [255]})

    min_bp, max_bp = find_bp_range_for_domains(df_exons, domains_in_region)

    assert min_bp <= max_bp
    assert (min_bp, max_bp) == (24573688, 24575313)


def test_find_bp_range_for_domains_returns_none_when_no_real_domains():
    """A domains_in_region with only AA_start == 0 placeholder rows -> (None, None)."""
    df_exons = pd.DataFrame({
        'abs_start_CDS': [1], 'abs_end_CDS': [100], 'genomic_start_tx': [1000],
    })
    domains_in_region = pd.DataFrame({'AA_start': [0], 'AA_end': [0]})
    assert find_bp_range_for_domains(df_exons, domains_in_region) == (None, None)


# ---------------------------------------------------------------------------
# find_matching_junction_indices - strand handling and tolerance
# ---------------------------------------------------------------------------

def _two_exon_df(exon1, exon2):
    """Build a minimal 2-exon DataFrame: each exon is (order, start, end)."""
    orders, starts, ends = zip(exon1, exon2)
    return pd.DataFrame({
        'order_in_transcript': orders,
        'genomic_start_tx': starts,
        'genomic_end_tx': ends,
    })


def _three_exon_df(exon1, exon2, exon3):
    """Same as _two_exon_df, for layouts needing an exon outside the junction range."""
    orders, starts, ends = zip(exon1, exon2, exon3)
    return pd.DataFrame({
        'order_in_transcript': orders,
        'genomic_start_tx': starts,
        'genomic_end_tx': ends,
    })


def test_find_matching_junction_indices_plus_strand():
    df_exons = _two_exon_df((1, 100, 200), (2, 300, 400))
    # Intron between the two exons: starts right after exon1 ends, ends right
    # before exon2 starts.
    junctions = [(200, 300)]
    assert find_matching_junction_indices(df_exons, junctions, strand='+') == {0}


def test_find_matching_junction_indices_minus_strand():
    # Same exon layout, but on the minus strand the roles of start/end are
    # swapped per the function's own docstring.
    df_exons = _two_exon_df((1, 300, 400), (2, 100, 200))
    junctions = [(200, 300)]
    assert find_matching_junction_indices(df_exons, junctions, strand='-') == {0}


def test_find_matching_junction_indices_respects_1bp_tolerance():
    df_exons = _two_exon_df((1, 100, 200), (2, 300, 400))
    # 1bp off on both sides should still match.
    assert find_matching_junction_indices(df_exons, [(201, 299)], strand='+') == {0}
    # 2bp off should not.
    assert find_matching_junction_indices(df_exons, [(202, 298)], strand='+') == set()


def test_find_matching_junction_indices_non_adjacent_exons_do_not_match():
    """Exons that aren't adjacent in transcript order shouldn't match even if
    their genomic boundaries happen to line up with a junction."""
    df_exons = pd.DataFrame({
        'order_in_transcript': [1, 2, 3],
        'genomic_start_tx': [100, 300, 500],
        'genomic_end_tx': [200, 400, 600],
    })
    # This junction's boundaries match exon 1's end and exon 3's start, but
    # those exons are two apart (order 1 and 3), not adjacent.
    assert find_matching_junction_indices(df_exons, [(200, 500)], strand='+') == set()


# ---------------------------------------------------------------------------
# classify_domain_change - pure decision-matrix function
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('c_count,t_count,c_length,t_length,expected', [
    (0, 1, None, 50, 'added_domain'),
    (1, 0, 50, None, 'dropped domain'),
    (1, 1, 50, 50, 'same'),
    (1, 1, 50, 80, 'longer'),
    (1, 1, 80, 50, 'shorter'),
    (1, 2, 50, 90, 'split domain'),
    (2, 1, 90, 50, 'merged domain'),
    (2, 2, 100, 100, 'same_domains'),
    (2, 2, 100, 150, 'longer_domains'),
    (2, 2, 150, 100, 'shorter_domains'),
    (2, 3, 100, 100, 'increased_domain_number'),
    (3, 2, 100, 100, 'reduced_domain_number'),
])
def test_classify_domain_change(c_count, t_count, c_length, t_length, expected):
    assert classify_domain_change(c_count, t_count, c_length, t_length) == expected


# ---------------------------------------------------------------------------
# collapse_contained_domains
# ---------------------------------------------------------------------------

def _domain_row(aa_start, aa_end, **names):
    row = {'AA_start': aa_start, 'AA_end': aa_end}
    row.update({col: names.get(col) for col in DOMAIN_NAME_COLUMNS})
    return row


def test_collapse_contained_domains_merges_contained_domain():
    """A domain fully contained within another (within tolerance) is dropped,
    and its names are merged into the surviving (longer) domain."""
    df = pd.DataFrame([
        _domain_row(10, 100, interpro='IPR001'),
        _domain_row(12, 98, pfam='PF001'),  # contained within the first, within tolerance=2
    ])
    result = collapse_contained_domains(df, tolerance=2)
    assert len(result) == 1
    assert result.iloc[0]['AA_start'] == 10 and result.iloc[0]['AA_end'] == 100
    assert result.iloc[0]['interpro'] == 'IPR001'
    assert result.iloc[0]['pfam'] == 'PF001'


def test_collapse_contained_domains_keeps_non_overlapping_domains_separate():
    df = pd.DataFrame([
        _domain_row(10, 50, interpro='IPR001'),
        _domain_row(200, 250, interpro='IPR002'),
    ])
    result = collapse_contained_domains(df, tolerance=2)
    assert len(result) == 2


def test_collapse_contained_domains_single_row_is_unchanged():
    df = pd.DataFrame([_domain_row(10, 50, interpro='IPR001')])
    result = collapse_contained_domains(df, tolerance=2)
    assert len(result) == 1


# ---------------------------------------------------------------------------
# group_by_shared_names - transitive grouping
# ---------------------------------------------------------------------------

def test_group_by_shared_names_groups_transitively():
    """A shares a name with B, B shares a different name with C - A and C
    share nothing directly, but should still end up in the same group."""
    items = [
        ('A', {'x'}),
        ('B', {'x', 'y'}),
        ('C', {'y'}),
        ('D', {'z'}),  # unrelated, should stay in its own group
    ]
    groups = group_by_shared_names(items)
    groups_as_sets = [set(g) for g in groups]

    assert {'A', 'B', 'C'} in groups_as_sets
    assert {'D'} in groups_as_sets
    assert len(groups) == 2


def test_group_by_shared_names_no_overlap_keeps_singletons():
    items = [('A', {'x'}), ('B', {'y'}), ('C', {'z'})]
    groups = group_by_shared_names(items)
    assert sorted(groups, key=lambda g: g[0]) == [['A'], ['B'], ['C']]


# ---------------------------------------------------------------------------
# _load_junctions_data - input validation
#
# _load_junctions_data() only ever takes an already-loaded DataFrame now -
# reading junctions from a file (plain CSV, hadas-format Excel, IOE, ...) is
# alternative_splicing.py's job, not junction_analisys.py's.
# ---------------------------------------------------------------------------

def _bare_analysis():
    """A JunctionsAnalysis instance with no real DB connection - fine for the
    paths under test here, which never touch self.con."""
    return JunctionsAnalysis.__new__(JunctionsAnalysis)


def test_load_junctions_data_raises_on_missing_required_column():
    """Every missing required column is named at once, rather than only the first."""
    analysis = _bare_analysis()
    df = pd.DataFrame({'gene_ensembl_id': ['G1'], 'start_position': [1]})  # missing end_position, cluster_name
    with pytest.raises(ValueError, match="required") as excinfo:
        analysis._load_junctions_data(df)
    message = str(excinfo.value)
    assert 'end_position' in message and 'cluster_name' in message, message


def test_load_junctions_data_renames_cluster_to_cluster_name():
    analysis = _bare_analysis()
    df = pd.DataFrame({'gene_ensembl_id': ['G1'], 'start_position': [1],
                        'end_position': [2], 'cluster': ['c1']})
    result = analysis._load_junctions_data(df)
    assert 'cluster_name' in result.columns
    assert result['cluster_name'].iloc[0] == 'c1'


# ---------------------------------------------------------------------------
# select_most_like_canonical - the most-like-canonical tie-break rule
# ---------------------------------------------------------------------------

def _exon_df(*exons):
    """Build a DataFrame with only the columns select_most_like_canonical() needs.
    Each exon is a (start, end) pair."""
    starts, ends = zip(*exons)
    return pd.DataFrame({'genomic_start_tx': starts, 'genomic_end_tx': ends})


def test_select_most_like_canonical_prefers_longest_cds_among_qualifying():
    """Both candidates have exons outside the junction range [150, 350] that exactly
    match canonical's - so both qualify. Qualifying candidates are then separated by CDS
    length alone; how closely they match canonical *inside* the range is not considered,
    so A wins on CDS length despite B matching canonical's inside exon."""
    transcript_exons = {
        'CANON': _exon_df((1, 100), (200, 300), (400, 500)),
        'A': _exon_df((1, 100), (220, 280), (400, 500)),  # inside exon differs from canonical
        'B': _exon_df((1, 100), (200, 300), (400, 500)),  # inside exon also matches canonical
    }
    cds_length_by_transcript = {'A': 1000, 'B': 10}

    chosen = select_most_like_canonical(
        ['A', 'B'], 'CANON', transcript_exons, junctions=[(150, 350)],
        cds_length_by_transcript=cds_length_by_transcript,
    )
    assert chosen == 'A'


def test_select_most_like_canonical_excludes_candidate_with_different_outside_exons():
    """A candidate whose exons outside the junction range don't exactly match
    canonical's doesn't qualify, even if it matches everything inside the range and has
    the longer CDS - qualifying is a hard filter, not a preference traded off against
    CDS length."""
    transcript_exons = {
        'CANON': _exon_df((1, 100), (200, 300), (400, 500)),
        'DISQUALIFIED': _exon_df((1, 90), (200, 300), (400, 500)),  # outside exon (1,90) != (1,100)
        'QUALIFIES': _exon_df((1, 100), (220, 280), (400, 500)),  # outside exons match exactly
    }
    cds_length_by_transcript = {'DISQUALIFIED': 1000, 'QUALIFIES': 10}

    chosen = select_most_like_canonical(
        ['DISQUALIFIED', 'QUALIFIES'], 'CANON', transcript_exons, junctions=[(150, 350)],
        cds_length_by_transcript=cds_length_by_transcript,
    )
    assert chosen == 'QUALIFIES'


def test_select_most_like_canonical_returns_none_when_none_qualify():
    """If no candidate's outside-range exons exactly match canonical's, nothing is
    most like canonical and the flag goes unset - it does NOT fall back to the
    longest-CDS candidate, which is_longest_cds already marks on its own."""
    transcript_exons = {
        'CANON': _exon_df((1, 100), (200, 300), (400, 500)),
        'SHORT': _exon_df((1, 90), (200, 300), (400, 500)),
        'LONG': _exon_df((1, 80), (200, 300), (400, 500)),
    }
    cds_length_by_transcript = {'SHORT': 10, 'LONG': 1000}

    chosen = select_most_like_canonical(
        ['SHORT', 'LONG'], 'CANON', transcript_exons, junctions=[(150, 350)],
        cds_length_by_transcript=cds_length_by_transcript,
    )
    assert chosen is None


def test_select_most_like_canonical_tiebreak_is_deterministic():
    """Two qualifying candidates with an equal CDS length must resolve
    deterministically (not depend on dict/set iteration order)."""
    transcript_exons = {
        'CANON': _exon_df((1, 100), (200, 300), (400, 500)),
        'ENST_BBB': _exon_df((1, 100), (220, 280), (400, 500)),
        'ENST_AAA': _exon_df((1, 100), (230, 270), (400, 500)),
    }
    cds_length_by_transcript = {'ENST_BBB': 10, 'ENST_AAA': 10}

    chosen = select_most_like_canonical(
        ['ENST_BBB', 'ENST_AAA'], 'CANON', transcript_exons, junctions=[(150, 350)],
        cds_length_by_transcript=cds_length_by_transcript,
    )
    # Tied on CDS length -> max() on the (cds_length, tid) key picks the
    # lexicographically last id, the same tie-break is_longest_cds uses.
    assert chosen == 'ENST_BBB'


# ---------------------------------------------------------------------------
# analyze: every comparable transcript is compared, tagged with
# is_longest_cds / is_most_like_canonical
# ---------------------------------------------------------------------------

def test_analyze_junction_compares_all_and_tags_longest_cds_deterministically():
    """
    When two transcripts both qualify for comparison and have EQUAL CDS
    length, both are compared to canonical (none are skipped), but exactly
    one is tagged is_longest_cds - and which one must be deterministic (not
    dependent on dict/set iteration order or the per-process string hash seed).
    """
    df_gene_transcripts = pd.DataFrame({
        'transcript_ensembl_id': ['CANON', 'ENST_AAA', 'ENST_BBB'],
        'transcript_refseq_id': [None, None, None],
        'canonical': [1, 0, 0],
        'cds_start': [0, 0, 0],
        'cds_end': [1000, 1000, 1000],  # equal CDS length for both candidates
    })

    # Canonical exons -> matches junction 0 only.
    canonical_exons = _two_exon_df((1, 0, 100), (2, 200, 300))
    # ENST_AAA and ENST_BBB share the same exon shape -> both match junction 1 only.
    candidate_exons = _two_exon_df((1, 250, 300), (2, 400, 500))
    candidate_exons['abs_start_CDS'] = [1, 51]
    candidate_exons['abs_end_CDS'] = [50, 100]
    canonical_exons['abs_start_CDS'] = [1, 51]
    canonical_exons['abs_end_CDS'] = [50, 100]

    exons_by_id = {'CANON': canonical_exons, 'ENST_AAA': candidate_exons, 'ENST_BBB': candidate_exons}
    exon_lookup = lambda tid: exons_by_id[tid]

    empty_domains = pd.DataFrame(columns=['AA_start', 'AA_end'] + DOMAIN_NAME_COLUMNS)
    domain_lookup = lambda tid: empty_domains

    cluster_result = ClusterAnalysisResult('cluster_1', 'ENSG_TEST', 'TESTGENE')
    cluster_result.junctions = [(100, 200), (300, 400)]

    cluster_result.analyze(
        df_gene_transcripts, canonical_transcript_ids={'CANON'},
        exon_lookup=exon_lookup, domain_lookup=domain_lookup,
    )

    df_results = cluster_result.get_results_df()
    compared = sorted(df_results['transcript_id'].dropna().unique().tolist())
    longest_cds_tagged = sorted(df_results.loc[df_results['is_longest_cds'] == True, 'transcript_id'].unique().tolist())

    # Both candidates are compared - nothing is skipped.
    assert compared == ['ENST_AAA', 'ENST_BBB']
    # max() with key=(cds_length, tid) picks the lexicographically larger id on a tie.
    assert longest_cds_tagged == ['ENST_BBB']


def test_analyze_junction_leaves_most_like_canonical_unset_when_none_qualify():
    """A cluster whose sole comparable transcript fails the outside-exon gate is still
    compared and still tagged is_longest_cds, but is_most_like_canonical stays unset
    for every row - the flag is not backfilled with the longest-CDS transcript."""
    df_gene_transcripts = pd.DataFrame({
        'transcript_ensembl_id': ['CANON', 'ENST_OFF'],
        'transcript_refseq_id': [None, None],
        'canonical': [1, 0],
        'cds_start': [0, 0],
        'cds_end': [1000, 500],
    })

    # Junction range is [100, 400]. Canonical's only exon lying entirely outside it is
    # (500, 600); the candidate's is (700, 800) - so the outside-exon sets differ and
    # the candidate does not qualify.
    canonical_exons = _three_exon_df((1, 0, 100), (2, 200, 300), (3, 500, 600))
    candidate_exons = _three_exon_df((1, 250, 300), (2, 400, 450), (3, 700, 800))
    for df in (canonical_exons, candidate_exons):
        df['abs_start_CDS'] = [1, 51, 101]
        df['abs_end_CDS'] = [50, 100, 150]

    exons_by_id = {'CANON': canonical_exons, 'ENST_OFF': candidate_exons}
    exon_lookup = lambda tid: exons_by_id[tid]

    empty_domains = pd.DataFrame(columns=['AA_start', 'AA_end'] + DOMAIN_NAME_COLUMNS)
    domain_lookup = lambda tid: empty_domains

    cluster_result = ClusterAnalysisResult('cluster_1', 'ENSG_TEST', 'TESTGENE')
    cluster_result.junctions = [(100, 200), (300, 400)]

    cluster_result.analyze(
        df_gene_transcripts, canonical_transcript_ids={'CANON'},
        exon_lookup=exon_lookup, domain_lookup=domain_lookup,
    )

    df_results = cluster_result.get_results_df()
    # The transcript is still compared - the gate tags, it doesn't filter.
    assert df_results['transcript_id'].dropna().unique().tolist() == ['ENST_OFF']
    assert df_results['is_longest_cds'].iloc[0] == True
    assert not (df_results['is_most_like_canonical'] == True).any()


def test_analyze_junction_tags_sole_comparable_transcript_as_both():
    """A cluster with only one comparable transcript trivially tags it as both
    the longest-CDS pick and the most-like-canonical pick - there's no tie to
    break."""
    df_gene_transcripts = pd.DataFrame({
        'transcript_ensembl_id': ['CANON', 'ENST_ONLY'],
        'transcript_refseq_id': [None, None],
        'canonical': [1, 0],
        'cds_start': [0, 0],
        'cds_end': [1000, 500],
    })

    canonical_exons = _two_exon_df((1, 0, 100), (2, 200, 300))
    candidate_exons = _two_exon_df((1, 250, 300), (2, 400, 500))
    candidate_exons['abs_start_CDS'] = [1, 51]
    candidate_exons['abs_end_CDS'] = [50, 100]
    canonical_exons['abs_start_CDS'] = [1, 51]
    canonical_exons['abs_end_CDS'] = [50, 100]

    exons_by_id = {'CANON': canonical_exons, 'ENST_ONLY': candidate_exons}
    exon_lookup = lambda tid: exons_by_id[tid]

    empty_domains = pd.DataFrame(columns=['AA_start', 'AA_end'] + DOMAIN_NAME_COLUMNS)
    domain_lookup = lambda tid: empty_domains

    cluster_result = ClusterAnalysisResult('cluster_1', 'ENSG_TEST', 'TESTGENE')
    cluster_result.junctions = [(100, 200), (300, 400)]

    cluster_result.analyze(
        df_gene_transcripts, canonical_transcript_ids={'CANON'},
        exon_lookup=exon_lookup, domain_lookup=domain_lookup,
    )

    df_results = cluster_result.get_results_df()
    assert df_results['transcript_id'].dropna().unique().tolist() == ['ENST_ONLY']
    assert df_results['is_longest_cds'].iloc[0] == True
    assert df_results['is_most_like_canonical'].iloc[0] == True


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__, '-v']))


# ---------------------------------------------------------------------------
# utils.ioe2junctions - SUPPA .ioe parsing (previously untested)
# ---------------------------------------------------------------------------

def _write_ioe(tmp_path, event_ids):
    """Minimal SUPPA .ioe file: one row per event_id, seqname taken from the id."""
    lines = ['seqname\tgene_id\tevent_id\talternative_transcripts\ttotal_transcripts']
    for event_id in event_ids:
        gene = event_id.split(';')[0]
        seqname = event_id.split(':')[1]
        lines.append(f'{seqname}\t{gene}\t{event_id}\tENSMUST1\tENSMUST1,ENSMUST2')
    path = tmp_path / 'events_strict.ioe'
    path.write_text('\n'.join(lines) + '\n')
    return str(path)


def test_ioe2junctions_emits_ri_as_a_feature_pair(tmp_path):
    """A SUPPA retained intron is written as RI:<chr>:<s1>:<e1>-<s2>:<e2>:<strand>,
    holding a single junction because the retained isoform is defined by that
    junction's ABSENCE. It is emitted twice at the same coordinates - once as a
    junction, once as a retained-intron feature - so each isoform matches the one
    it actually carries. Emitted as a lone junction the event was unanalysable."""
    import utils
    ioe = _write_ioe(tmp_path, [
        'ENSMUSG1;SE:1:100-200:300-400:+',
        'ENSMUSG2;RI:1:500:600-700:800:+',
    ])
    df = utils.ioe2junctions(ioe)

    ri = df[df['event_type'] == 'RI']
    assert len(ri) == 2, f"RI must emit two features, got {len(ri)}"
    assert set(ri[utils.FEATURE_TYPE_COLUMN]) == {utils.FEATURE_JUNCTION,
                                                  utils.FEATURE_RETAINED_INTRON}
    # Both features describe the SAME interval - only the type distinguishes them.
    assert ri[['start_position', 'end_position']].drop_duplicates().shape[0] == 1
    assert list(ri['start_position'])[0] == 600 and list(ri['end_position'])[0] == 700

    # Non-RI events stay plain junctions.
    se = df[df['event_type'] == 'SE']
    assert set(se[utils.FEATURE_TYPE_COLUMN]) == {utils.FEATURE_JUNCTION}


def test_ioe2junctions_ri_only_file_is_analysable(tmp_path):
    """An all-RI file used to yield nothing. It must now produce a usable frame
    carrying the feature-type column."""
    import utils
    ioe = _write_ioe(tmp_path, ['ENSMUSG2;RI:1:500:600-700:800:+'])
    df = utils.ioe2junctions(ioe)

    assert len(df) == 2
    assert utils.FEATURE_TYPE_COLUMN in df.columns
    assert df['cluster_name'].nunique() == 1, "both features belong to one event"


def test_ioe2junctions_keeps_both_forms_of_alt_splice_site_events(tmp_path):
    """A3/A5 events carry two explicit junctions, so - unlike the rMATS path, which
    names exons by transcript role - they need no strand handling and must not
    collapse to one junction."""
    import utils
    ioe = _write_ioe(tmp_path, [
        'ENSMUSG3;A3:1:100-200:100-250:-',
        'ENSMUSG4;A5:1:300-500:400-500:-',
    ])
    df = utils.ioe2junctions(ioe)

    for cluster, rows in df.groupby('cluster_name'):
        distinct = rows[['start_position', 'end_position']].drop_duplicates()
        assert len(distinct) == 2, f"{cluster} collapsed to {len(distinct)} junction(s)"


# ---------------------------------------------------------------------------
# utils._rmats_event_junctions - strand handling for the alt-splice-site types
# ---------------------------------------------------------------------------

def _a5ss_row(strand, long_start, long_end, short_start, short_end, flank_start, flank_end):
    return pd.Series({
        'strand': strand,
        'longExonStart_0base': long_start, 'longExonEnd': long_end,
        'shortES': short_start, 'shortEE': short_end,
        'flankingES': flank_start, 'flankingEE': flank_end,
    })


def test_rmats_a5ss_minus_strand_yields_two_distinct_junctions():
    """Real minus-strand A5SS geometry: the long and short forms of the exon share
    their genomic END and differ at their genomic START, because the alternative
    donor of a minus-strand transcript is the genomically lower boundary. Reading
    the plus-strand boundary here made both forms produce the same junction,
    collapsing the event to one junction - which can never yield a comparable
    transcript, so every minus-strand A5SS event was silently unanalysable."""
    import utils
    # Coordinates taken from a real fixture row (ENSG on chr7, minus strand).
    row = _a5ss_row('-', 45283351, 45283527, 45283382, 45283527, 45282348, 45282459)

    junctions = utils._rmats_event_junctions('A5SS', row)

    assert len(set(junctions)) == 2, f"expected two distinct junctions, got {junctions}"
    # Both junctions run from the flanking exon's upper boundary to the alternative
    # donor, which is what differs between the two forms.
    assert set(junctions) == {(45282459, 45283351), (45282459, 45283382)}


def test_rmats_a5ss_plus_strand_unchanged():
    """The plus-strand construction must be untouched by the minus-strand fix."""
    import utils
    row = _a5ss_row('+', 1000, 1500, 1000, 1400, 2000, 2100)

    junctions = utils._rmats_event_junctions('A5SS', row)

    assert set(junctions) == {(1500, 2000), (1400, 2000)}


def test_rmats_a3ss_strands_are_mirror_images():
    """A3SS is the mirror of A5SS: the varying boundary is the genomic START on the
    plus strand and the genomic END on the minus strand. Both must give two
    distinct junctions."""
    import utils
    plus = _a5ss_row('+', 1000, 1500, 1100, 1500, 500, 600)
    minus = _a5ss_row('-', 1000, 1500, 1000, 1400, 2000, 2100)

    plus_junctions = utils._rmats_event_junctions('A3SS', plus)
    minus_junctions = utils._rmats_event_junctions('A3SS', minus)

    assert set(plus_junctions) == {(600, 1000), (600, 1100)}
    assert set(minus_junctions) == {(1500, 2000), (1400, 2000)}


def test_rmats_mxe_does_not_require_alt_splice_site_columns():
    """MXE rows carry none of the longExon*/short*/flanking* columns. Building those
    coordinate pairs unconditionally raised KeyError, which rmats2junctions() swallows
    via `except (KeyError, ValueError, TypeError): continue` - silently dropping every
    MXE event. Guard against that regression."""
    import utils
    row = pd.Series({
        'strand': '-',
        'upstreamEE': 100, 'downstreamES': 900,
        '1stExonStart_0base': 200, '1stExonEnd': 300,
        '2ndExonStart_0base': 400, '2ndExonEnd': 500,
    })

    junctions = utils._rmats_event_junctions('MXE', row)

    assert junctions == [(100, 200), (300, 900), (100, 400), (500, 900)]


# ---------------------------------------------------------------------------
# Retained-intron event features (FEATURE_RETAINED_INTRON)
# ---------------------------------------------------------------------------

def _spliced_exons():
    """Two exons abutting the intron (1100, 1300) - the spliced isoform."""
    return _two_exon_df((1, 1000, 1100), (2, 1300, 1400))


def _retained_exons():
    """One exon containing the intron (1100, 1300) - the retained isoform."""
    return pd.DataFrame({'order_in_transcript': [1],
                         'genomic_start_tx': [1000], 'genomic_end_tx': [1400]})


INTRON = (1100, 1300)


def test_retained_intron_feature_matrix():
    """The four-way matrix that makes an RI event analysable: each isoform must
    match its own feature type and only that one."""
    from junction_analisys import FEATURE_JUNCTION, FEATURE_RETAINED_INTRON

    cases = {
        ('spliced', FEATURE_JUNCTION): {0},
        ('spliced', FEATURE_RETAINED_INTRON): set(),
        ('retained', FEATURE_JUNCTION): set(),
        ('retained', FEATURE_RETAINED_INTRON): {0},
    }
    exons = {'spliced': _spliced_exons(), 'retained': _retained_exons()}
    for (isoform, feature_type), expected in cases.items():
        got = find_matching_junction_indices(
            exons[isoform], [INTRON], feature_types=[feature_type])
        assert got == expected, f"{isoform} vs {feature_type}: expected {expected}, got {got}"


def test_retained_intron_event_makes_both_isoforms_comparable():
    """An RI event emits both features at the SAME coordinates. Comparability is
    `matched - canonical_matched`, and it must be non-empty in BOTH directions:
    whichever isoform is canonical, the other one differs from it within the event."""
    from junction_analisys import FEATURE_JUNCTION, FEATURE_RETAINED_INTRON

    features = [INTRON, INTRON]
    types = [FEATURE_JUNCTION, FEATURE_RETAINED_INTRON]
    spliced = find_matching_junction_indices(_spliced_exons(), features, feature_types=types)
    retained = find_matching_junction_indices(_retained_exons(), features, feature_types=types)

    assert spliced == {0}, spliced
    assert retained == {1}, retained
    assert retained - spliced == {1}, "retained isoform must be comparable to a spliced canonical"
    assert spliced - retained == {0}, "spliced isoform must be comparable to a retained canonical"


def test_retained_intron_matching_ignores_strand():
    """Containment is a genomic test, so a retained intron matches identically on
    both strands - unlike the junction predicate, whose boundary roles swap."""
    from junction_analisys import FEATURE_RETAINED_INTRON

    for strand in ('+', '-'):
        got = find_matching_junction_indices(
            _retained_exons(), [INTRON], strand=strand,
            feature_types=[FEATURE_RETAINED_INTRON])
        assert got == {0}, f"strand {strand}: {got}"


def test_retained_intron_requires_a_single_containing_exon():
    """An exon that only partly overlaps the intron does not count - otherwise any
    transcript whose exon merely reaches into the intron would look retained."""
    from junction_analisys import FEATURE_RETAINED_INTRON

    partial = pd.DataFrame({'order_in_transcript': [1, 2],
                            'genomic_start_tx': [1000, 1250],
                            'genomic_end_tx': [1200, 1400]})  # neither exon spans 1100-1300
    got = find_matching_junction_indices(partial, [INTRON],
                                         feature_types=[FEATURE_RETAINED_INTRON])
    assert got == set(), got


def test_feature_types_none_is_backward_compatible():
    """Omitting feature_types must behave exactly as before the parameter existed -
    every feature treated as a junction. Junctions frames on disk have no such
    column, so this is the path all existing input takes."""
    features = [INTRON]
    assert find_matching_junction_indices(_spliced_exons(), features) == {0}
    assert find_matching_junction_indices(_retained_exons(), features) == set()


def test_analyze_treats_retained_intron_transcript_as_comparable():
    """End-to-end: an RI event's two features make the retained-intron transcript
    comparable to a spliced canonical. Before feature types existed this cluster
    produced no comparable transcript at all - the spliced form satisfied both
    identically-positioned features, so nothing held a feature canonical lacked."""
    from junction_analisys import FEATURE_JUNCTION, FEATURE_RETAINED_INTRON

    df_gene_transcripts = pd.DataFrame({
        'transcript_ensembl_id': ['CANON', 'ENST_RETAINED'],
        'transcript_refseq_id': [None, None],
        'canonical': [1, 0],
        'cds_start': [0, 0],
        'cds_end': [1000, 900],
    })
    # canonical splices the intron (1100, 1300); the other retains it in one exon.
    canonical_exons = _two_exon_df((1, 1000, 1100), (2, 1300, 1400))
    retained_exons = pd.DataFrame({'order_in_transcript': [1],
                                   'genomic_start_tx': [1000], 'genomic_end_tx': [1400]})
    for df in (canonical_exons, retained_exons):
        df['abs_start_CDS'] = [1] * len(df)
        df['abs_end_CDS'] = [100] * len(df)

    exons_by_id = {'CANON': canonical_exons, 'ENST_RETAINED': retained_exons}
    empty_domains = pd.DataFrame(columns=['AA_start', 'AA_end'] + DOMAIN_NAME_COLUMNS)

    cluster_result = ClusterAnalysisResult('cluster_ri', 'ENSG_TEST', 'TESTGENE')
    cluster_result.junctions = [(1100, 1300), (1100, 1300)]
    cluster_result.feature_types = [FEATURE_JUNCTION, FEATURE_RETAINED_INTRON]

    cluster_result.analyze(
        df_gene_transcripts, canonical_transcript_ids={'CANON'},
        exon_lookup=lambda tid: exons_by_id[tid], domain_lookup=lambda tid: empty_domains,
    )

    df_results = cluster_result.get_results_df()
    skipped = {'transcript_doesnt_have_features', 'no_unique_features',
               'no_canonical_features', 'feature_not_mapped', 'gene_not_in_db'}
    comparable = set(df_results.loc[~df_results['event'].isin(skipped), 'transcript_id'].dropna())
    assert comparable == {'ENST_RETAINED'}, f"expected the retained transcript, got {comparable}"


# ---------------------------------------------------------------------------
# utils.voila2junctions - MAJIQ intron retention
# ---------------------------------------------------------------------------

def _write_voila(tmp_path, junctions_coords, ir_coords):
    header = ('#Gene Name\tGene ID\tLSV ID\tchr\tstrand\tA5SS\tA3SS\tES\t'
              'Junctions coords\tIR coords')
    row = ('BRD9\tENSG00000028310\tENSG00000028310:s:1\tchr5\t-\tFalse\tFalse\tTrue\t'
           f'{junctions_coords}\t{ir_coords}')
    path = tmp_path / 'voila.tsv'
    path.write_text(header + '\n' + row + '\n')
    return str(path)


def test_voila_emits_retained_intron_once_as_a_span():
    """MAJIQ quantifies intron retention as one of the LSV's edges, so the retained
    intron is listed in 'Junctions coords' as well, with 'IR coords' naming which
    edge it is. That edge is not a splice junction - the real junction for the same
    intron is listed separately, at (a-1, b+1) here - so the interval must appear
    ONCE, as a containment feature. Emitting it as a junction too asked the opposite
    question of it and duplicated the real junction, which the 1bp tolerance already
    matches."""
    import utils, tempfile, pathlib
    with tempfile.TemporaryDirectory() as d:
        tsv = _write_voila(pathlib.Path(d),
                           '883450-883938;881182-883938;883451-883937',
                           '883451-883937')
        df = utils.voila2junctions(tsv)

    ir_rows = df[(df.start_position == 883451) & (df.end_position == 883937)]
    assert len(ir_rows) == 1, f"the IR interval must appear once, got {len(ir_rows)}"
    assert ir_rows[utils.FEATURE_TYPE_COLUMN].iat[0] == utils.FEATURE_RETAINED_INTRON

    # The real splice junction for that intron is still present, as a junction.
    real = df[(df.start_position == 883450) & (df.end_position == 883938)]
    assert len(real) == 1 and real[utils.FEATURE_TYPE_COLUMN].iat[0] == utils.FEATURE_JUNCTION

    # No interval may carry both types.
    assert (df.groupby(['start_position', 'end_position'])[utils.FEATURE_TYPE_COLUMN]
            .nunique() == 1).all()


def test_voila_lsv_without_intron_retention_is_all_junctions():
    """An LSV with no IR coords must be unaffected by the retained-intron path."""
    import utils, tempfile, pathlib
    with tempfile.TemporaryDirectory() as d:
        tsv = _write_voila(pathlib.Path(d), '100-200;300-400', '')
        df = utils.voila2junctions(tsv)

    assert len(df) == 2
    assert set(df[utils.FEATURE_TYPE_COLUMN]) == {utils.FEATURE_JUNCTION}


# ---------------------------------------------------------------------------
# Selection step 1: protein-coding, and the CDS-length metric
# ---------------------------------------------------------------------------

def _selection_fixture(protein_ids, cds_offsets):
    """Two comparable candidates plus canonical. `protein_ids` maps transcript id to
    its protein_ensembl_id ('' for non-coding); `cds_offsets` maps it to the exons'
    largest abs_end_CDS, i.e. the coding length the selection sees."""
    ids = ['CANON', 'ENST_A', 'ENST_B']
    df_gene_transcripts = pd.DataFrame({
        'transcript_ensembl_id': ids,
        'transcript_refseq_id': [None] * 3,
        'canonical': [1, 0, 0],
        'cds_start': [0, 0, 0],
        'cds_end': [1000, 1000, 1000],
        'protein_ensembl_id': [protein_ids.get(i, 'ENSP1') for i in ids],
        'protein_refseq_id': [None] * 3,
    })
    canonical_exons = _two_exon_df((1, 0, 100), (2, 200, 300))
    canonical_exons['abs_start_CDS'] = [1, 51]
    canonical_exons['abs_end_CDS'] = [50, cds_offsets.get('CANON', 100)]

    exons_by_id = {'CANON': canonical_exons}
    for tid in ('ENST_A', 'ENST_B'):
        # both candidates match junction 1 only -> both comparable
        exons = _two_exon_df((1, 250, 300), (2, 400, 500))
        exons['abs_start_CDS'] = [1, 51]
        exons['abs_end_CDS'] = [50, cds_offsets[tid]]
        exons_by_id[tid] = exons
    return df_gene_transcripts, exons_by_id


def _run_selection(df_gene_transcripts, exons_by_id):
    empty_domains = pd.DataFrame(columns=['AA_start', 'AA_end'] + DOMAIN_NAME_COLUMNS)
    cluster_result = ClusterAnalysisResult('cluster_1', 'ENSG_TEST', 'TESTGENE')
    cluster_result.junctions = [(100, 200), (300, 400)]
    cluster_result.analyze(
        df_gene_transcripts, canonical_transcript_ids={'CANON'},
        exon_lookup=lambda tid: exons_by_id[tid], domain_lookup=lambda tid: empty_domains,
    )
    df = cluster_result.get_results_df()
    def tagged(col):
        ids = set(df.loc[df[col] == True, 'transcript_id'].dropna())
        return ids.pop() if ids else None
    return tagged('is_longest_cds'), tagged('is_most_like_canonical')


def _run_without_canonical(df_gene_transcripts, exons_by_id):
    """Same fixture, but with no transcript flagged canonical - so the longest-CDS
    fallback has to choose one. Returns the ClusterAnalysisResult."""
    empty_domains = pd.DataFrame(columns=['AA_start', 'AA_end'] + DOMAIN_NAME_COLUMNS)
    cluster_result = ClusterAnalysisResult('cluster_1', 'ENSG_TEST', 'TESTGENE')
    cluster_result.junctions = [(100, 200), (300, 400)]
    cluster_result.analyze(
        df_gene_transcripts, canonical_transcript_ids=set(),
        exon_lookup=lambda tid: exons_by_id[tid], domain_lookup=lambda tid: empty_domains,
    )
    return cluster_result


def test_no_canonical_falls_back_to_longest_cds():
    """With no canonical transcript flagged for the gene, the longest-CDS transcript
    stands in, so the cluster is analyzed instead of dropped."""
    df, exons = _selection_fixture(protein_ids={},
                                   cds_offsets={'ENST_A': 900, 'ENST_B': 300})
    cluster_result = _run_without_canonical(df, exons)
    assert cluster_result.canonical_transcript_id == 'ENST_A'
    assert 'no_canonical_transcript' not in {e[0] for e in cluster_result.events}


def test_no_canonical_fallback_prefers_protein_coding():
    """The fallback follows the same priority as the longest-CDS tag: a non-coding
    transcript does not become the stand-in canonical just for being longer. It has
    no domains, so every canonical domain would trivially read as dropped."""
    df, exons = _selection_fixture(protein_ids={'ENST_A': ''},        # A non-coding
                                   cds_offsets={'ENST_A': 900, 'ENST_B': 300})
    cluster_result = _run_without_canonical(df, exons)
    assert cluster_result.canonical_transcript_id == 'ENST_B', \
        f"non-coding ENST_A must not become canonical on length, got {cluster_result.canonical_transcript_id}"


def test_canonical_flag_wins_over_longer_cds():
    """The fallback is only for genes with no canonical at all - where one IS
    flagged, it stays canonical even though another transcript has a longer CDS."""
    df, exons = _selection_fixture(protein_ids={},
                                   cds_offsets={'CANON': 100, 'ENST_A': 900, 'ENST_B': 300})
    empty_domains = pd.DataFrame(columns=['AA_start', 'AA_end'] + DOMAIN_NAME_COLUMNS)
    cluster_result = ClusterAnalysisResult('cluster_1', 'ENSG_TEST', 'TESTGENE')
    cluster_result.junctions = [(100, 200), (300, 400)]
    cluster_result.analyze(
        df, canonical_transcript_ids={'CANON'},
        exon_lookup=lambda tid: exons[tid], domain_lookup=lambda tid: empty_domains,
    )
    assert cluster_result.canonical_transcript_id == 'CANON'


def test_selection_prefers_protein_coding_over_longer_cds():
    """Step 1 of the priority: a non-coding candidate is not selected even when it
    has the longer CDS. A transcript with no protein has no domains, so comparing it
    trivially reports every canonical domain as dropped."""
    df, exons = _selection_fixture(protein_ids={'ENST_A': ''},          # A non-coding
                                   cds_offsets={'ENST_A': 900, 'ENST_B': 300})
    longest, most_like = _run_selection(df, exons)
    assert longest == 'ENST_B', f"non-coding ENST_A must not win on length, got {longest}"
    assert most_like == 'ENST_B', most_like


def test_selection_falls_back_when_no_candidate_is_coding():
    """Where no candidate is coding they tie on step 1 and selection falls through,
    so the cluster still resolves to one transcript rather than none."""
    df, exons = _selection_fixture(protein_ids={'ENST_A': '', 'ENST_B': ''},
                                   cds_offsets={'ENST_A': 900, 'ENST_B': 300})
    longest, most_like = _run_selection(df, exons)
    assert longest == 'ENST_A', f"expected the longest of the non-coding pair, got {longest}"
    assert most_like is not None


def test_selection_without_protein_columns_does_not_filter():
    """Absent protein columns mean "unknown", not "non-coding" - a caller building
    df_gene_transcripts by hand must not silently lose every candidate."""
    df, exons = _selection_fixture(protein_ids={}, cds_offsets={'ENST_A': 900, 'ENST_B': 300})
    df = df.drop(columns=['protein_ensembl_id', 'protein_refseq_id'])
    longest, _ = _run_selection(df, exons)
    assert longest == 'ENST_A', longest


def test_cds_length_is_coding_length_not_genomic_span():
    """The longest-CDS tag ranks by coding length (largest abs_end_CDS), not by
    cds_end - cds_start, which is a genomic span counting the introns between the
    first and last coding exon - a median of ~10x the coding length. Here both
    candidates carry the SAME cds_start/cds_end, so a genomic-span metric could not
    separate them at all; only the exonic offsets can."""
    df, exons = _selection_fixture(protein_ids={}, cds_offsets={'ENST_A': 120, 'ENST_B': 800})
    assert (df.cds_end - df.cds_start).nunique() == 1, "fixture: genomic spans are equal"
    longest, _ = _run_selection(df, exons)
    assert longest == 'ENST_B', f"expected the longer coding length, got {longest}"


# ---------------------------------------------------------------------------
# Gene identity: gene_ensembl_id OR gene_GeneID_id
# ---------------------------------------------------------------------------

def test_combined_gene_ids_falls_back_to_geneid():
    """13% of DoChaP genes carry no gene_ensembl_id, only a gene_GeneID_id. The key
    a lookup uses must fall back to it, mirroring the transcript_ensembl_id /
    transcript_refseq_id fallback already used for transcripts."""
    import utils
    df = pd.DataFrame({'gene_ensembl_id': ['ENSG1', None, 'ENSG3'],
                       'gene_GeneID_id': ['111', '222', None]})
    assert list(utils.combined_gene_ids(df)) == ['ENSG1', '222', 'ENSG3']


def test_combined_gene_ids_without_geneid_column():
    """A frame carrying only gene_ensembl_id must pass through unchanged."""
    import utils
    df = pd.DataFrame({'gene_ensembl_id': ['ENSG1', 'ENSG2']})
    assert list(utils.combined_gene_ids(df)) == ['ENSG1', 'ENSG2']


def test_gene_id_clause_matches_either_column():
    import utils
    clause, params = utils.gene_id_clause(['ENSG1', '222'])
    for column in utils.GENE_ID_COLUMNS:
        assert f'{column} IN (?,?)' in clause, clause
    assert ' OR ' in clause
    # ids repeated once per column, so the same params bind both IN lists
    assert params == ['ENSG1', '222', 'ENSG1', '222']


def test_transcripts_grouped_by_combined_gene_id_keeps_geneid_only_genes():
    """The regression this fixes: grouping on gene_ensembl_id alone silently drops
    every gene without one, because groupby() skips NaN keys - so the cluster
    reported gene_not_in_db even though the gene is present in the database."""
    import utils
    df_transcripts = pd.DataFrame({
        'transcript_ensembl_id': ['ENST1', None, 'ENST3'],
        'transcript_refseq_id': [None, 'NM_2', None],
        'gene_ensembl_id': ['ENSG_A', None, None],
        'gene_GeneID_id': ['111', '222', '222'],
    })
    dropped = {gid: g for gid, g in df_transcripts.groupby('gene_ensembl_id')}
    kept = {gid: g for gid, g in df_transcripts.groupby(utils.combined_gene_ids(df_transcripts))}

    assert '222' not in dropped, "precondition: the old grouping loses the GeneID-only gene"
    assert set(kept) == {'ENSG_A', '222'}
    assert len(kept['222']) == 2, "both transcripts of the GeneID-only gene are reachable"


# ---------------------------------------------------------------------------
# LeafCutter gene annotation: one cluster per named gene, and no_gene_specified
# ---------------------------------------------------------------------------

def test_leafcutter_gene_symbols_parsing():
    """A significance row can name several genes, or none. Missing must give [] -
    it previously went through str(genes).split(','), which turns NaN into the
    literal string 'nan' and carried that into the output's gene_symbol column."""
    from alternative_splicing import _leafcutter_gene_symbols
    assert _leafcutter_gene_symbols('CDK11B,CDK11A,AL031282.2') == ['CDK11B', 'CDK11A', 'AL031282.2']
    assert _leafcutter_gene_symbols('PFKP') == ['PFKP']
    assert _leafcutter_gene_symbols(float('nan')) == []
    assert _leafcutter_gene_symbols(None) == []
    assert _leafcutter_gene_symbols('') == []
    assert _leafcutter_gene_symbols('.') == []
    # whitespace trimmed, duplicates collapsed, order preserved
    assert _leafcutter_gene_symbols(' A , B ,A') == ['A', 'B']


def test_missing_gene_id_detection():
    """A reader may leave an unnamed gene as None, NaN or a placeholder string."""
    from junction_analisys import _is_missing_gene_id
    for value in (None, float('nan'), '', '  ', 'nan', 'None', 'NA', '.'):
        assert _is_missing_gene_id(value), value
    for value in ('ENSG00000123456', '79400', 'CDK11B'):
        assert not _is_missing_gene_id(value), value


def test_cluster_with_no_gene_reports_no_gene_specified():
    """An event naming no gene is distinct from one naming a gene that isn't in the
    database: no lookup was possible, so gene_not_in_db claimed a failure that never
    happened. LeafCutter builds clusters annotation-free, so this is an expected
    outcome and needs its own label to keep the coverage audit honest."""
    cluster_result = ClusterAnalysisResult('chr1:clu_1', float('nan'), None, specie='H_sapiens')
    cluster_result.junctions = [(100, 200)]
    cluster_result.analyze(
        df_gene_transcripts=None, canonical_transcript_ids=set(),
        exon_lookup=lambda t: pd.DataFrame(), domain_lookup=lambda t: pd.DataFrame(),
    )
    events = [e[0] for e in cluster_result.events]
    assert events == ['no_gene_specified'], events


def test_named_but_absent_gene_still_reports_gene_not_in_db():
    """The new label must not swallow the old one: a gene that WAS named and is
    genuinely missing from the database still reports gene_not_in_db."""
    cluster_result = ClusterAnalysisResult('chr1:clu_2', 'ENSG99999999', 'FAKE', specie='H_sapiens')
    cluster_result.junctions = [(100, 200)]
    cluster_result.analyze(
        df_gene_transcripts=None, canonical_transcript_ids=set(),
        exon_lookup=lambda t: pd.DataFrame(), domain_lookup=lambda t: pd.DataFrame(),
    )
    events = [e[0] for e in cluster_result.events]
    assert events == ['gene_not_in_db'], events


def test_named_gene_with_unresolved_symbol_is_not_no_gene_specified():
    """A gene that WAS named but whose symbol did not resolve to a database id - a
    lncRNA clone id such as AL451164.2, absent from DoChaP - arrives with no gene
    id, exactly like an unannotated cluster. It is not the same thing: a lookup was
    possible and failed, so it must stay gene_not_in_db. Keying only on the missing
    id mislabels 1,582 clusters of the leafcutter fixture as no_gene_specified."""
    cluster_result = ClusterAnalysisResult('chr10:clu_1:AL451164.2', float('nan'),
                                           'AL451164.2', specie='H_sapiens')
    cluster_result.junctions = [(100, 200)]
    cluster_result.analyze(
        df_gene_transcripts=None, canonical_transcript_ids=set(),
        exon_lookup=lambda t: pd.DataFrame(), domain_lookup=lambda t: pd.DataFrame(),
    )
    events = [e[0] for e in cluster_result.events]
    assert events == ['gene_not_in_db'], events


# ---------------------------------------------------------------------------
# The junctions-frame contract (utils.normalize_junctions_frame)
# ---------------------------------------------------------------------------

def test_normalize_junctions_frame_fills_optional_columns():
    """Every optional column is present with its declared default afterwards, so
    downstream code can read it directly instead of carrying its own fallback."""
    import utils
    df = pd.DataFrame({'gene_ensembl_id': ['ENSG1'], 'start_position': [1],
                       'end_position': [2], 'cluster_name': ['c1']})
    out = utils.normalize_junctions_frame(df)
    for column in utils.REQUIRED_JUNCTION_COLUMNS:
        assert column in out.columns
    assert out['feature_type'].iat[0] == utils.FEATURE_JUNCTION
    assert out['event_type'].iat[0] is None
    assert 'gene_symbol' in out.columns and 'junction_name' in out.columns


def test_normalize_junctions_frame_does_not_overwrite_supplied_values():
    """A reader that already sets a column keeps its value."""
    import utils
    df = pd.DataFrame({'gene_ensembl_id': ['ENSG1'], 'start_position': [1],
                       'end_position': [2], 'cluster_name': ['c1'],
                       'specie': ['mouse'], 'event_type': ['SE'],
                       'feature_type': [utils.FEATURE_RETAINED_INTRON]})
    out = utils.normalize_junctions_frame(df)
    assert out['specie'].iat[0] == 'mouse', "a human ENSG id must not override a stated specie"
    assert out['event_type'].iat[0] == 'SE'
    assert out['feature_type'].iat[0] == utils.FEATURE_RETAINED_INTRON


def test_specie_derived_for_the_tool_formats():
    """rMATS, MAJIQ and SUPPA embed an Ensembl gene id but no species, so their
    frames had no specie column at all - and clusters are grouped by
    (specie, cluster_name), so a multi-species run merged same-named clusters."""
    import utils
    df = pd.DataFrame({'gene_ensembl_id': ['ENSG1', 'ENSMUSG2', 'ENSRNOG3'],
                       'start_position': [1, 1, 1], 'end_position': [2, 2, 2],
                       'cluster_name': ['c', 'c', 'c']})
    out = utils.normalize_junctions_frame(df)
    assert list(out['specie']) == ['human', 'mouse', 'rat']


def test_retired_species_are_rejected_not_silently_analysed():
    """zebrafish and frog were dropped (no MANE/RefSeq Select, so no canonical to
    compare against). Their ids must still be derivable: an underived species is
    indistinguishable from a GeneID-keyed gene, which normalize lets pass, so
    dropping the prefixes would turn a hard stop into a silent wrong-species run."""
    import utils
    assert utils.specie_from_gene_id('ENSDARG3') == 'zebrafish'
    assert utils.specie_from_gene_id('ENSXETG4') == 'frog'

    def frame(gene_id):
        return pd.DataFrame({'gene_ensembl_id': [gene_id], 'start_position': [1],
                             'end_position': [2], 'cluster_name': ['c']})

    # derived from the ids, with no -specie stated
    with pytest.raises(ValueError, match='no longer supported'):
        utils.normalize_junctions_frame(frame('ENSDARG3'))
    # asked for explicitly
    with pytest.raises(ValueError, match='no longer supported'):
        utils.normalize_junctions_frame(frame('ENSG1'), specie='frog')
    # stated as a supported species, but the ids say otherwise
    with pytest.raises(ValueError, match='does not match'):
        utils.normalize_junctions_frame(frame('ENSDARG3'), specie='human')


def test_cluster_grouping_keeps_rows_whose_specie_cannot_be_derived():
    """specie is derived from the Ensembl gene id, so a GeneID-keyed or otherwise
    non-Ensembl input leaves it None. groupby() drops NaN keys by default, which
    would silently discard those clusters - the same trap that lost every gene
    without a gene_ensembl_id."""
    import utils
    analysis = _bare_analysis()
    df = utils.normalize_junctions_frame(pd.DataFrame({
        'gene_ensembl_id': ['79400', '79400'],   # a GeneID, not an Ensembl id
        'start_position': [1, 3], 'end_position': [2, 4],
        'cluster_name': ['clu_a', 'clu_a'],
    }))
    assert df['specie'].isna().all(), "precondition: specie cannot be derived here"
    groups = analysis._prepare_cluster_groups(df)
    assert len(groups) == 1, f"the cluster must survive grouping, got {len(groups)} groups"
    assert len(groups[0][1]) == 2


# ---------------------------------------------------------------------------
# -specie: stated by the user, cross-checked against the data
# ---------------------------------------------------------------------------

def _specie_frame(gene_ids):
    return pd.DataFrame({'gene_ensembl_id': gene_ids,
                         'start_position': list(range(len(gene_ids))),
                         'end_position': [x + 1 for x in range(len(gene_ids))],
                         'cluster_name': ['c'] * len(gene_ids)})


def test_stated_specie_is_authoritative():
    import utils
    out = utils.normalize_junctions_frame(_specie_frame(['ENSG1', 'ENSG2']), specie='human')
    assert list(out['specie']) == ['human', 'human']


def test_contradicting_gene_ids_abort():
    """A wrong -specie must stop the run, not silently mislabel every row."""
    import utils
    with pytest.raises(ValueError, match="does not match"):
        utils.normalize_junctions_frame(_specie_frame(['ENSG1', 'ENSG2']), specie='mouse')


def test_undeterminable_gene_ids_do_not_abort():
    """A GeneID names no species, so it cannot contradict the stated one. Treating
    silence as a mismatch would reject exactly the GeneID-only genes DOMAS supports."""
    import utils
    out = utils.normalize_junctions_frame(_specie_frame(['79400', '434223']), specie='mouse')
    assert list(out['specie']) == ['mouse', 'mouse']


def test_unknown_specie_rejected():
    import utils
    with pytest.raises(ValueError, match="Unknown specie"):
        utils.normalize_junctions_frame(_specie_frame(['ENSG1']), specie='dog')


def test_no_stated_specie_falls_back_to_derivation():
    """Library callers that predate the flag keep working."""
    import utils
    out = utils.normalize_junctions_frame(_specie_frame(['ENSG1', 'ENSMUSG2']))
    assert list(out['specie']) == ['human', 'mouse']


def test_database_specie_mismatch_aborts():
    """The database knows the species of every gene it holds, including GeneID-keyed
    ones the Ensembl prefix cannot read - so it catches a wrong -specie the prefix
    check is blind to."""
    from junction_analisys import _assert_specie_matches_database
    df = _specie_frame(['79400'])
    df['specie'] = 'mouse'
    with pytest.raises(ValueError, match="Species mismatch"):
        _assert_specie_matches_database(df, {'79400': 'H_sapiens'})


def test_database_specie_check_ignores_genes_it_does_not_hold():
    """A gene absent from the database says nothing about the species; it is left to
    the gene_not_in_db path rather than aborting the run."""
    from junction_analisys import _assert_specie_matches_database
    df = _specie_frame(['ENSG_UNKNOWN'])
    df['specie'] = 'human'
    _assert_specie_matches_database(df, {'79400': 'H_sapiens'})  # must not raise


# ---------------------------------------------------------------------------
# enrichment.py CLI - provisioning and -enrich argument handling
# ---------------------------------------------------------------------------

def test_enrichment_cli_requires_something_to_do():
    import enrichment
    with pytest.raises(SystemExit):
        enrichment.parse_args([])


def test_enrichment_cli_enrich_requires_dochap():
    import enrichment
    with pytest.raises(SystemExit):
        enrichment.parse_args(['-enrich', 'results.csv'])


def test_enrichment_cli_scores_requires_enrich():
    """-scores is a second pass over the enriched table, so it is meaningless alone."""
    import enrichment
    with pytest.raises(SystemExit):
        enrichment.parse_args(['-build', '-scores'])


def test_enrichment_cli_setup_defaults(tmp_path):
    """-enrichment_db and -pfam_hmm default to where -build leaves them."""
    import enrichment
    args = enrichment.parse_args(['-build', '-data_dir', str(tmp_path)])
    assert args.enrichment_db == os.path.join(str(tmp_path), 'enrichment.sqlite')
    assert args.pfam_hmm == os.path.join(str(tmp_path), 'pfam', 'Pfam-A.hmm')


def test_enrichment_cli_derives_output_name(tmp_path):
    """-out defaults beside the input rather than overwriting it."""
    import enrichment
    results = tmp_path / 'results.csv'
    results.write_text('event_type\n')
    db = tmp_path / 'enrichment.sqlite'; db.write_text('')
    dochap = tmp_path / 'dochap.sqlite'; dochap.write_text('')
    pfam = tmp_path / 'pfam'; pfam.mkdir(); (pfam / 'Pfam-A.hmm').write_text('')

    args = enrichment.parse_args(['-enrich', str(results), '-dochap', str(dochap),
                                  '-data_dir', str(tmp_path)])
    assert args.out == str(tmp_path / 'results.enriched.csv')
    assert args.out != str(results), "must not overwrite the input by default"


def test_enrichment_cli_reports_missing_inputs_before_running(tmp_path):
    """Paths are checked up front: an hmmsearch can run for minutes before it would
    otherwise discover a missing database."""
    import enrichment
    results = tmp_path / 'results.csv'; results.write_text('event_type\n')
    dochap = tmp_path / 'dochap.sqlite'; dochap.write_text('')
    with pytest.raises(SystemExit):
        enrichment.parse_args(['-enrich', str(results), '-dochap', str(dochap),
                               '-data_dir', str(tmp_path)])  # no enrichment.sqlite built


# ---------------------------------------------------------------------------
# filter_representative_domains: Site/PTM removal
#
# This branch is reachable with real data - 3,041 Site/PTM rows in the database,
# on 1,007 human transcripts - but no cluster in any fixture happens to sample one,
# so it had no coverage at all. The frames below are synthetic for that reason:
# they pin the rule rather than a particular gene.
# ---------------------------------------------------------------------------

def _typed_domain_row(domain_id, entry_type, start, end):
    return {'domain_id': domain_id, 'type': entry_type, 'AA_start': start, 'AA_end': end,
            **{c: None for c in DOMAIN_NAME_COLUMNS}}


def test_site_and_ptm_entries_are_removed():
    """Site and PTM entries are residue-level features, not structural units, so they
    are dropped outright - not merely demoted like a redundant Family entry."""
    from junction_analisys import filter_representative_domains
    df = pd.DataFrame([
        _typed_domain_row('IPR000001', 'Domain', 10, 200),
        _typed_domain_row('IPR000002', 'Conserved_site', 50, 62),
        _typed_domain_row('IPR000003', 'Binding_site', 70, 75),
        _typed_domain_row('IPR000004', 'PTM', 90, 92),
    ])
    kept = filter_representative_domains(df)
    assert list(kept['domain_id']) == ['IPR000001'], list(kept['domain_id'])


def test_site_entries_removed_even_where_no_domain_covers_them():
    """Removal is unconditional. The >50% coverage rule that spares an otherwise
    redundant Family or member entry - 'demote, don't delete' - does not apply here:
    a residue-level feature is never the evidence a domain comparison rests on, so a
    Conserved_site alone in its region is still dropped rather than kept as the only
    annotation."""
    from junction_analisys import filter_representative_domains
    df = pd.DataFrame([
        _typed_domain_row('IPR000001', 'Domain', 10, 100),
        _typed_domain_row('IPR000002', 'Conserved_site', 500, 512),   # nothing else near it
    ])
    kept = filter_representative_domains(df)
    assert list(kept['domain_id']) == ['IPR000001'], list(kept['domain_id'])


def _same_id_kept(span_a, span_b):
    """AA spans kept after filter_representative_domains() collapses two entries that
    share one accession."""
    from junction_analisys import filter_representative_domains
    df = pd.DataFrame([
        _typed_domain_row('IPR000001', 'Domain', *span_a),
        _typed_domain_row('IPR000001', 'Domain', *span_b),
    ])
    kept = filter_representative_domains(df)
    return sorted(zip(kept['AA_start'], kept['AA_end']))


def test_same_id_collapsed_when_overlap_is_a_majority_of_the_shorter():
    """Two entries of one accession overlapping by >= 50% of the shorter are one
    physical domain - the longer is kept."""
    # shorter is 100-160 (61 aa); overlap 100-160 is all of it
    assert _same_id_kept((10, 200), (100, 160)) == [(10, 200)]
    # shorter is 150-250 (101 aa); overlap 150-200 is 51 aa = 50.5%
    assert _same_id_kept((10, 200), (150, 250)) == [(10, 200)]


def test_same_id_kept_apart_when_overlap_is_below_half_the_shorter():
    """Below the threshold both survive: two tandem instances of a repeat that share
    a few boundary residues are two domains, not one. This is what the previous
    any-overlap rule collapsed."""
    # shorter is 190-300 (111 aa); overlap 190-200 is 11 aa = 9.9%
    assert _same_id_kept((10, 200), (190, 300)) == [(10, 200), (190, 300)]
    # a single shared residue
    assert _same_id_kept((10, 200), (200, 400)) == [(10, 200), (200, 400)]


def test_same_id_disjoint_positions_are_two_domains():
    """Unchanged by the threshold: no overlap at all stays two domains."""
    assert _same_id_kept((10, 100), (300, 400)) == [(10, 100), (300, 400)]


def test_non_site_tiers_survive_alongside_site_removal():
    """Removing the sites must not disturb the tier ladder: a Family entry that no
    Domain covers is still kept, while the sites go."""
    from junction_analisys import filter_representative_domains
    df = pd.DataFrame([
        _typed_domain_row('IPR000001', 'Domain', 10, 100),
        _typed_domain_row('IPR000002', 'Family', 400, 500),           # no Domain covers this
        _typed_domain_row('IPR000003', 'PTM', 45, 47),
    ])
    kept = filter_representative_domains(df)
    assert set(kept['domain_id']) == {'IPR000001', 'IPR000002'}, list(kept['domain_id'])


# ---------------------------------------------------------------------------
# Empty transcript set: a run where no event resolves a gene
# ---------------------------------------------------------------------------

def test_get_exons_for_transcripts_with_no_ids_returns_empty_frame(db_path):
    """A run whose events all fail to resolve a gene has no transcript ids at all,
    so the chunk loop never runs and there is nothing to concatenate. pd.concat([])
    raises "No objects to concatenate", which killed the run with an opaque pandas
    error instead of letting it finish and report each event's reason.

    It happens for real: every LeafCutter cluster unannotated (82 in the fixture), or
    a gene set absent from the database.
    """
    import sqlite3
    import utils
    if not os.path.exists(db_path):
        pytest.skip(f"DoChaP database not found at {db_path}")
    con = sqlite3.connect(db_path)
    try:
        df = utils.get_exons_for_transcripts(con, set())
    finally:
        con.close()

    assert len(df) == 0
    # The table's own shape, not a bare DataFrame(): downstream code reads these
    # columns, and an empty frame without them fails on attribute access instead.
    for column in ('transcript_ensembl_id', 'transcript_refseq_id',
                   'genomic_start_tx', 'genomic_end_tx', 'order_in_transcript'):
        assert column in df.columns, f"{column} missing from the empty result"


def test_get_exons_for_transcripts_empty_matches_populated_shape(db_path):
    """The empty result must carry exactly the columns a populated one does."""
    import sqlite3
    import utils
    if not os.path.exists(db_path):
        pytest.skip(f"DoChaP database not found at {db_path}")
    con = sqlite3.connect(db_path)
    try:
        empty = utils.get_exons_for_transcripts(con, set())
        populated = utils.get_exons_for_transcripts(con, {'ENST00000388866.8'})
    finally:
        con.close()

    assert len(populated) > 0, "fixture transcript should have exons"
    assert list(empty.columns) == list(populated.columns)


# ---------------------------------------------------------------------------
# Canonical precedence when a gene has more than one canonical transcript.
# ---------------------------------------------------------------------------

def _canonical_pick(canonical_ids, canonical_rank):
    """Which transcript analyze() adopts as the gene's canonical one."""
    df_gene_transcripts = pd.DataFrame({
        'transcript_ensembl_id': ['ENST_AAA', 'ENST_BBB', 'ENST_CCC'],
        'transcript_refseq_id': [None, None, None],
        'canonical': [canonical_rank.get(t, 0) for t in ('ENST_AAA', 'ENST_BBB', 'ENST_CCC')],
        'cds_start': [0, 0, 0],
        'cds_end': [1000, 1000, 1000],
    })
    exons = _two_exon_df((1, 0, 100), (2, 200, 300))
    exons['abs_start_CDS'] = [1, 51]
    exons['abs_end_CDS'] = [50, 100]
    empty_domains = pd.DataFrame(columns=['AA_start', 'AA_end'] + DOMAIN_NAME_COLUMNS)

    cluster_result = ClusterAnalysisResult('cluster_1', 'ENSG_TEST', 'TESTGENE')
    cluster_result.junctions = [(100, 200)]
    cluster_result.analyze(
        df_gene_transcripts, canonical_transcript_ids=set(canonical_ids),
        exon_lookup=lambda tid: exons, domain_lookup=lambda tid: empty_domains,
        canonical_rank=canonical_rank,
    )
    return cluster_result.canonical_transcript_id


def test_canonical_precedence_prefers_both_then_ensembl_then_refseq():
    """DoChaP flags 1=RefSeq, 2=Ensembl, 3=both. Reading RefSeq Select gave ~4,800
    mouse and ~7,100 rat genes two canonical transcripts, so the pick must be
    ordered, not alphabetical: both-sources first, then Ensembl, then RefSeq.
    The lowest id would pick ENST_AAA in every case below."""
    # RefSeq-flagged AAA vs Ensembl-flagged BBB -> Ensembl wins
    assert _canonical_pick(['ENST_AAA', 'ENST_BBB'],
                           {'ENST_AAA': 1, 'ENST_BBB': 2}) == 'ENST_BBB'
    # Ensembl-flagged AAA vs both-flagged BBB -> both wins
    assert _canonical_pick(['ENST_AAA', 'ENST_BBB'],
                           {'ENST_AAA': 2, 'ENST_BBB': 3}) == 'ENST_BBB'
    # all three ranks present -> the 3 wins
    assert _canonical_pick(['ENST_AAA', 'ENST_BBB', 'ENST_CCC'],
                           {'ENST_AAA': 1, 'ENST_BBB': 2, 'ENST_CCC': 3}) == 'ENST_CCC'


def test_canonical_precedence_ties_and_missing_rank_keep_lowest_id():
    """Equal flags, and the no-rank path used by callers that pass only a set,
    both keep the previous behaviour: the lowest id, never hash order."""
    assert _canonical_pick(['ENST_AAA', 'ENST_BBB'],
                           {'ENST_AAA': 2, 'ENST_BBB': 2}) == 'ENST_AAA'
    assert _canonical_pick(['ENST_BBB', 'ENST_AAA'],
                           {'ENST_AAA': 3, 'ENST_BBB': 3}) == 'ENST_AAA'
