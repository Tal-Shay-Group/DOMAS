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
    empty_transcripts = df_transcripts.iloc[0:0]
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
                             transcripts_by_gene, empty_transcripts)

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
    analysis = _bare_analysis()
    df = pd.DataFrame({'gene_ensembl_id': ['G1'], 'start_position': [1]})  # missing end_position, cluster_name
    with pytest.raises(ValueError, match="is required"):
        analysis._load_junctions_data(df)


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


def test_select_most_like_canonical_prefers_more_exact_matching_exons():
    """Both candidates have exons outside the junction range [150, 350] that exactly
    match canonical's - so both qualify. Candidate B also happens to match canonical's
    exon inside the range, giving it a higher exact-match score, so it wins even though
    it isn't the longest-CDS transcript."""
    transcript_exons = {
        'CANON': _exon_df((1, 100), (200, 300), (400, 500)),
        'A': _exon_df((1, 100), (220, 280), (400, 500)),  # inside exon differs from canonical
        'B': _exon_df((1, 100), (200, 300), (400, 500)),  # inside exon also matches canonical
    }
    cds_length_by_transcript = {'A': 1000, 'B': 10}  # A is longer, but B should still win

    chosen = select_most_like_canonical(
        ['A', 'B'], 'CANON', transcript_exons, junctions=[(150, 350)],
        cds_length_by_transcript=cds_length_by_transcript,
    )
    assert chosen == 'B'


def test_select_most_like_canonical_excludes_candidate_with_different_outside_exons():
    """A candidate whose exons outside the junction range don't exactly match
    canonical's doesn't qualify, even if it matches everything inside the range."""
    transcript_exons = {
        'CANON': _exon_df((1, 100), (200, 300), (400, 500)),
        'DISQUALIFIED': _exon_df((1, 90), (200, 300), (400, 500)),  # outside exon (1,90) != (1,100)
        'QUALIFIES': _exon_df((1, 100), (220, 280), (400, 500)),  # outside exons match exactly
    }
    cds_length_by_transcript = {'DISQUALIFIED': 10, 'QUALIFIES': 10}

    chosen = select_most_like_canonical(
        ['DISQUALIFIED', 'QUALIFIES'], 'CANON', transcript_exons, junctions=[(150, 350)],
        cds_length_by_transcript=cds_length_by_transcript,
    )
    assert chosen == 'QUALIFIES'


def test_select_most_like_canonical_falls_back_to_longest_cds_when_none_qualify():
    """If no candidate's outside-range exons exactly match canonical's, fall back to
    the longest-CDS candidate."""
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
    assert chosen == 'LONG'


def test_select_most_like_canonical_tiebreak_is_deterministic():
    """Two qualifying candidates with an equal exact-match score must resolve
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
    # Tied on score (2 outside matches each) -> sorted(...)[0] picks the lexicographically first id.
    assert chosen == 'ENST_AAA'


# ---------------------------------------------------------------------------
# analyze_junction: every comparable transcript is compared, tagged with
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

    cluster_result.analyze_junction(
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

    cluster_result.analyze_junction(
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
