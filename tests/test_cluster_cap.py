"""-max_clusters keeps the most significant clusters, not an arbitrary slice.

A capped run analyses only part of its input, so which part it picks decides
what the user sees. Selecting the alphabetically-first N was quietly wrong on a
genome-wide input: the 100 first cluster names of a LeafCutter run are a
contiguous stretch of one chromosome, and on a real file only two of those 100
carried a gene annotation at all - the web page reported one event where the
same input, pre-filtered by the user, gave nineteen.

Every reader that can supply one now carries the tool's own statistic on
utils.SIGNIFICANCE_COLUMN, ordered ascending (smaller is more significant), and
this pins down both halves of that contract: the ordering here, and each
reader's translation into it.
"""
import os
import sys

import pandas as pd
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.normpath(os.path.join(TESTS_DIR, '..', 'code'))
sys.path.insert(0, CODE_DIR)

import alternative_splicing  # noqa: E402
import utils  # noqa: E402

SIGNIFICANCE = utils.SIGNIFICANCE_COLUMN


def junctions(rows, with_column=True):
    """A minimal junctions frame: (cluster_name, significance) per row."""
    frame = {'cluster_name': [name for name, _ in rows]}
    if with_column:
        frame[SIGNIFICANCE] = [value for _, value in rows]
    return pd.DataFrame(frame)


def kept(frame, max_clusters):
    return list(alternative_splicing._limit_clusters(frame, max_clusters)['cluster_name'])


def test_keeps_the_most_significant_not_the_first_by_name():
    # Were the cap still taking names in order it would keep aaa and mmm.
    frame = junctions([('zzz', 0.001), ('aaa', 0.5), ('mmm', 0.01)])
    assert sorted(kept(frame, 2)) == ['mmm', 'zzz']


def test_a_cluster_with_no_statistic_ranks_after_every_cluster_with_one():
    # LeafCutter writes no significance row for a cluster it could not test, and
    # 'untested' must not outrank a weak-but-measured result.
    frame = junctions([('aaa', float('nan')), ('zzz', 0.9)])
    assert kept(frame, 1) == ['zzz']


def test_ties_break_on_cluster_name_so_the_selection_is_reproducible():
    frame = junctions([('ccc', 0.1), ('aaa', 0.1), ('bbb', 0.1)])
    assert sorted(kept(frame, 2)) == ['aaa', 'bbb']


@pytest.mark.parametrize('with_column', [False, True], ids=['absent', 'all_null'])
def test_falls_back_to_cluster_name_where_the_format_names_no_statistic(with_column):
    # SUPPA's .ioe and a plain junctions CSV carry none; the column is then
    # missing or wholly null and the previous behaviour stands.
    frame = junctions([('ccc', float('nan')), ('aaa', float('nan')), ('bbb', float('nan'))],
                      with_column=with_column)
    assert sorted(kept(frame, 2)) == ['aaa', 'bbb']


def test_a_cluster_is_kept_whole_and_its_split_genes_rank_together():
    # _leafcutter_attach_genes() renames a multi-gene cluster to cluster:SYMBOL,
    # so one cluster arrives here as several. They share a statistic and must be
    # ranked once, with every junction of a kept cluster retained.
    frame = pd.DataFrame({'cluster_name': ['c:A', 'c:A', 'c:B', 'd', 'd'],
                          SIGNIFICANCE: [0.2, 0.2, 0.2, 0.9, 0.9]})
    assert kept(frame, 2) == ['c:A', 'c:A', 'c:B']


def test_no_cap_leaves_the_frame_untouched():
    frame = junctions([('a', 0.1)])
    assert len(alternative_splicing._limit_clusters(frame, 0)) == 1
    assert alternative_splicing._limit_clusters(None, 5) is None


# --- each reader's translation into the ascending column ---------------------

def test_rmats_carries_its_fdr():
    frame = utils.rmats2junctions(os.path.join(TESTS_DIR, 'rmats'))
    assert frame[SIGNIFICANCE].notna().all()
    assert frame[SIGNIFICANCE].between(0, 1).all()


def test_majiq_probability_is_inverted_so_it_reads_ascending():
    """MAJIQ's P(|dPSI|>=0.20) is larger for a stronger event, the opposite of an
    adjusted p-value, so the reader stores 1 - p. A round trip back through that
    inversion must land on the LSV's strongest junction."""
    path = os.path.join(TESTS_DIR, 'majiq', 'NveB_Mono_voila.txt')
    frame = utils.voila2junctions(path)
    assert frame[SIGNIFICANCE].notna().all()

    source = pd.read_csv(path, sep='\t', dtype=str)
    column = next(c for c in source.columns if c.startswith('P(|dPSI|>='))
    expected = {row['LSV ID'].strip():
                1.0 - max(float(v) for v in str(row[column]).split(';'))
                for _, row in source.iterrows()}
    got = dict(zip(frame['cluster_name'], frame[SIGNIFICANCE]))
    for lsv, value in got.items():
        assert value == pytest.approx(expected[lsv]), lsv


# --- the cap bounds events, not cluster-gene pairs ---------------------------

def test_a_multi_gene_event_takes_one_slot_and_brings_all_its_genes():
    """_leafcutter_attach_genes() splits a cluster naming several genes into one
    entry per gene. Counting those as separate clusters spent the cap on them and
    silently dropped other events: on a real CD4/monocyte run, -max_clusters 100
    analysed 94 events, not 100."""
    frame = pd.DataFrame({
        'cluster_name':              ['c1:A', 'c1:B', 'c1:C', 'c2', 'c3'],
        utils.SOURCE_CLUSTER_COLUMN: ['c1',   'c1',   'c1',   'c2', 'c3'],
        SIGNIFICANCE:                [0.01,   0.01,   0.01,   0.02, 0.03],
    })
    # Two events asked for: the 3-gene one and the next best - all four rows.
    assert kept(frame, 2) == ['c1:A', 'c1:B', 'c1:C', 'c2']
    # One event asked for: its three genes, and nothing else.
    assert kept(frame, 1) == ['c1:A', 'c1:B', 'c1:C']


def test_the_cap_falls_back_to_cluster_name_where_no_event_is_recorded():
    # Formats that never split keep the column empty (or lack it entirely), and
    # the pair IS the event.
    frame = pd.DataFrame({'cluster_name': ['c1', 'c2', 'c3'],
                          utils.SOURCE_CLUSTER_COLUMN: [None, None, None],
                          SIGNIFICANCE: [0.03, 0.01, 0.02]})
    assert kept(frame, 2) == ['c2', 'c3']


def test_the_summary_counts_events_and_pairs_separately():
    """One event naming three genes is three cluster-gene pairs and one event,
    and its junctions are counted once distinct and three times as rows."""
    import junction_analisys as ja
    frame = pd.DataFrame({
        'cluster_name':              ['c1:A', 'c1:B', 'c1:C', 'c2'],
        utils.SOURCE_CLUSTER_COLUMN: ['c1',   'c1',   'c1',   'c2'],
        'gene_symbol':               ['A',    'B',    'C',    'D'],
        'gene_ensembl_id':           ['A',    'B',    'C',    'D'],
        'specie':                    ['human'] * 4,
        'chromosome':                ['1'] * 4,
        'start_position':            [100, 100, 100, 300],
        'end_position':              [200, 200, 200, 400],
        utils.FEATURE_TYPE_COLUMN:   [utils.FEATURE_JUNCTION] * 4,
    })
    groups = list(frame.groupby(['specie', 'cluster_name'], dropna=False))
    summary = ja.RunSummary()
    summary.seed(frame, groups)
    assert summary.input_clusters == 2, summary.input_clusters
    assert summary.input_pairs == 4, summary.input_pairs
    assert summary.input_junctions == 4
    assert summary.input_junctions_distinct == 2, summary.input_junctions_distinct
    # and the breakdown sees the 3-gene event as one event naming three genes
    assert summary.genes_per_cluster == {3: 1, 1: 1}, summary.genes_per_cluster
