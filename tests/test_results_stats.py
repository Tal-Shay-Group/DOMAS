"""
End-to-end coverage for the Hadas-analysis -> stats-report flow: Hadas Excel
input -> alternative_splicing.analyze_hadas_input() -> results.csv ->
results_stats._run_pdf_report() with a custom label (the same flow used for
one-off named datasets, e.g. "hadas2"). Nothing in test_flags.py exercises
results_stats.py's reporting layer - its tests are scoped to
JunctionsAnalysis.analyze_junctions() itself - so this is the first coverage
for that half of the pipeline.
"""
import os
import sqlite3
import sys

import matplotlib
matplotlib.use('Agg')  # headless rendering for PDF generation in tests

import pandas as pd
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.normpath(os.path.join(TESTS_DIR, '..', 'code'))
sys.path.insert(0, CODE_DIR)

from alternative_splicing import analyze_hadas_input  # noqa: E402
import results_stats  # noqa: E402

HADAS_XLSX = os.path.join(TESTS_DIR, 'short_H_vs_M_HN6.xlsx')


@pytest.fixture(scope='module')
def con(db_path):
    if not os.path.exists(db_path):
        pytest.skip(f"DoChaP database not found at {db_path}")
    connection = sqlite3.connect(db_path)
    yield connection
    connection.close()


def test_hadas_analysis_and_stats_report_flow(con, tmp_path):
    """
    Full flow: analyze_hadas_input(use_representative_domains=True,
    create_pdf=False, num_workers=2) against a real Hadas fixture, then feed
    the resulting results.csv straight into _run_pdf_report() under a custom
    label - exactly how a one-off named dataset (e.g. "hadas2") gets
    analyzed and reported on.

    num_workers=2 deliberately exercises analyze_junctions()'s
    ProcessPoolExecutor path with more than one worker, which the rest of
    this suite never does (test_flags.py always uses num_workers=1) - safe
    here since pytest is itself a properly-guarded multiprocessing entry
    point, unlike a bare script with module-level code and no
    `if __name__ == "__main__":` guard (which hangs/crashes under macOS's
    spawn-based multiprocessing - the actual root cause the one time this
    flow was run from an ad-hoc script instead of through pytest).
    """
    results_csv = str(tmp_path / 'results.csv')
    analyze_hadas_input(
        con, HADAS_XLSX, results_csv,
        use_representative_domains=True,
        create_pdf=False,
        num_workers=2,
    )

    assert os.path.exists(results_csv), "analyze_hadas_input should have written a results.csv"
    df = pd.read_csv(results_csv)
    assert len(df) > 0, "Expected at least one analyzed row"
    assert 'event_type' in df.columns

    stats_pdf = str(tmp_path / 'stats_report.pdf')
    label_dfs = results_stats._run_pdf_report(
        [("test_hadas", results_csv, None)],
        stats_pdf,
        "test_hadas Full Report",
    )

    assert "test_hadas" in label_dfs
    assert len(label_dfs["test_hadas"]) == len(df)
    assert os.path.exists(stats_pdf), "_run_pdf_report should have written a PDF"

    PyPDF2 = pytest.importorskip('PyPDF2')
    reader = PyPDF2.PdfReader(stats_pdf)
    assert len(reader.pages) > 1, "Expected more than just a title page in the stats report"
