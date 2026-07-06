"""
Tests for the use_longest_cds and restrict_pdf_to_comparable flags on
JunctionsAnalysis.analyze_junctions(), using the real fixture files
ioe_example_junctions.csv (plain CSV format) and short_H_vs_M_HN6.xlsx
(hadas format) against the local DoChaP database.
"""
import glob
import json
import os
import shutil
import sqlite3
import sys

import matplotlib
matplotlib.use('Agg')  # headless rendering for PDF generation in tests

import pandas as pd
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.normpath(os.path.join(TESTS_DIR, '..', 'code'))
sys.path.insert(0, CODE_DIR)

from junction_analisys import JunctionsAnalysis, ClusterAnalysisResult  # noqa: E402
from pdf_text_utils import (  # noqa: E402
    PDF_TEXT_MANIFEST_FILENAME, build_pdf_text_manifest,
)

IOE_CSV = os.path.join(TESTS_DIR, 'ioe_example_junctions.csv')
HADAS_XLSX = os.path.join(TESTS_DIR, 'short_H_vs_M_HN6.xlsx')

# Mirrors JunctionsAnalysis._SKIPPED_TRANSCRIPT_EVENTS plus the cluster-level
# events that carry no real transcript id.
_SKIPPED_EVENTS = {
    'transcript_doesnt_have_junctions', 'no_unique_junctions', 'skipped_not_longest_cds',
    'no_canonical_transcript', 'only_one_transcript', 'no_canonical_junctions', 'junction_not_mapped',
}

FLAG_COMBINATIONS = [
    (False, False),
    (False, True),
    (True, False),
    (True, True),
]

INPUT_FILES = [
    ('ioe_csv', IOE_CSV, False),
    ('hadas_xlsx', HADAS_XLSX, True),
]

REFERENCE_OUTPUTS_DIR = os.path.join(TESTS_DIR, 'reference_outputs')
GENERATED_OUTPUTS_DIR = os.path.join(TESTS_DIR, 'generated_outputs')


@pytest.fixture(scope='module')
def con(db_path):
    if not os.path.exists(db_path):
        pytest.skip(f"DoChaP database not found at {db_path}")
    connection = sqlite3.connect(db_path)
    yield connection
    connection.close()


def _compared_transcript_ids(cluster_result):
    """Transcript ids that were actually compared to the canonical transcript (not skipped)."""
    df = cluster_result.get_results_df()
    if len(df) == 0:
        return set()
    mask = ~df['event'].isin(_SKIPPED_EVENTS)
    return set(df.loc[mask, 'transcript_id'].dropna())


def _assert_longest_cds_invariant(results, use_longest_cds):
    """When use_longest_cds is True, no cluster should compare more than one transcript."""
    if not use_longest_cds:
        return
    for cluster_result in results:
        compared = _compared_transcript_ids(cluster_result)
        assert len(compared) <= 1, (
            f"Cluster {cluster_result.cluster_name} compared {len(compared)} transcripts "
            f"({compared}) although use_longest_cds=True"
        )


def _run_analysis(con, tmp_path, junctions_csv, hadas_format, use_longest_cds, restrict_pdf_to_comparable):
    """Run analyze_junctions with cwd set to tmp_path (PDFs are written relative to cwd)."""
    analysis = JunctionsAnalysis(con)
    output_path = str(tmp_path / 'results.csv')
    cwd_before = os.getcwd()
    os.chdir(tmp_path)
    try:
        results = analysis.analyze_junctions(
            junctions_csv=junctions_csv,
            output_path=output_path,
            hadas_format=hadas_format,
            create_pdf=True,
            num_workers=1,
            use_longest_cds=use_longest_cds,
            restrict_pdf_to_comparable=restrict_pdf_to_comparable,
        )
    finally:
        os.chdir(cwd_before)
    return results, output_path


@pytest.mark.parametrize('use_longest_cds,restrict_pdf_to_comparable', FLAG_COMBINATIONS)
def test_ioe_csv_all_flag_combinations(con, tmp_path, use_longest_cds, restrict_pdf_to_comparable):
    """analyze_junctions runs end-to-end on ioe_example_junctions.csv for every combination of the two flags."""
    results, output_path = _run_analysis(
        con, tmp_path, IOE_CSV, hadas_format=False,
        use_longest_cds=use_longest_cds, restrict_pdf_to_comparable=restrict_pdf_to_comparable,
    )

    assert len(results) > 0, "Expected at least one cluster to be analyzed"
    assert os.path.exists(output_path), "Results CSV should have been written"

    pdf_files = glob.glob(str(tmp_path / '*_junction_comparison.pdf'))
    assert len(pdf_files) > 0, "Expected at least one PDF to be generated"

    _assert_longest_cds_invariant(results, use_longest_cds)


@pytest.mark.parametrize('use_longest_cds,restrict_pdf_to_comparable', FLAG_COMBINATIONS)
def test_hadas_xlsx_all_flag_combinations(con, tmp_path, use_longest_cds, restrict_pdf_to_comparable):
    """analyze_junctions runs end-to-end on short_H_vs_M_HN6.xlsx for every combination of the two flags."""
    results, output_path = _run_analysis(
        con, tmp_path, HADAS_XLSX, hadas_format=True,
        use_longest_cds=use_longest_cds, restrict_pdf_to_comparable=restrict_pdf_to_comparable,
    )

    assert len(results) > 0, "Expected at least one cluster to be analyzed"
    assert os.path.exists(output_path), "Results CSV should have been written"

    pdf_files = glob.glob(str(tmp_path / '*_junction_comparison.pdf'))
    assert len(pdf_files) > 0, "Expected at least one PDF to be generated"

    _assert_longest_cds_invariant(results, use_longest_cds)


# ---------------------------------------------------------------------------
# Automatic comparison against the committed golden reference outputs
# ---------------------------------------------------------------------------

def _case_name(label, use_longest_cds, restrict_pdf_to_comparable):
    return f"{label}__longestcds_{use_longest_cds}__restrict_{restrict_pdf_to_comparable}"


def _run_case_to_dir(con, case_dir, junctions_csv, hadas_format, use_longest_cds, restrict_pdf_to_comparable):
    """Run analyze_junctions, writing results.csv + PDFs into a freshly-created case_dir."""
    if os.path.exists(case_dir):
        shutil.rmtree(case_dir)
    os.makedirs(case_dir)

    output_path = os.path.join(case_dir, 'results.csv')
    cwd_before = os.getcwd()
    os.chdir(case_dir)
    try:
        analysis = JunctionsAnalysis(con)
        analysis.analyze_junctions(
            junctions_csv=junctions_csv,
            output_path=output_path,
            hadas_format=hadas_format,
            create_pdf=True,
            num_workers=1,
            use_longest_cds=use_longest_cds,
            restrict_pdf_to_comparable=restrict_pdf_to_comparable,
        )
    finally:
        os.chdir(cwd_before)
    return output_path


def _compare_csv_to_reference(generated_csv, reference_csv):
    """Assert generated_csv has the same rows as reference_csv, ignoring row order."""
    if not os.path.exists(reference_csv):
        pytest.skip(f"No reference output committed to compare against at {reference_csv}")

    df_generated = pd.read_csv(generated_csv).fillna('')
    df_reference = pd.read_csv(reference_csv).fillna('')
    sort_columns = list(df_reference.columns)
    df_generated = df_generated.sort_values(sort_columns).reset_index(drop=True)
    df_reference = df_reference.sort_values(sort_columns).reset_index(drop=True)
    pd.testing.assert_frame_equal(df_generated, df_reference, check_dtype=False)


def _compare_pdf_text_to_reference(output_dir, reference_dir):
    """Assert every PDF's extracted text in output_dir matches the golden manifest
    in reference_dir, comparing text content rather than raw PDF bytes (which
    differ run-to-run due to matplotlib's embedded CreationDate).
    """
    pytest.importorskip('PyPDF2')
    reference_manifest_path = os.path.join(reference_dir, PDF_TEXT_MANIFEST_FILENAME)
    if not os.path.exists(reference_manifest_path):
        pytest.skip(f"No PDF text reference committed to compare against at {reference_manifest_path}")

    with open(reference_manifest_path) as f:
        reference_manifest = json.load(f)
    generated_manifest = build_pdf_text_manifest(output_dir)

    missing = sorted(set(reference_manifest) - set(generated_manifest))
    extra = sorted(set(generated_manifest) - set(reference_manifest))
    mismatched = sorted(
        name for name in reference_manifest.keys() & generated_manifest.keys()
        if reference_manifest[name] != generated_manifest[name]
    )
    assert not missing and not extra and not mismatched, (
        f"PDF text mismatch vs {reference_manifest_path}: "
        f"missing={missing}, extra={extra}, mismatched_content={mismatched}"
    )


@pytest.mark.parametrize('label,junctions_csv,hadas_format', INPUT_FILES)
@pytest.mark.parametrize('use_longest_cds,restrict_pdf_to_comparable', FLAG_COMBINATIONS)
def test_compare_against_reference_outputs(con, keep_test_output, label, junctions_csv, hadas_format,
                                            use_longest_cds, restrict_pdf_to_comparable):
    """Run analyze_junctions for this flag combination, writing its CSV/PDF output to
    tests/generated_outputs/<case_name>/, then compare the resulting results.csv
    against the golden reference under tests/reference_outputs/<case_name>/.

    The generated output directory is deleted after a successful comparison by
    default. Pass --keep-test-output on the pytest command line to keep it
    (e.g. to inspect the PDFs by hand).
    """
    case_name = _case_name(label, use_longest_cds, restrict_pdf_to_comparable)
    output_dir = os.path.join(GENERATED_OUTPUTS_DIR, case_name)
    reference_dir = os.path.join(REFERENCE_OUTPUTS_DIR, case_name)
    reference_csv = os.path.join(reference_dir, 'results.csv')

    generated_csv = _run_case_to_dir(
        con, output_dir, junctions_csv, hadas_format, use_longest_cds, restrict_pdf_to_comparable,
    )

    pdf_files = glob.glob(os.path.join(output_dir, '*_junction_comparison.pdf'))
    assert len(pdf_files) > 0, f"Expected at least one PDF to be generated in {output_dir}"

    _compare_csv_to_reference(generated_csv, reference_csv)
    _compare_pdf_text_to_reference(output_dir, reference_dir)

    if not keep_test_output:
        shutil.rmtree(output_dir)


# ---------------------------------------------------------------------------
# Focused, DB-independent unit tests for the two new pieces of logic
# ---------------------------------------------------------------------------

def test_comparable_transcript_ids_excludes_skipped_events():
    """JunctionsAnalysis._comparable_transcript_ids keeps the canonical id plus only the
    transcripts that were actually compared, dropping every skip-event transcript."""
    analysis = JunctionsAnalysis.__new__(JunctionsAnalysis)  # no real db connection needed

    cluster_result = ClusterAnalysisResult('cluster_1', 'ENSG00001', 'GENE1')
    cluster_result.canonical_transcript_id = 'ENST_CANON'
    cluster_result.add_event('added_domain', transcript_id='ENST_COMPARED')
    cluster_result.add_event('no_domains_in_region', transcript_id='ENST_COMPARED_2')
    cluster_result.add_event('transcript_doesnt_have_junctions', transcript_id='ENST_NO_JUNCTIONS')
    cluster_result.add_event('no_unique_junctions', transcript_id='ENST_NOT_UNIQUE')
    cluster_result.add_event('skipped_not_longest_cds', transcript_id='ENST_SHORTER_CDS')

    comparable_ids = analysis._comparable_transcript_ids(cluster_result)

    assert comparable_ids == {'ENST_CANON', 'ENST_COMPARED', 'ENST_COMPARED_2'}


def test_comparable_transcript_ids_handles_no_events():
    """With no events recorded at all, only the canonical transcript id is comparable."""
    analysis = JunctionsAnalysis.__new__(JunctionsAnalysis)
    cluster_result = ClusterAnalysisResult('cluster_1', 'ENSG00001', 'GENE1')
    cluster_result.canonical_transcript_id = 'ENST_CANON'

    assert analysis._comparable_transcript_ids(cluster_result) == {'ENST_CANON'}


def test_transcript_matches_ids_checks_both_id_columns():
    """GeneVisualization._transcript_matches_ids matches on either ensembl or refseq id."""
    from generate_gene_pdf import GeneVisualization

    viz = GeneVisualization(conn=None, gene_name='GENE1')

    ensembl_only = {'info': pd.Series({'transcript_ensembl_id': 'ENST1', 'transcript_refseq_id': None})}
    refseq_only = {'info': pd.Series({'transcript_ensembl_id': None, 'transcript_refseq_id': 'NM_1'})}
    neither = {'info': pd.Series({'transcript_ensembl_id': 'ENST2', 'transcript_refseq_id': 'NM_2'})}

    wanted_ids = {'ENST1', 'NM_1'}
    assert viz._transcript_matches_ids(ensembl_only, wanted_ids) is True
    assert viz._transcript_matches_ids(refseq_only, wanted_ids) is True
    assert viz._transcript_matches_ids(neither, wanted_ids) is False


def test_create_pdf_transcript_ids_filters_transcripts(tmp_path, capsys):
    """create_pdf(transcript_ids=...) with no matching transcripts produces no PDF
    and a clear message instead of raising."""
    from generate_gene_pdf import GeneVisualization

    gene_data = pd.Series({
        'gene_ensembl_id': 'ENSG00001', 'gene_symbol': 'GENE1',
        'chromosome': '1', 'strand': '+', 'specie': 'H_sapiens',
    })
    preloaded = {'gene_data': gene_data, 'transcripts': []}
    viz = GeneVisualization(conn=None, gene_name='GENE1', preloaded=preloaded)
    viz.load_gene_data()

    output_file = str(tmp_path / 'no_match.pdf')
    viz.create_pdf(output_file, transcript_ids={'NOT_A_REAL_TRANSCRIPT'})

    captured = capsys.readouterr()
    assert 'None of the specified transcripts were found' in captured.out
    assert not os.path.exists(output_file), "No PDF should be written when there are no valid transcripts"


def test_create_pdf_transcript_ids_filters_real_gene(con, tmp_path):
    """create_pdf(transcript_ids=...) against a real, multi-transcript gene only draws
    the requested transcript(s) - verified by the resulting PDF's page count.
    """
    PyPDF2 = pytest.importorskip('PyPDF2')
    from generate_gene_pdf import GeneVisualization, prepare_gene_data_bulk

    gene_ensembl_id = 'ENSG00000174456'  # C12orf76 - has 12 ensembl transcripts in the DB
    one_transcript_id = 'ENST00000615315.2'

    preloaded_all = prepare_gene_data_bulk(con, [gene_ensembl_id])
    assert gene_ensembl_id in preloaded_all, f"Fixture gene {gene_ensembl_id} not found in the local DB"
    num_real_transcripts = len(preloaded_all[gene_ensembl_id]['transcripts'])
    assert num_real_transcripts > 4, (
        "Fixture gene is expected to have more transcripts than fit on one page "
        "(transcripts_per_page default is 4); the local DB content may have changed."
    )

    unrestricted_pdf = str(tmp_path / 'unrestricted.pdf')
    viz_all = GeneVisualization(con, 'C12orf76', preloaded=preloaded_all[gene_ensembl_id])
    viz_all.create_pdf(unrestricted_pdf)
    unrestricted_pages = len(PyPDF2.PdfReader(unrestricted_pdf).pages)

    restricted_pdf = str(tmp_path / 'restricted.pdf')
    preloaded_again = prepare_gene_data_bulk(con, [gene_ensembl_id])
    viz_one = GeneVisualization(con, 'C12orf76', preloaded=preloaded_again[gene_ensembl_id])
    viz_one.create_pdf(restricted_pdf, transcript_ids={one_transcript_id})
    restricted_pages = len(PyPDF2.PdfReader(restricted_pdf).pages)

    assert unrestricted_pages > 1, "Sanity check: unrestricted PDF should span more than one page"
    assert restricted_pages == 1, "Restricting to a single transcript should always fit on one page"
    assert restricted_pages < unrestricted_pages


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-v']))
