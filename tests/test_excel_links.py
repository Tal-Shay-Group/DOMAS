"""The workbook's gene cells link to DoChaP, and no row is written twice.

Two things the CSV cannot carry and the Excel copy can: a link behind the gene
symbol, the way the web GUI's results table does it, and - since the same
oversight produced both - rows that actually differ from one another.
"""
import os
import sys

import pandas as pd
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.normpath(os.path.join(TESTS_DIR, '..', 'code'))
sys.path.insert(0, CODE_DIR)

import junction_analisys as ja  # noqa: E402


def test_the_url_matches_what_the_web_gui_builds():
    # domasController.geneHref: origin + '#!/results/' + db specie + '/' + gene,
    # then the row's two transcripts, comma-joined and encoded as one segment.
    assert ja.dochap_gene_url('MBD2', 'human', 'ENST00000256429.8', 'ENST00000578272.1') == (
        'https://dochap.bgu.ac.il/#!/results/H_sapiens/MBD2/'
        'ENST00000256429.8%2CENST00000578272.1')


def test_the_species_is_translated_to_the_database_spelling():
    assert '/M_musculus/' in ja.dochap_gene_url('Mbd2', 'mouse')
    # An unknown species falls back to 'all', which the server reads as
    # "do not filter by species" - not to a species that would match nothing.
    assert '/all/' in ja.dochap_gene_url('MBD2', 'martian')


def test_a_row_naming_no_transcript_still_links_to_the_gene():
    assert ja.dochap_gene_url('TNFAIP8L2-SCNM1', 'human') == (
        'https://dochap.bgu.ac.il/#!/results/H_sapiens/TNFAIP8L2-SCNM1')


@pytest.mark.parametrize('gene', [None, '', '   ', 'nan', 'NaN'], ids=
                         ['none', 'empty', 'blank', 'nan', 'NaN'])
def test_a_row_naming_no_gene_gets_no_link(gene):
    # no_gene_specified rows are a real LeafCutter outcome, not bad input, and a
    # link to the empty string would 404.
    assert ja.dochap_gene_url(gene, 'human', 'a', 'b') is None


def test_the_workbook_keeps_the_symbol_as_its_text(tmp_path):
    openpyxl = pytest.importorskip('openpyxl')
    csv = tmp_path / 'r.csv'
    pd.DataFrame({
        'event': ['c1', 'c2'],
        'gene_symbol': ['MBD2', ''],
        'species': ['human', 'human'],
        'canonical_transcript_id': ['ENST1', ''],
        'alternative_transcript_id': ['ENST2', ''],
    }).to_csv(csv, index=False)

    target = ja.write_excel_copy(str(csv))
    sheet = openpyxl.load_workbook(target).active
    column = [c.value for c in sheet[1]].index('gene_symbol') + 1

    linked = sheet.cell(row=2, column=column)
    assert linked.value == 'MBD2', 'the cell must read as the gene, not a formula'
    assert linked.hyperlink.target.endswith('/MBD2/ENST1%2CENST2')
    # and the row naming no gene is left alone
    assert sheet.cell(row=3, column=column).hyperlink is None


def test_a_linked_cell_is_drawn_as_a_link(tmp_path):
    """Underlined, in an explicit rgb, at the size of the row it sits in.

    openpyxl's builtin 'Hyperlink' named style - which this used to apply - gives
    a cell `<color theme="10"/>` with no underline at size 12. In Excel that is
    blue and a point too large; in Numbers and Google Sheets, which do not
    reliably resolve the theme's hyperlink slot, it is indistinguishable from
    plain text, and the links were reported as missing.
    """
    openpyxl = pytest.importorskip('openpyxl')
    csv = tmp_path / 'r.csv'
    pd.DataFrame({
        'event': ['c1'],
        'gene_symbol': ['MBD2'],
        'species': ['human'],
        'canonical_transcript_id': ['ENST1'],
        'alternative_transcript_id': ['ENST2'],
    }).to_csv(csv, index=False)

    sheet = openpyxl.load_workbook(ja.write_excel_copy(str(csv))).active
    column = [c.value for c in sheet[1]].index('gene_symbol') + 1
    cell = sheet.cell(row=2, column=column)

    assert cell.font.underline == 'single'
    assert cell.font.color.rgb == ja.EXCEL_HYPERLINK_RGB
    # 'rgb', not 'theme': a theme reference is what the reader has to resolve
    assert cell.font.color.type == 'rgb'
    # copy()'d from the cell, so the link matches the row rather than resizing it
    neighbour = sheet.cell(row=2, column=1)
    assert cell.font.size == neighbour.font.size
    assert cell.font.name == neighbour.font.name
    # the header is not a link and must keep its own styling
    assert sheet.cell(row=1, column=column).font.underline is None


def test_a_sheet_over_excels_link_ceiling_is_written_unlinked(tmp_path, monkeypatch):
    """Excel refuses to open a worksheet with too many links, so the workbook is
    written without them rather than written broken."""
    pytest.importorskip('openpyxl')
    monkeypatch.setattr(ja, 'EXCEL_MAX_HYPERLINKS', 2)
    frame = pd.DataFrame({'gene_symbol': ['A', 'B', 'C'], 'species': ['human'] * 3})
    assert ja._link_gene_cells(str(tmp_path / 'unused.xlsx'), frame) is None


def test_every_novel_junction_row_names_its_own_junction():
    """One row per unmapped junction, and the row says which. Before it did not,
    so a cluster with six unmapped junctions wrote six identical rows - which
    read as duplication and lost the only thing those rows had to say."""
    result = ja.ClusterAnalysisResult('c1', 'ENSG1', 'G1', specie='human', strand='+')
    result.junctions = [(100, 200), (300, 400)]
    result.matched_features = {}
    for junction in result.junctions:
        low, high = min(junction), max(junction)
        result.add_event('novel_junction', None, alternative_junctions=f'{low}-{high}')

    written = [event[3] for event in result.events]
    assert written == ['100-200', '300-400']
    assert len(set(written)) == len(written), 'the rows must differ'
