"""Tests for the per-row projection_status contract on
utils._reproject_inherited_domains().

The point of the column is that a caller can tell a repositioned row from one
the projection could not validate, and both from a row whose provenance was
never established at all - all three come back carrying coordinates.
"""
import os
import sqlite3
import sys

import pandas as pd
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.normpath(os.path.join(TESTS_DIR, '..', 'code'))
sys.path.insert(0, CODE_DIR)

from utils import (  # noqa: E402
    _reproject_inherited_domains,
    PROJECTION_OWN,
    PROJECTION_PROJECTED,
    PROJECTION_TRUNCATED,
    PROJECTION_UNVALIDATED,
    PROJECTION_UNKNOWN,
    PROJECTION_COMPARABLE,
)
from junction_analisys import domain_coordinates_comparable  # noqa: E402


def _frame(flags, accessions=None):
    n = len(flags)
    return pd.DataFrame({
        'interpro_domains_are_own': flags,
        'protein_interpro_id': accessions or [f'ACC{i}' for i in range(n)],
        'transcript_ensembl_id_version': [f'ENST{i}.1' for i in range(n)],
        'AA_start': [10] * n,
        'AA_end': [50] * n,
        'length': [200] * n,
    })


def _empty_db():
    """Enough schema for the borrowed path to run and find no reference."""
    con = sqlite3.connect(':memory:')
    con.execute("CREATE TABLE Proteins (protein_interpro_id TEXT, "
                "transcript_ensembl_id TEXT, transcript_refseq_id TEXT, "
                "interpro_domains_are_own INTEGER)")
    con.execute("CREATE TABLE Genes (gene_ensembl_id TEXT, gene_GeneID_id TEXT, "
                "strand TEXT)")
    con.execute("CREATE TABLE Transcripts (transcript_ensembl_id TEXT, "
                "transcript_refseq_id TEXT, cds_start INTEGER, cds_end INTEGER, "
                "gene_ensembl_id TEXT, gene_GeneID_id TEXT)")
    return con


def test_comparable_set_includes_truncated():
    """A truncated domain's coordinates are as established as a projected one's;
    the partial encoding is the finding, not a doubt about the measurement."""
    assert PROJECTION_COMPARABLE == {PROJECTION_OWN, PROJECTION_PROJECTED,
                                     PROJECTION_TRUNCATED}
    for status in (PROJECTION_UNVALIDATED, PROJECTION_UNKNOWN):
        assert status not in PROJECTION_COMPARABLE


def test_status_values_are_distinct():
    values = {PROJECTION_OWN, PROJECTION_PROJECTED, PROJECTION_TRUNCATED,
              PROJECTION_UNVALIDATED, PROJECTION_UNKNOWN}
    assert len(values) == 5


def test_missing_column_is_unknown_not_own():
    """A database predating the isoform step knows nothing about provenance."""
    df = _frame([1, 1]).drop(columns=['interpro_domains_are_own'])
    out = _reproject_inherited_domains(None, df)
    assert list(out['projection_status']) == [PROJECTION_UNKNOWN] * 2


def test_empty_frame_still_has_the_column():
    out = _reproject_inherited_domains(None, _frame([]).iloc[0:0])
    assert 'projection_status' in out.columns
    assert out.empty


def test_null_flag_is_unknown_and_one_is_own():
    """The crux: a NULL flag must not be absorbed into 'own'."""
    df = _frame([1, None, 1])
    out = _reproject_inherited_domains(None, df)
    assert list(out['projection_status']) == [
        PROJECTION_OWN, PROJECTION_UNKNOWN, PROJECTION_OWN]


def test_null_flag_survives_the_borrowed_path():
    """With a borrowed row present the function takes the long path; the NULL
    row must still come out 'unknown' rather than 'own'."""
    con = _empty_db()
    df = _frame([1, None, 0])
    out = _reproject_inherited_domains(con, df)
    assert list(out['projection_status']) == [
        PROJECTION_OWN, PROJECTION_UNKNOWN, PROJECTION_UNVALIDATED]


def test_borrowed_without_reference_is_unvalidated_and_keeps_coordinates():
    con = _empty_db()
    df = _frame([0])
    out = _reproject_inherited_domains(con, df)
    assert list(out['projection_status']) == [PROJECTION_UNVALIDATED]
    assert out['AA_start'].iloc[0] == 10 and out['AA_end'].iloc[0] == 50


@pytest.mark.parametrize('flags', [[1, 1], [None, None], [0], [1, None, 0]])
def test_column_is_always_present(flags):
    con = _empty_db()
    out = _reproject_inherited_domains(con, _frame(flags))
    assert 'projection_status' in out.columns
    assert out['projection_status'].notna().all()


# --- the per-transcript comparison gate ------------------------------------

def _lookup(frames):
    return lambda transcript_id: frames.get(transcript_id, pd.DataFrame())


@pytest.mark.parametrize('statuses,expected', [
    ([PROJECTION_OWN, PROJECTION_PROJECTED, PROJECTION_TRUNCATED], True),
    ([PROJECTION_OWN, PROJECTION_UNVALIDATED], False),
    ([PROJECTION_OWN, PROJECTION_UNKNOWN], False),
    ([PROJECTION_TRUNCATED], True),
])
def test_gate_reads_every_status(statuses, expected):
    frame = pd.DataFrame({'projection_status': statuses,
                          'AA_start': [1] * len(statuses),
                          'AA_end': [9] * len(statuses)})
    assert domain_coordinates_comparable(_lookup({'T1': frame}), 'T1') is expected


def test_gate_passes_a_transcript_with_no_domains():
    assert domain_coordinates_comparable(_lookup({}), 'T1') is True


def test_gate_passes_a_frame_without_the_column():
    """The DomainEvent / DomainType source carries no provenance column."""
    frame = pd.DataFrame({'AA_start': [1], 'AA_end': [9]})
    assert domain_coordinates_comparable(_lookup({'T1': frame}), 'T1') is True


def test_gate_tolerates_a_nan_status():
    """A run mixing RepresentativeDomains with the legacy source leaves NaN on
    the legacy rows; those are not failures."""
    frame = pd.DataFrame({'projection_status': [PROJECTION_OWN, None],
                          'AA_start': [1, 1], 'AA_end': [9, 9]})
    assert domain_coordinates_comparable(_lookup({'T1': frame}), 'T1') is True


def test_gate_sees_a_domain_outside_any_window():
    """The reason the gate reads the whole frame: an unvalidated domain sitting
    at inherited coordinates far from the junction joins no identity group, so a
    per-group check would never see it."""
    frame = pd.DataFrame({'projection_status': [PROJECTION_OWN, PROJECTION_UNVALIDATED],
                          'AA_start': [10, 4000], 'AA_end': [50, 4100]})
    assert domain_coordinates_comparable(_lookup({'T1': frame}), 'T1') is False
