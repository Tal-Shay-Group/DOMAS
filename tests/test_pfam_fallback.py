"""A domain row whose ext_id names no Pfam id takes one from its DomainType.

The same region is annotated by several sources, and which ids SpeciesDB
collapsed into one row varies between isoforms of the same gene. AKAP13 carries
PE/DAG-bd on both its canonical and its alternative protein, at 47 aa on each -
but the alternative's row reads 'pfam00130; IPR002219; smart00109' and the
canonical's reads 'smart00109; IPR002219'. Both are type_id 117, whose
DomainType names pfam00130. Reading ext_id alone, the domain existed on one side
only, and DOMAS reported domain_gain for a domain that never changed.

The fallback has to stay narrow: this source's value is that a Pfam-only set is
non-redundant, and resolving every cd/smart row through its DomainType would put
that redundancy straight back.
"""
import os
import sqlite3
import sys

import pandas as pd
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.normpath(os.path.join(TESTS_DIR, '..', 'code'))
sys.path.insert(0, CODE_DIR)

import utils  # noqa: E402


@pytest.fixture
def con():
    """The cases that matter, as the two tables actually spell them."""
    connection = sqlite3.connect(':memory:')
    connection.executescript("""
        CREATE TABLE DomainEvent (protein_ensembl_id TEXT, protein_refseq_id TEXT,
                                  type_id INTEGER, AA_start INTEGER, AA_end INTEGER,
                                  ext_id TEXT);
        CREATE TABLE DomainType (type_id INTEGER, pfam TEXT);
        INSERT INTO DomainType VALUES
            (117,  'pfam00130'),                                -- one accession
            (151,  'pfam00621'),                                -- one accession
            (11,   'pfam00169; pfam16652; pfam16457'),          -- several: ambiguous
            (9999, NULL);                                       -- none to fall back on
        INSERT INTO DomainEvent VALUES
            -- AKAP13: the alternative states Pfam, the canonical does not.
            ('ALT', NULL, 117,   54,  100, 'pfam00130; IPR002219; smart00109'),
            ('CAN', NULL, 117, 1792, 1838, 'smart00109; IPR002219'),
            -- one region, three sources, under a type that already has a Pfam row:
            -- only the stated one may survive, or one domain becomes three.
            ('CAN', NULL, 151, 1998, 2188, 'pfam00621; IPR000219'),
            ('CAN', NULL, 151, 1995, 2189, 'cd00160; IPR000219'),
            ('CAN', NULL, 151, 1998, 2190, 'smart00325; IPR000219'),
            -- a type naming several Pfam accessions: nothing to choose between them
            ('CAN', NULL,  11, 2232, 2335, 'smart00233; IPR001849'),
            -- a type naming none
            ('CAN', NULL, 9999, 10,   20,  'cd20878'),
            -- two non-Pfam rows describing ONE region: both are added, and
            -- _collapse_duplicate_spans() merges them by the 50% rule downstream.
            ('WIDE', NULL, 117,  10,   20, 'cd11111'),
            ('WIDE', NULL, 117,   5,   95, 'smart22222'),
            -- a tandem REPEAT: two occurrences of one family, far apart. Both
            -- must survive - collapsing these is what lost 32,631 rows.
            ('REP',  NULL, 117, 100,  150, 'smart00109'),
            ('REP',  NULL, 117, 500,  550, 'smart00109'),
            -- a repeat sharing ONE residue with a stated instance. The majority
            -- rule keeps it; "any overlap" would suppress a whole instance.
            ('EDGE', NULL, 151, 100,  200, 'pfam00621; IPR000219'),
            ('EDGE', NULL, 151, 200,  300, 'smart00325; IPR000219');
    """)
    connection.commit()
    return connection


def events(con):
    frame = utils._pfam_domain_events(con)
    return {(row.protein_ensembl_id, row.type_id, row.AA_start, row.AA_end): row.domain_id
            for row in frame.itertuples()}


def test_a_row_stating_its_pfam_id_keeps_it(con):
    assert events(con)[('ALT', 117, 54, 100)] == 'pfam00130'


def test_the_akap13_case_recovers_the_accession_from_the_domain_type(con):
    # Without this the canonical had no PE/DAG-bd row at all, and a 47 aa domain
    # present on both sides was reported as domain_gain.
    assert events(con)[('CAN', 117, 1792, 1838)] == 'pfam00130'


def test_a_region_a_stated_row_covers_gains_no_others(con):
    """cd00160 and smart00325 describe the same region as pfam00621 and are
    suppressed, so one RhoGEF domain does not become three. Suppressed HERE
    rather than left to the collapse, which keeps the LONGER span: that would
    report CDD's 187 residues in place of Pfam's own 184."""
    found = [key for key in events(con) if key[0] == 'CAN' and key[1] == 151]
    assert found == [('CAN', 151, 1998, 2188)], found


def test_a_type_naming_several_pfam_accessions_is_left_alone(con):
    # Choosing one would assert an identity the row never stated.
    assert not [key for key in events(con) if key[1] == 11]


def test_a_type_naming_no_pfam_accession_is_left_alone(con):
    assert not [key for key in events(con) if key[1] == 9999]


def test_rows_describing_one_region_are_left_for_the_overlap_rule(con):
    """Both are added here. Deciding they are one domain is
    _collapse_duplicate_spans()'s job, by the same _SAME_ID_OVERLAP majority the
    rest of DOMAS uses - this function must not apply a second, different rule."""
    found = sorted(key for key in events(con) if key[0] == 'WIDE')
    assert found == [('WIDE', 117, 5, 95), ('WIDE', 117, 10, 20)], found


def test_a_tandem_repeat_survives_as_two_instances(con):
    """The whole point of the positional test. Keyed on type_id instead, a
    protein kept one occurrence of a repeated family and lost the rest: PAPPA2
    reported 1 NL repeat where its canonical has 2, and 2 CCP where it has 4."""
    found = sorted(key for key in events(con) if key[0] == 'REP')
    assert found == [('REP', 117, 100, 150), ('REP', 117, 500, 550)], found


def test_sharing_one_residue_with_a_stated_row_is_not_enough_to_suppress(con):
    """_SAME_ID_OVERLAP, not any overlap: two tandem copies can share a boundary,
    and suppressing an instance over one residue is what the majority rule exists
    to prevent. 100-200 is stated; 200-300 overlaps it by a single residue."""
    found = sorted(key for key in events(con) if key[0] == 'EDGE')
    assert found == [('EDGE', 151, 100, 200), ('EDGE', 151, 200, 300)], found


def test_the_recovered_rows_carry_the_same_columns_as_the_stated_ones(con):
    frame = utils._pfam_domain_events(con)
    assert 'type_pfam' not in frame.columns, 'the working column must not leak out'
    assert set(frame.columns) >= {'protein_ensembl_id', 'protein_refseq_id', 'type_id',
                                  'AA_start', 'AA_end', 'ext_id', 'domain_id'}
    assert frame['domain_id'].notna().all()


# --- _overlaps_a_stated_instance, directly -----------------------------------
# The rule above is decided by one function, and it answers every fallback row
# against every stated row sharing its protein and accession as a single joined
# pass. The cases below are the ones that pass through its branches rather than
# its arithmetic, and which a fixture-driven test reaches only by accident.

def _frame(rows):
    """rows: (ensembl, refseq, domain_id, AA_start, AA_end)"""
    return pd.DataFrame(rows, columns=['protein_ensembl_id', 'protein_refseq_id',
                                       'domain_id', 'AA_start', 'AA_end'])


def test_nothing_is_suppressed_when_the_two_sides_share_no_protein_or_accession():
    """The join is empty - a different branch from "joined, but nothing overlaps",
    and the one a run hits whenever the stated set covers none of the fallback's
    proteins."""
    fallback = _frame([('P1', None, 'pfam00001', 10, 100),
                       ('P2', None, 'pfam00002', 10, 100)])
    # same coordinates, but another protein and another accession
    direct = _frame([('P9', None, 'pfam00001', 10, 100),
                     ('P1', None, 'pfam00009', 10, 100)])
    mask = utils._overlaps_a_stated_instance(fallback, direct)
    assert not mask.any()
    assert mask.index.equals(fallback.index)


def test_a_null_coordinate_suppresses_nothing():
    """AA_start is only cast to int downstream, so a null reaches this function.
    It must not suppress, and must not raise."""
    fallback = _frame([('P1', None, 'pfam00001', None, 100),
                       ('P1', None, 'pfam00001', 10, 100)])
    direct = _frame([('P1', None, 'pfam00001', 10, 100)])
    mask = utils._overlaps_a_stated_instance(fallback, direct)
    assert list(mask) == [False, True], 'the null row stays, the real overlap goes'


def test_one_overlapping_instance_is_enough_among_many():
    """A protein carrying a repeat has several stated rows for one accession. The
    fallback row is suppressed if ANY of them covers it, not all."""
    fallback = _frame([('P1', None, 'pfam00001', 300, 400)])
    direct = _frame([('P1', None, 'pfam00001', 10, 100),
                     ('P1', None, 'pfam00001', 150, 200),
                     ('P1', None, 'pfam00001', 310, 420)])   # this one covers it
    assert bool(utils._overlaps_a_stated_instance(fallback, direct).iloc[0])


def test_the_refseq_key_is_honoured_when_there_is_no_ensembl_id():
    """224,530 Pfam rows are reachable only through the RefSeq id. Keying on the
    Ensembl column alone would silently stop suppressing for all of them."""
    fallback = _frame([(None, 'XP_1', 'pfam00001', 10, 100)])
    direct = _frame([(None, 'XP_1', 'pfam00001', 10, 100)])
    assert bool(utils._overlaps_a_stated_instance(fallback, direct).iloc[0])
    # ...and a different RefSeq protein must not be confused with it
    other = _frame([(None, 'XP_2', 'pfam00001', 10, 100)])
    assert not bool(utils._overlaps_a_stated_instance(fallback, other).iloc[0])


def test_an_empty_side_is_not_an_error():
    empty = _frame([])
    populated = _frame([('P1', None, 'pfam00001', 10, 100)])
    assert not utils._overlaps_a_stated_instance(populated, empty).any()
    assert utils._overlaps_a_stated_instance(empty, populated).empty


# --- restricting the read to a run's own proteins -----------------------------
# DomainEvent holds 2.4M rows and a run uses a few thousand, so both queries take
# the run's protein ids where it has them. The answer must not depend on that:
# these check the restricted read against the unrestricted one, which is the
# behaviour every other test in this file describes.

def _proteins(*ids, refseq=()):
    return pd.DataFrame({'protein_ensembl_id': list(ids) + [None] * len(refseq),
                         'protein_refseq_id': [None] * len(ids) + list(refseq)})


def _rows(frame):
    return sorted(frame[['protein_ensembl_id', 'protein_refseq_id', 'type_id',
                         'AA_start', 'AA_end', 'domain_id']]
                  .itertuples(index=False, name=None))


def test_restricting_to_a_protein_gives_that_protein_s_rows_unchanged(con):
    """The whole point: same answer, less work. Every row the unrestricted read
    returns for CAN must come back, and nothing else."""
    everything = utils._pfam_domain_events(con)
    restricted = utils._pfam_domain_events(con, df_protein=_proteins('CAN'))
    expected = everything[everything['protein_ensembl_id'] == 'CAN']
    assert _rows(restricted) == _rows(expected)
    assert list(restricted.columns) == list(everything.columns), 'no working column may leak'
    assert '_rowid' not in restricted.columns


def test_the_akap13_recovery_survives_a_restricted_read(con):
    """The fallback's reason for existing has to work on the restricted path too -
    it is a JOIN against DomainType, and the restriction is added to it as well."""
    restricted = utils._pfam_domain_events(con, df_protein=_proteins('CAN', 'ALT'))
    recovered = [r for r in _rows(restricted) if r[0] == 'CAN' and r[5] == 'pfam00130']
    assert recovered == [('CAN', None, 117, 1792, 1838, 'pfam00130')], recovered


def test_a_repeat_still_survives_as_two_instances_when_restricted(con):
    assert len([r for r in _rows(utils._pfam_domain_events(con, df_protein=_proteins('REP')))
                if r[5] == 'pfam00130']) == 2


def test_no_protein_frame_reads_the_whole_table(con):
    """The default, and what the DB-inspection callers and these tests rely on."""
    assert _rows(utils._pfam_domain_events(con, df_protein=None)) == \
           _rows(utils._pfam_domain_events(con))
    assert _rows(utils._pfam_domain_events(con, df_protein=pd.DataFrame())) == \
           _rows(utils._pfam_domain_events(con))


def test_a_protein_named_in_neither_id_space_returns_nothing(con):
    assert utils._pfam_domain_events(con, df_protein=_proteins('NOT_A_PROTEIN')).empty


def test_ids_spanning_several_batches_are_not_returned_twice(monkeypatch):
    """A row comes back from every batch holding either of its ids, and a protein's
    two ids need not land in the same batch. With one id per batch, a row keyed on
    both would be read twice and read downstream as two instances of one domain."""
    connection = sqlite3.connect(':memory:')
    connection.executescript("""
        CREATE TABLE DomainEvent (protein_ensembl_id TEXT, protein_refseq_id TEXT,
                                  type_id INTEGER, AA_start INTEGER, AA_end INTEGER,
                                  ext_id TEXT);
        CREATE TABLE DomainType (type_id INTEGER, pfam TEXT);
        INSERT INTO DomainType VALUES (117, 'pfam00130');
        INSERT INTO DomainEvent VALUES ('ENSP1', 'NP_1', 117, 10, 90, 'pfam00130');
    """)
    connection.commit()
    monkeypatch.setattr(utils, '_PROTEIN_ID_BATCH', 1)
    # Each id list is sorted and sliced independently, so these put the row's
    # ensembl id in batch 0 (ENSP1 before ENSP_OTHER) and its refseq id in batch 1
    # (NP_0 before NP_1) - the arrangement that reads the row twice.
    frame = utils._pfam_domain_events(
        connection,
        df_protein=pd.DataFrame({'protein_ensembl_id': ['ENSP1', 'ENSP_OTHER'],
                                 'protein_refseq_id': ['NP_0', 'NP_1']}))
    assert len(frame) == 1, f'the row was read once per batch holding an id: {frame}'


def test_too_many_proteins_falls_back_to_reading_the_table(con, monkeypatch):
    """Past the limit an indexed lookup per protein costs more than the scan, so the
    restriction is dropped - and dropping it must not change the answer."""
    monkeypatch.setattr(utils, '_PROTEIN_FILTER_LIMIT', 1)
    # two proteins against a limit of one, so the restriction is dropped and the
    # read returns every protein in the table, not just these two
    assert _rows(utils._pfam_domain_events(con, df_protein=_proteins('CAN', 'ALT'))) == \
           _rows(utils._pfam_domain_events(con)), 'the scan must return everything'
    # ...and under the limit it is still restricted, so the test above is testing
    # the limit rather than a path that never restricts
    monkeypatch.setattr(utils, '_PROTEIN_FILTER_LIMIT', 250_000)
    assert _rows(utils._pfam_domain_events(con, df_protein=_proteins('CAN', 'ALT'))) != \
           _rows(utils._pfam_domain_events(con))
