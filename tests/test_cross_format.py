"""Cross-format consistency: the same event, written in each tool's syntax, must
normalise to the same set of features.

Every reader funnels into one analysis, so if two inputs produce the same junctions
frame the results are identical by construction. What this pins down is the part
that is NOT shared - each reader's translation from its own file format into
(chromosome, start, end, feature_type). A disagreement here means the same biology
analysed two ways gives two answers, which is exactly the class of bug that made
every minus-strand rMATS A5SS/A3SS event unanalysable: SUPPA described those events
correctly while rMATS collapsed them, and nothing compared the two.

The events below are stated once as intron coordinates and then rendered into each
format, rather than lifted from a fixture, so the expected answer is written down
independently of any reader.
"""
import os
import sqlite3
import sys
import tempfile

import pandas as pd
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.normpath(os.path.join(TESTS_DIR, '..', 'code'))
sys.path.insert(0, CODE_DIR)

import alternative_splicing  # noqa: E402
import utils  # noqa: E402


# --- writers: one event, rendered into each format -------------------------

def _rmats_dir(tmp, event_type, columns, strand, chrom='21'):
    """A one-row [EventType].MATS.JC.txt carrying the event's exon coordinates."""
    row = {'ID': 1, 'GeneID': 'ENSG00000000001', 'geneSymbol': 'GENEX',
           'chr': f'chr{chrom}', 'strand': strand}
    row.update(columns)
    path = os.path.join(tmp, utils._RMATS_EVENT_FILES[event_type])
    pd.DataFrame([row]).to_csv(path, sep='\t', index=False)
    return utils.rmats2junctions(tmp)


def _suppa(tmp, event_id, chrom='21'):
    path = os.path.join(tmp, 'events.ioe')
    with open(path, 'w') as handle:
        handle.write("seqname\tgene_id\tevent_id\talternative_transcripts\ttotal_transcripts\n")
        handle.write(f"{chrom}\tENSG00000000001\t{event_id}\tT1\tT1,T2\n")
    return utils.ioe2junctions(path)


def _majiq(tmp, junction_coords, ir_coords='', strand='+', chrom='21'):
    path = os.path.join(tmp, 'voila.tsv')
    with open(path, 'w') as handle:
        handle.write("#Gene Name\tGene ID\tLSV ID\tchr\tstrand\tA5SS\tA3SS\tES\t"
                     "Junctions coords\tIR coords\n")
        handle.write(f"GENEX\tENSG00000000001\tENSG00000000001:s:1\tchr{chrom}\t{strand}\t"
                     f"False\tFalse\tTrue\t{junction_coords}\t{ir_coords}\n")
    return utils.voila2junctions(path)


def _leafcutter(tmp, introns, strand='+', chrom='21'):
    """The leafcutter_ds pair for one cluster of `introns` [(start, end), ...].

    Unlike the other three readers this one takes a DB connection, to resolve the
    significance file's gene symbol - so it gets a throwaway in-memory Genes table
    rather than the real database, keeping the comparison reader-only. Its frame
    carries no feature-type column (LeafCutter reports introns and nothing else),
    which normalize_junctions_frame fills in - the same boundary analyze_junctions
    puts it through.
    """
    sig_path = os.path.join(tmp, 'leafcutter_ds_cluster_significance.txt')
    effect_path = os.path.join(tmp, 'leafcutter_ds_effect_sizes.txt')
    with open(sig_path, 'w') as handle:
        handle.write("cluster\tstatus\tloglr\tdf\tp\tp.adjust\tgenes\n")
        handle.write(f"chr{chrom}:clu_1_{strand}\tSuccess\t4.2\t2\t0.001\t0.01\tGENEX\n")
    with open(effect_path, 'w') as handle:
        handle.write("intron\tlogef\tGENEX\tdeltapsi\n")
        for start, end in introns:
            handle.write(f"chr{chrom}:{start}:{end}:clu_1_{strand}\t0.03\t0.4\t0.01\n")

    con = sqlite3.connect(':memory:')
    try:
        con.execute("CREATE TABLE Genes (gene_GeneID_id TEXT, gene_ensembl_id TEXT, "
                    "gene_symbol TEXT, specie TEXT)")
        con.execute("INSERT INTO Genes VALUES ('1', 'ENSG00000000001', 'GENEX', ?)",
                    (utils.SPECIE_DB_NAME['human'],))
        df = alternative_splicing.leafcutter_read_input_files(
            con, significance_file=sig_path, effect_sizes_file=effect_path, specie='human')
    finally:
        con.close()
    return utils.normalize_junctions_frame(df, specie='human')


def _features(df):
    """The comparable part of a junctions frame: what each feature is, not what the
    tool called it or how the cluster was named."""
    return set(zip(df['start_position'], df['end_position'], df[utils.FEATURE_TYPE_COLUMN]))


J = utils.FEATURE_JUNCTION
R = utils.FEATURE_RETAINED_INTRON


# --- the events, stated independently of any reader ------------------------
#
# Skipped exon: upstream exon ends at E1, skipped exon spans S2..E2, downstream
# exon starts at S3. Introns: E1-S2 (inclusion, upstream), E2-S3 (inclusion,
# downstream), E1-S3 (skipping).
E1, S2, E2, S3 = 1000, 2000, 2100, 3000
SE_FEATURES = {(E1, S2, J), (E2, S3, J), (E1, S3, J)}

# Retained intron: one intron, present as a junction in the spliced isoform and
# contained in an exon in the retained one.
RI_START, RI_END = 1500, 2500
RI_FEATURES = {(RI_START, RI_END, J), (RI_START, RI_END, R)}


@pytest.mark.parametrize('strand', ['+', '-'])
def test_skipped_exon_agrees_between_rmats_and_majiq(strand):
    """rMATS names the exons by genomic position for SE, and MAJIQ lists the
    junctions outright, so both must give the same three introns on either strand."""
    with tempfile.TemporaryDirectory() as tmp:
        rmats = _features(_rmats_dir(tmp, 'SE', {
            'upstreamES': 500, 'upstreamEE': E1,
            'exonStart_0base': S2, 'exonEnd': E2,
            'downstreamES': S3, 'downstreamEE': 3500}, strand))
    with tempfile.TemporaryDirectory() as tmp:
        majiq = _features(_majiq(tmp, f'{E1}-{S2};{E2}-{S3};{E1}-{S3}', strand=strand))

    assert rmats == SE_FEATURES, f"rMATS SE ({strand}): {sorted(rmats)}"
    assert majiq == SE_FEATURES, f"MAJIQ SE ({strand}): {sorted(majiq)}"
    assert rmats == majiq


@pytest.mark.parametrize('strand', ['+', '-'])
def test_retained_intron_agrees_between_rmats_and_suppa(strand):
    """Both must emit the interval twice - once to be found spliced out, once to be
    found retained - or the retaining isoform is invisible."""
    with tempfile.TemporaryDirectory() as tmp:
        rmats = _features(_rmats_dir(tmp, 'RI', {
            'riExonStart_0base': 900, 'riExonEnd': 3100,
            'upstreamES': 900, 'upstreamEE': RI_START,
            'downstreamES': RI_END, 'downstreamEE': 3100}, strand))
    with tempfile.TemporaryDirectory() as tmp:
        suppa = _features(_suppa(tmp, f'ENSG00000000001;RI:21:900:{RI_START}-{RI_END}:3100:{strand}'))

    assert rmats == RI_FEATURES, f"rMATS RI ({strand}): {sorted(rmats)}"
    assert suppa == RI_FEATURES, f"SUPPA RI ({strand}): {sorted(suppa)}"
    assert rmats == suppa


@pytest.mark.parametrize('strand', ['+', '-'])
def test_alt_splice_site_yields_two_distinct_junctions_in_both_formats(strand):
    """The alternative-splice-site types are where the formats differ most: SUPPA
    writes both junctions explicitly, while rMATS names the long and short forms of
    one exon and leaves the reader to work out which boundary varies - which is
    strand-dependent. Reading it wrongly collapsed the event to a single junction,
    and this is the check that would have caught it.
    """
    # Minus-strand geometry: long and short forms share their genomic END and differ
    # at their START; the flanking exon lies below. Plus-strand is the mirror.
    if strand == '-':
        columns = {'longExonStart_0base': 1000, 'longExonEnd': 1500,
                   'shortES': 1100, 'shortEE': 1500,
                   'flankingES': 500, 'flankingEE': 800}
        expected = {(800, 1000, J), (800, 1100, J)}
    else:
        columns = {'longExonStart_0base': 1000, 'longExonEnd': 1500,
                   'shortES': 1000, 'shortEE': 1400,
                   'flankingES': 2000, 'flankingEE': 2500}
        expected = {(1500, 2000, J), (1400, 2000, J)}

    with tempfile.TemporaryDirectory() as tmp:
        rmats = _features(_rmats_dir(tmp, 'A5SS', columns, strand))

    assert len({(s, e) for s, e, _ in rmats}) == 2, \
        f"A5SS ({strand}) collapsed to one junction: {sorted(rmats)}"
    assert rmats == expected, f"A5SS ({strand}): {sorted(rmats)}"


@pytest.mark.parametrize('strand', ['+', '-'])
def test_skipped_exon_agrees_between_rmats_and_leafcutter(strand):
    """LeafCutter names introns directly rather than exons, so it is the one format
    that states the answer in the same terms the event is defined in - which makes it
    the check that the exon-based readers put their boundaries on the right side."""
    with tempfile.TemporaryDirectory() as tmp:
        rmats = _features(_rmats_dir(tmp, 'SE', {
            'upstreamES': 500, 'upstreamEE': E1,
            'exonStart_0base': S2, 'exonEnd': E2,
            'downstreamES': S3, 'downstreamEE': 3500}, strand))
    with tempfile.TemporaryDirectory() as tmp:
        leafcutter = _features(_leafcutter(tmp, [(E1, S2), (E2, S3), (E1, S3)], strand=strand))

    assert leafcutter == SE_FEATURES, f"LeafCutter SE ({strand}): {sorted(leafcutter)}"
    assert rmats == leafcutter, f"rMATS {sorted(rmats)} != LeafCutter {sorted(leafcutter)}"


@pytest.mark.parametrize('strand', ['+', '-'])
def test_alt_splice_site_agrees_between_rmats_and_leafcutter(strand):
    """The A5SS geometry that rMATS states as long/short exon forms, stated by
    LeafCutter as the two introns it actually is. Same expectation as the rMATS-only
    case above - written here from the other side of the translation."""
    if strand == '-':
        columns = {'longExonStart_0base': 1000, 'longExonEnd': 1500,
                   'shortES': 1100, 'shortEE': 1500,
                   'flankingES': 500, 'flankingEE': 800}
        introns = [(800, 1000), (800, 1100)]
    else:
        columns = {'longExonStart_0base': 1000, 'longExonEnd': 1500,
                   'shortES': 1000, 'shortEE': 1400,
                   'flankingES': 2000, 'flankingEE': 2500}
        introns = [(1500, 2000), (1400, 2000)]
    expected = {(s, e, J) for s, e in introns}

    with tempfile.TemporaryDirectory() as tmp:
        rmats = _features(_rmats_dir(tmp, 'A5SS', columns, strand))
    with tempfile.TemporaryDirectory() as tmp:
        leafcutter = _features(_leafcutter(tmp, introns, strand=strand))

    assert leafcutter == expected, f"LeafCutter A5SS ({strand}): {sorted(leafcutter)}"
    assert rmats == leafcutter, f"rMATS {sorted(rmats)} != LeafCutter {sorted(leafcutter)}"


@pytest.mark.parametrize('strand', ['+', '-'])
def test_skipped_exon_agrees_between_rmats_and_suppa(strand):
    """SUPPA used to emit 2 junctions here where rMATS emits 3: its event id names
    both inclusion junctions, but the downstream one supplied only its right
    coordinate to the constructed skipping junction and was never added in its own
    right. A transcript carrying only that downstream junction was comparable when
    the event came from rMATS and not when it came from SUPPA."""
    with tempfile.TemporaryDirectory() as tmp:
        rmats = _features(_rmats_dir(tmp, 'SE', {
            'upstreamES': 500, 'upstreamEE': E1,
            'exonStart_0base': S2, 'exonEnd': E2,
            'downstreamES': S3, 'downstreamEE': 3500}, strand))
    with tempfile.TemporaryDirectory() as tmp:
        suppa = _features(_suppa(tmp, f'ENSG00000000001;SE:21:{E1}-{S2}:{E2}-{S3}:{strand}'))

    assert suppa == SE_FEATURES, f"SUPPA SE ({strand}): {sorted(suppa)}"
    assert rmats == suppa, f"rMATS {sorted(rmats)} != SUPPA {sorted(suppa)}"
