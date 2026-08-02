"""Generate the review_cases PDFs for the three formats whose tests run with
create_pdf=False, so every input format has something that can be looked at.

The PDFs themselves are not committed - .gitignore excludes *.pdf across the repo,
and they are 5 MB of regenerable output - so this script is what makes the set
reproducible. See review_cases/INDEX.md for what each case demonstrates.

Run from anywhere:

    python3 tests/generate_review_pdfs.py

Must be a real file rather than a heredoc: ProcessPoolExecutor spawns workers that
re-import __main__, which fails for '<stdin>'.
"""
import os
import shutil
import sys
import sqlite3
import tempfile
import warnings

warnings.filterwarnings('ignore')

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(REPO, 'code'))

import pandas as pd  # noqa: E402
import utils  # noqa: E402
from junction_analisys import JunctionsAnalysis, NON_COMPARISON_EVENTS  # noqa: E402
from alternative_splicing import leafcutter_read_input_files  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from conftest import DEFAULT_DB_PATH as DB  # noqa: E402
OUT_ROOT = os.path.join(REPO, 'review_cases')


def head_rows(path, n=3):
    """First n data rows of a tab file - keeps rmats parsing cheap."""
    return pd.read_csv(path, sep='\t', nrows=n)


def _clusters(df, *names):
    """Every cluster of `df` named by one of `names`.

    A LeafCutter cluster naming more than one gene is analysed once per gene and
    carries a `:SYMBOL` suffix, so an exact match is not enough; a prefix match
    alone would also catch chr14:clu_94280, hence the explicit separator.
    """
    return [c for c in df.cluster_name.unique()
            if any(c == name or c.startswith(name + ':') for name in names)]


# (event type, strand) pairs whose first event resolves no comparable transcript,
# so several candidates are read and the one that compares something is drawn - see
# comparable_by_cluster(). SE + is deliberately absent: the fixture holds 3
# plus-strand SE clusters and the only one that compares anything reaches no domain,
# so D3a stays the documented gap.
WIDEN = {('A3SS', '-'), ('A5SS', '+'), ('SE', '-'), ('MXE', '-')}
CANDIDATE_ROWS = 30
MAJIQ_CANDIDATES = 24


def gene_strands(con):
    """gene_ensembl_id -> strand, for grouping rMATS clusters by the strand whose
    handling they exercise. The junctions frame does not carry it - the analysis
    reads it from the database per gene."""
    return dict(pd.read_sql_query(
        "SELECT gene_ensembl_id, strand FROM Genes WHERE specie = 'H_sapiens'", con).values)


def rmats_subset(tmp_dir):
    """A couple of events per rMATS type, written to a small directory so the
    parser does not have to walk all 98k clusters.

    One event per strand is enough where that event resolves a comparison; for the
    pairs in WIDEN it does not, so more candidates are read and the choice is made
    on the outcome instead of on position in the file.
    """
    os.makedirs(tmp_dir, exist_ok=True)
    for event_type, filename in utils._RMATS_EVENT_FILES.items():
        src = os.path.join(REPO, 'tests', 'rmats', filename)
        if not os.path.exists(src):
            continue
        df = pd.read_csv(src, sep='\t')
        # one plus-strand and one minus-strand event per type, so the
        # strand-dependent A5SS/A3SS construction is visible
        picks = [g.head(CANDIDATE_ROWS if (event_type, strand) in WIDEN else 1)
                 for strand, g in df.groupby('strand')]
        pd.concat(picks).to_csv(os.path.join(tmp_dir, filename), sep='\t', index=False)
    return utils.rmats2junctions(tmp_dir)


def comparable_by_cluster(con, df_junctions, specie='human'):
    """{cluster name: {n_comparable, has_domain_call}}.

    A first pass that draws nothing, so the clusters worth drawing are chosen by
    what they resolve rather than by their position in the input file.
    """
    with tempfile.TemporaryDirectory() as tmp:
        results = JunctionsAnalysis(con).analyze_junctions(
            df_junctions=df_junctions, output_path=os.path.join(tmp, 'results.csv'),
            specie=specie, create_pdf=False, num_workers=1,
            use_representative_domains=True,
            write_all_comparable=True,
        )
    out = {}
    for result in results:
        df = result.get_results_df()
        compared = df[~df.event.isin(NON_COMPARISON_EVENTS)]
        out[result.cluster_name] = {
            'n_comparable': compared.transcript_id.dropna().nunique(),
            # A comparison reaching no domain beats none, but shows no domain
            # call - so it wins only where nothing else exists.
            'has_domain_call': bool(len(compared[compared.event != 'no_domains_in_region'])),
        }
    return out


def best_per_group(scores, group_of, sizes):
    """One cluster per group: a real domain call first, then most transcripts
    compared, then the smallest gene, then the name - deterministic and quick to
    read. A cluster whose group is unknown (a gene with no strand in the database)
    is skipped rather than forming a group of its own."""
    best = {}
    for cluster, score in scores.items():
        group = group_of(cluster)
        if group is None or None in group:
            continue
        key = (score['has_domain_call'], score['n_comparable'],
               -sizes.get(cluster, 0), str(cluster))
        if group not in best or key > best[group][0]:
            best[group] = (key, cluster)
    return [cluster for _, cluster in best.values()]


def run(name, df_junctions, con, specie='human', restrict=True):
    """Analyse `df_junctions` and write its PDFs to review_cases/<name>/.

    restrict: draw only the canonical transcript and the ones compared to it, which
    is what almost every case wants - a figure about intron retention or the domain
    ladder should not bury its two relevant tracks among twenty. The A series is the
    exception: its subject is which candidate was chosen, so they must be visible.
    """
    out = os.path.join(OUT_ROOT, name)
    # Cleared, not merely created: the counter in each file name is the cluster's
    # position in the run, so changing the kept set renumbers everything and files
    # from an earlier run would survive under different names. collect_review_pdfs.py
    # resolves a case with sorted(glob)[0] and would pick whichever sorts first.
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(out)
    try:
        JunctionsAnalysis(con).analyze_junctions(
            df_junctions=df_junctions,
            output_path=os.path.join(out, 'results.csv'),
            specie=specie, create_pdf=True, num_workers=1,
            restrict_pdf_to_comparable=restrict,
            use_representative_domains=True,
            write_all_comparable=True,
        )
    finally:
        os.chdir(cwd)
    pdfs = [f for f in os.listdir(out) if f.endswith('.pdf')]
    print(f'{name}: {df_junctions.cluster_name.nunique()} clusters -> {len(pdfs)} PDFs', flush=True)


def main():
    con = sqlite3.connect(DB)

    # The cut-down rMATS directory is scratch, so it goes to a temp dir rather than
    # being left behind inside review_cases/ looking like one of the case sets.
    with tempfile.TemporaryDirectory() as tmp:
        rmats = rmats_subset(tmp)
    strand_of = gene_strands(con)
    sizes = rmats.groupby('cluster_name').size().to_dict()

    def type_and_strand(cluster):
        gene = str(cluster).split('_')[1]
        return str(cluster).split('_')[0], strand_of.get(gene)

    keep = best_per_group(comparable_by_cluster(con, rmats), type_and_strand, sizes)
    run('rmats', rmats[rmats.cluster_name.isin(keep)].copy(), con)

    majiq = utils.voila2junctions(os.path.join(REPO, 'tests/majiq/NveB_Mono_voila.txt'))
    # A wider slice than the four drawn, for the same reason as WIDEN: the first
    # four by name include LSVs that resolve no comparable transcript.
    candidates = sorted(majiq.cluster_name.unique())[:MAJIQ_CANDIDATES]
    scored = comparable_by_cluster(con, majiq[majiq.cluster_name.isin(candidates)].copy())
    keep = [c for c, _ in sorted(
        scored.items(),
        key=lambda kv: (not kv[1]['has_domain_call'], -kv[1]['n_comparable'], str(kv[0])))[:4]]
    run('majiq', majiq[majiq.cluster_name.isin(keep)].copy(), con)

    lc = leafcutter_read_input_files(
        con,
        os.path.join(REPO, 'tests/leafcutter/leafcutter_ds_cluster_significance.txt'),
        os.path.join(REPO, 'tests/leafcutter/leafcutter_ds_effect_sizes.txt'))
    # a multi-gene cluster - one event, two genes, analysed once per gene. Both
    # halves of chr15:clu_17411 (HYPK/SERF2) compare something, which the second-gene
    # case needs to show anything at all. 24 clusters of the full run qualify.
    multi = _clusters(lc, 'chr15:clu_17411')
    # ... and the only reviewable source of `longer_domains`, the one Table S4
    # outcome absent from every other fixture. Of the 83 clusters carrying it, this
    # is the one with a single comparable transcript - one comparison to follow.
    longer_domains = [c for c in lc.cluster_name.unique() if c.startswith('chr5:clu_3349')]
    # ... and the clusters the A series needs. It wants a gene with 3-5 transcripts
    # and more than one comparable candidate, which no reference fixture holds: the
    # ioe one is 31 two-transcript genes, the category one jumps from 2 to 6. These
    # three carry a real choice (generate_review_index.py derives which is which):
    #   chr14:clu_9428  TMED8   the two rules land on different transcripts
    #   chr12:clu_6316  GLIPR1  three coding candidates, the gate rejects all of them
    #   chr20:clu_4032  FAM209A a non-coding candidate holds the longest CDS
    a_series = _clusters(lc, 'chr14:clu_9428', 'chr12:clu_6316', 'chr20:clu_4032')
    # ... and the clusters showing the duplicate-collapse step of
    # filter_representative_domains(), which nothing else in the review set reaches.
    # It removes an entry sharing an accession with a longer kept one that it
    # overlaps by as little as ONE residue, collapsing adjacent repeat instances:
    #   chr3:clu_2426    TOPBP1   IPR036420 102-199 + 199-302, 1 aa
    #   chr19:clu_15038  ZNF317   IPR036236 x3 abutting by 2 aa, two discarded
    #   chrX:clu_319     DNASE1L1 IPR016202 2-273 vs 3-274 - 271 aa, the case where
    #                             collapsing IS right, for contrast
    # SORL1 (IPR003961, 5 aa) is already drawn by the category_examples reference.
    ladder_dup = _clusters(lc, 'chr3:clu_2426', 'chr19:clu_15038', 'chrX:clu_319')
    keep = (sorted(lc.cluster_name.unique())[:4] + multi + longer_domains
            + a_series + ladder_dup)
    run('leafcutter', lc[lc.cluster_name.isin(keep)].copy(), con)
    # The A clusters again, every transcript drawn. A second directory rather than
    # switching the run above: other cases in this fixture (E1, B6/B7) read better
    # restricted, and the flag is per run.
    run('leafcutter_all_transcripts', lc[lc.cluster_name.isin(a_series)].copy(),
        con, restrict=False)

    con.close()


if __name__ == '__main__':
    main()
