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
import sys
import sqlite3
import tempfile
import warnings

warnings.filterwarnings('ignore')

REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(REPO, 'code'))

import pandas as pd  # noqa: E402
import utils  # noqa: E402
from junction_analisys import JunctionsAnalysis  # noqa: E402
from alternative_splicing import leafcutter_read_input_files  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from conftest import DEFAULT_DB_PATH as DB  # noqa: E402
OUT_ROOT = os.path.join(REPO, 'review_cases')


def head_rows(path, n=3):
    """First n data rows of a tab file - keeps rmats parsing cheap."""
    return pd.read_csv(path, sep='\t', nrows=n)


def rmats_subset(tmp_dir):
    """A couple of events per rMATS type, written to a small directory so the
    parser does not have to walk all 98k clusters."""
    os.makedirs(tmp_dir, exist_ok=True)
    for event_type, filename in utils._RMATS_EVENT_FILES.items():
        src = os.path.join(REPO, 'tests', 'rmats', filename)
        if not os.path.exists(src):
            continue
        df = pd.read_csv(src, sep='\t')
        # one plus-strand and one minus-strand event per type, so the
        # strand-dependent A5SS/A3SS construction is visible
        picks = [g.head(1) for _, g in df.groupby('strand')]
        pd.concat(picks).to_csv(os.path.join(tmp_dir, filename), sep='\t', index=False)
    return utils.rmats2junctions(tmp_dir)


def run(name, df_junctions, con, specie='human'):
    out = os.path.join(OUT_ROOT, name)
    os.makedirs(out, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(out)
    try:
        JunctionsAnalysis(con).analyze_junctions(
            df_junctions=df_junctions,
            output_path=os.path.join(out, 'results.csv'),
            specie=specie, create_pdf=True, num_workers=1,
            use_representative_domains=True)
    finally:
        os.chdir(cwd)
    pdfs = [f for f in os.listdir(out) if f.endswith('.pdf')]
    print(f'{name}: {df_junctions.cluster_name.nunique()} clusters -> {len(pdfs)} PDFs', flush=True)


def main():
    con = sqlite3.connect(DB)

    # The cut-down rMATS directory is scratch, so it goes to a temp dir rather than
    # being left behind inside review_cases/ looking like one of the case sets.
    with tempfile.TemporaryDirectory() as tmp:
        run('rmats', rmats_subset(tmp), con)

    majiq = utils.voila2junctions(os.path.join(REPO, 'tests/majiq/NveB_Mono_voila.txt'))
    keep = sorted(majiq.cluster_name.unique())[:4]
    run('majiq', majiq[majiq.cluster_name.isin(keep)].copy(), con)

    lc = leafcutter_read_input_files(
        con,
        os.path.join(REPO, 'tests/leafcutter/leafcutter_ds_cluster_significance.txt'),
        os.path.join(REPO, 'tests/leafcutter/leafcutter_ds_effect_sizes.txt'))
    # include the multi-gene cluster explicitly - one event, two genes
    multi = [c for c in lc.cluster_name.unique() if c.startswith('chr3:clu_10461')]
    # ... and the only reviewable source of `longer_domains`, the one Table S4
    # outcome absent from every other fixture. It occurs in 83 clusters of the full
    # leafcutter run, none with fewer than 11 transcripts; this one has just a single
    # comparable transcript, so there is one comparison to follow.
    longer_domains = [c for c in lc.cluster_name.unique() if c.startswith('chr5:clu_3349')]
    keep = sorted(lc.cluster_name.unique())[:4] + multi + longer_domains
    run('leafcutter', lc[lc.cluster_name.isin(keep)].copy(), con)

    con.close()


if __name__ == '__main__':
    main()
