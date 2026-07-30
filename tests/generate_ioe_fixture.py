"""Regenerate tests/ioe_example_junctions.csv from the real SUPPA .ioe files.

The fixture is a *sample* of real SUPPA events, frozen so the golden outputs have a
stable input. It used to be a by-product of analyze_ioe_files(), which writes the
same file as a side effect - which left no record of how it had been produced, and
no way to refresh it when the reader or the database moved on. Hence this script.

Two details that are not obvious and were lost when the original was made:

  * the `num_transcripts >= 2` filter. keep_min_transcript_clusters() keeps the
    clusters whose gene has the FEWEST transcripts, and analyze_ioe_files() has this
    filter commented out, so regenerating without it yields single-transcript genes -
    clusters that can only ever report `only_one_transcript`. The original fixture's
    num_transcripts starts at 2, so the filter was plainly active when it was built.
  * the species. These are the H_sapiens files, matching the default -ioe_pattern.

Run from the repository root:

    python3 tests/generate_ioe_fixture.py
"""
import os
import re
import sqlite3
import sys

import pandas as pd

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(TESTS_DIR, '..'))
sys.path.insert(0, os.path.join(REPO, 'code'))

import utils  # noqa: E402
from alternative_splicing import keep_min_transcript_clusters  # noqa: E402
from conftest import DEFAULT_DB_PATH  # noqa: E402

IOE_DIR = os.path.join(REPO, 'external_data', 'H_sapiens')
IOE_PATTERN = r"output_prefix_.*_strict.ioe"
OUT = os.path.join(TESTS_DIR, 'ioe_example_junctions.csv')
EXAMPLES_PER_EVENT = 5
MIN_TRANSCRIPTS = 2


def main():
    if not os.path.isdir(IOE_DIR):
        raise SystemExit(f"SUPPA input not found at {IOE_DIR}")
    con = sqlite3.connect(DEFAULT_DB_PATH)

    frames = [utils.ioe2junctions(os.path.join(IOE_DIR, name))
              for name in sorted(os.listdir(IOE_DIR)) if re.match(IOE_PATTERN, name)]
    if not frames:
        raise SystemExit(f"no files matching {IOE_PATTERN} in {IOE_DIR}")
    df = pd.concat(frames, ignore_index=True)

    gene_ids = df.gene_ensembl_id.unique().tolist()
    df['gene_symbol'] = df.gene_ensembl_id.map(utils.get_gene_symbols(con, gene_ids))
    df['num_transcripts'] = df.gene_ensembl_id.map(
        utils.get_genes_number_of_transcripts(con, gene_ids))
    df['num_canonical_exons'] = df.gene_ensembl_id.map(
        utils.get_canonical_exon_counts(con, gene_ids))

    # See the module docstring: without this the sample is all single-transcript genes.
    df = df[df.num_transcripts >= MIN_TRANSCRIPTS]

    examples = keep_min_transcript_clusters(df, examples_per_event=EXAMPLES_PER_EVENT)
    examples.to_csv(OUT, index=False)

    by_type = examples.assign(t=examples.cluster_name.str.split('_').str[0]) \
                      .groupby('t').cluster_name.nunique().to_dict()
    print(f"wrote {OUT}")
    print(f"  {len(examples)} rows, {examples.cluster_name.nunique()} clusters")
    print(f"  per event type: {by_type}")
    print(f"  num_transcripts {examples.num_transcripts.min():.0f}"
          f"-{examples.num_transcripts.max():.0f}")
    con.close()


if __name__ == '__main__':
    main()
