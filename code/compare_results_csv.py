"""
Compare two junction-analysis results.csv files for exact equality, ignoring
row order. Built for comparing a full-scale re-run's output against a prior
baseline (e.g. tests/manual_review/ioe_full_representative/results.csv) after
a performance-only refactor of junction_analisys.py, where the two files can
be tens of millions of rows / a couple GB - too large to comfortably
pd.read_csv() + pd.concat() the naive way (see the OOM issues that approach
caused for results_stats.py on the same kind of file).

Usage:
    python3 compare_results_csv.py <old.csv> <new.csv>

Exit code 0 and "MATCH" printed if identical (order-independent); exit code 1
and a sample of differing rows printed otherwise.
"""
import sys

import numpy as np
import pandas as pd

# Low-cardinality columns, read as category (and the length/count pair as
# float32): a 20M-row file then costs a few hundred MB instead of several GB.
CATEGORY_COLS = [
    'event', 'gene_symbol', 'specie', 'event_type', 'canonical_transcript_id',
    'transcript_id', 'domain_name', 'is_longest_cds', 'is_most_like_canonical',
]
FLOAT_COLS = ['c_domain_length', 't_domain_length', 'c_domains_number', 't_domains_number']


def _read_compact(path, chunk_rows=1_000_000):
    header_cols = list(pd.read_csv(path, nrows=0).columns)
    dtype = {c: 'category' for c in CATEGORY_COLS if c in header_cols}
    dtype.update({c: 'float32' for c in FLOAT_COLS if c in header_cols})

    chunks = list(pd.read_csv(path, dtype=dtype, chunksize=chunk_rows))
    if len(chunks) == 1:
        return chunks[0]

    # Union the per-chunk categories explicitly: plain concat upcasts them back
    # to object whenever the chunks disagree, which they always do.
    cat_cols = [c for c in chunks[0].columns if isinstance(chunks[0][c].dtype, pd.CategoricalDtype)]
    other_cols = [c for c in chunks[0].columns if c not in cat_cols]
    data = {}
    for col in cat_cols:
        data[col] = pd.api.types.union_categoricals([chunk[col] for chunk in chunks], sort_categories=False)
    for col in other_cols:
        data[col] = np.concatenate([chunk[col].to_numpy() for chunk in chunks])
    return pd.DataFrame(data, columns=list(chunks[0].columns))


def main():
    if len(sys.argv) != 3:
        print(f"Usage: python3 {sys.argv[0]} <old.csv> <new.csv>")
        sys.exit(2)

    old_path, new_path = sys.argv[1], sys.argv[2]
    print(f"Reading {old_path} ...")
    df_old = _read_compact(old_path)
    print(f"  {len(df_old):,} rows")
    print(f"Reading {new_path} ...")
    df_new = _read_compact(new_path)
    print(f"  {len(df_new):,} rows")
    # NaN is normalised by the .astype(str) below (-> 'nan' in both frames);
    # fillna() is not an option, category columns reject an unknown value.

    if set(df_old.columns) != set(df_new.columns):
        print("COLUMN MISMATCH")
        print("  only in old:", sorted(set(df_old.columns) - set(df_new.columns)))
        print("  only in new:", sorted(set(df_new.columns) - set(df_old.columns)))
        sys.exit(1)

    sort_cols = list(df_old.columns)
    print("Sorting both for order-independent comparison...")
    # Sort on the string form. Unordered categories sort by code, assigned in
    # first-appearance order, so the two files could misalign on equal values.
    df_old = df_old.loc[df_old[sort_cols].astype(str).sort_values(sort_cols).index].reset_index(drop=True)
    df_new = df_new.loc[df_new[sort_cols].astype(str).sort_values(sort_cols).index].reset_index(drop=True)

    if len(df_old) != len(df_new):
        print(f"ROW COUNT MISMATCH: old={len(df_old):,} new={len(df_new):,}")
        sys.exit(1)

    # Compare as strings so a dtype difference between the two loads cannot
    # register as a mismatch; only the values matter.
    mismatch_mask = np.zeros(len(df_old), dtype=bool)
    for col in sort_cols:
        mismatch_mask |= (df_old[col].astype(str).to_numpy() != df_new[col].astype(str).to_numpy())

    n_mismatch = int(mismatch_mask.sum())
    if n_mismatch == 0:
        print(f"MATCH: {len(df_old):,} rows identical (order-independent).")
        sys.exit(0)

    print(f"MISMATCH: {n_mismatch:,} / {len(df_old):,} rows differ after sorting.")
    print("\nSample of differing rows (old vs new):")
    sample_idx = np.where(mismatch_mask)[0][:10]
    with pd.option_context('display.max_columns', None, 'display.width', 200):
        for i in sample_idx:
            print(f"\n--- row {i} ---")
            print("old:", df_old.iloc[i].to_dict())
            print("new:", df_new.iloc[i].to_dict())
    sys.exit(1)


if __name__ == '__main__':
    main()
