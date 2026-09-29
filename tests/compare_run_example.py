"""Compare a run_examples.sh result against its stored reference.

Row order is not part of the result - clusters are analysed in parallel and
returned as they finish - so both sides are sorted before comparing. Exits 0 on
a match, 1 otherwise, printing the first few differing rows.

    python3 tests/compare_run_example.py <results.csv> <reference.csv>
"""
import sys

import pandas as pd


def load(path):
    df = pd.read_csv(path).fillna('')
    return df.sort_values(list(df.columns)).reset_index(drop=True)


def main():
    if len(sys.argv) != 3:
        print(__doc__.strip(), file=sys.stderr)
        return 2

    generated, reference = load(sys.argv[1]), load(sys.argv[2])

    if list(generated.columns) != list(reference.columns):
        print(f"  columns differ:\n    got      {list(generated.columns)}\n"
              f"    expected {list(reference.columns)}", file=sys.stderr)
        return 1

    if len(generated) != len(reference):
        print(f"  row count differs: got {len(generated)}, expected {len(reference)}",
              file=sys.stderr)
        return 1

    unequal = (generated.astype(str) != reference.astype(str)).any(axis=1)
    if unequal.any():
        print(f"  {int(unequal.sum())} of {len(generated)} rows differ; first few:",
              file=sys.stderr)
        print(generated[unequal].head(3).to_string(), file=sys.stderr)
        return 1

    print(f"{sys.argv[1]}: {len(generated)} rows, matches {sys.argv[2]}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
