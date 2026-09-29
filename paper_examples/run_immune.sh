#!/usr/bin/env bash
#
# Use case 1 of the DOMAS paper: differential splicing between human monocytes and
# CD4 T cells, dataset GSE60424 (Linsley et al., 2014), called with LeafCutter.
#
# The DoChaP database is not in this repository - download it and pass its path:
#
#     ./run_immune.sh /path/to/DB_merged.sqlite [output_dir]
#
# -max_clusters 100 is what the paper reports and what the web server applies: the
#   100 events with the lowest p.adjust, not the first 100 in the file.
# -no_excel keeps the results as CSV, so they can be read and diffed as text. Drop
#   it and DOMAS writes .xlsx instead - its default - with every gene symbol linked
#   to that gene's DoChaP page.

set -e

DOCHAP=$1
OUT=${2:-results/immune}

if [ -z "$DOCHAP" ]; then
    echo "usage: $0 /path/to/DB_merged.sqlite [output_dir]" >&2
    exit 1
fi

# Run from this directory, so the paths the run summary reports are the ones in
# this repository rather than wherever the caller happened to be standing. A
# relative output_dir is therefore relative to paper_examples/.
cd "$(dirname "$0")"

python3 ../code/domas.py \
    -format leafcutter \
    -lc_sig    input/immune_leafcutter_significance.txt \
    -lc_effect input/immune_leafcutter_effect_sizes.txt \
    -species human \
    -dochap "$DOCHAP" \
    -max_clusters 100 \
    -output_csv "$OUT/compared.csv" \
    -no_excel \
    -num_workers 1

