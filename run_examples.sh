#!/usr/bin/env bash
#
# One DOMAS run per input format, against the fixtures in tests/.
# The DoChaP database is not in this repository - download it and pass its path:
#
#     ./run_examples.sh /path/to/DB_merged.sqlite [output_dir]
#
# Every command below is complete: copy one, swap in your own input, and run it.
# -species is required because none of these formats carries a species field.
# -max_clusters caps the three large fixtures so an example finishes in seconds.
# DOMAS creates the output directory, and writes the statistics report only with -stats.
#
# -no_excel is here for the comparison at the bottom of this file, which reads the
# stored CSV references. Drop it from a command you copy and the results are saved
# as .xlsx instead, which is DOMAS's default.
set -e

DOCHAP=$1
OUT=${2:-run_examples_output}

python3 code/domas.py -format ioe -input tests/ioe -species human -dochap "$DOCHAP" -output_csv "$OUT/ioe.csv" -no_excel

python3 code/domas.py -format rmats -input tests/rmats -species human -dochap "$DOCHAP" -output_csv "$OUT/rmats.csv" -max_clusters 200 -no_excel

python3 code/domas.py -format majiq -input tests/majiq/NveB_Mono_voila.txt -species human -dochap "$DOCHAP" -output_csv "$OUT/majiq.csv" -max_clusters 200 -no_excel

python3 code/domas.py -format leafcutter -lc_sig tests/leafcutter/leafcutter_ds_cluster_significance.txt -lc_effect tests/leafcutter/leafcutter_ds_effect_sizes.txt -species human -dochap "$DOCHAP" -output_csv "$OUT/leafcutter.csv" -max_clusters 200 -no_excel

# Each result should equal the stored reference:
for name in ioe rmats majiq leafcutter; do
    python3 tests/compare_run_example.py "$OUT/$name.csv" "tests/run_examples/$name.csv"
done
