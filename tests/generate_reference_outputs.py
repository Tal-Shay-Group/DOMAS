"""
Run JunctionsAnalysis.analyze_junctions() for every combination of
use_longest_cds / restrict_pdf_to_comparable, against both fixture files,
and save the resulting CSV + PDFs into persistent, descriptively-named
subdirectories under tests/reference_outputs/ for manual inspection.

These outputs are meant to be reviewed by hand once, then used as a
reference point for test_flags.py going forward.
"""
import os
import shutil
import sqlite3
import sys

import matplotlib
matplotlib.use('Agg')

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.normpath(os.path.join(TESTS_DIR, '..', 'code'))
sys.path.insert(0, CODE_DIR)

from junction_analisys import JunctionsAnalysis  # noqa: E402
from alternative_splicing import hadas_read_input_file, read_junctions_csv  # noqa: E402
from pdf_text_utils import write_pdf_text_manifest  # noqa: E402

DB_PATH = '/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-web/DB_merged.sqlite'
IOE_CSV = os.path.join(TESTS_DIR, 'ioe_example_junctions.csv')
HADAS_XLSX = os.path.join(TESTS_DIR, 'short_H_vs_M_HN6.xlsx')
OUTPUT_ROOT = os.path.join(TESTS_DIR, 'reference_outputs')

FLAG_COMBINATIONS = [
    (False, False),
    (False, True),
    (True, False),
    (True, True),
]

INPUT_FILES = [
    ('ioe_csv', IOE_CSV, False),
    ('hadas_xlsx', HADAS_XLSX, True),
]


def run_one(con, label, junctions_csv, hadas_format, use_longest_cds, restrict_pdf_to_comparable):
    case_name = f"{label}__longestcds_{use_longest_cds}__restrict_{restrict_pdf_to_comparable}"
    case_dir = os.path.join(OUTPUT_ROOT, case_name)
    if os.path.exists(case_dir):
        shutil.rmtree(case_dir)
    os.makedirs(case_dir)

    df_junctions = hadas_read_input_file(con, junctions_csv) if hadas_format else read_junctions_csv(junctions_csv)
    output_path = os.path.join(case_dir, 'results.csv')
    cwd_before = os.getcwd()
    os.chdir(case_dir)
    try:
        analysis = JunctionsAnalysis(con)
        results = analysis.analyze_junctions(
            df_junctions=df_junctions,
            output_path=output_path,
            create_pdf=True,
            num_workers=1,
            use_longest_cds=use_longest_cds,
            restrict_pdf_to_comparable=restrict_pdf_to_comparable,
        )
    finally:
        os.chdir(cwd_before)

    num_pdfs = len([f for f in os.listdir(case_dir) if f.endswith('.pdf')])
    write_pdf_text_manifest(case_dir)
    print(f"  {case_name}: {len(results)} clusters, {num_pdfs} PDFs -> {case_dir}")
    return case_dir


def main():
    if not os.path.exists(DB_PATH):
        print(f"ERROR: DoChaP database not found at {DB_PATH}")
        sys.exit(1)

    if os.path.exists(OUTPUT_ROOT):
        shutil.rmtree(OUTPUT_ROOT)
    os.makedirs(OUTPUT_ROOT)

    con = sqlite3.connect(DB_PATH)
    try:
        for label, junctions_csv, hadas_format in INPUT_FILES:
            print(f"\n=== {label} ({os.path.basename(junctions_csv)}) ===")
            for use_longest_cds, restrict_pdf_to_comparable in FLAG_COMBINATIONS:
                run_one(con, label, junctions_csv, hadas_format, use_longest_cds, restrict_pdf_to_comparable)
    finally:
        con.close()

    print(f"\nAll reference outputs written under {OUTPUT_ROOT}")


if __name__ == '__main__':
    main()
