"""
Frontend CLI for the DOMAS junction/domain analysis pipeline.

Parses and validates command-line arguments, then delegates to the actual
analysis logic in alternative_splicing.py (which wires the chosen input
format into JunctionsAnalysis.analyze_junctions() in junction_analisys.py).
This module owns no analysis logic itself and nothing else should import
from it.
"""
import argparse
import os
import sqlite3

import alternative_splicing


def parse_args():
    parser = argparse.ArgumentParser(description="Analyze junctions and detect domain changes across alternative transcripts.")
    parser.add_argument("-format", required=True, choices=["hadas", "ioe"],
                         help="Input file format: 'hadas' for a hadas-style comparative "
                              "splicing Excel file, 'ioe' for a SUPPA .ioe file (or a "
                              "directory of them, see -ioe_pattern)")
    parser.add_argument("-input", required=True, type=str,
                         help="Path to the input file (Excel for hadas; a single .ioe file "
                              "or a directory of .ioe files for ioe)")
    parser.add_argument("-dochap", required=True, type=str, help="Path to the DoChaP sqlite db")
    parser.add_argument("-output_csv", type=str, default="junctions_analysis.csv", help="Path to the output csv")
    parser.add_argument("-gene_ids", type=str, default=None,
                         help="Comma-separated list of gene symbols to generate PDFs for "
                              "(only honored with -format hadas; default is all genes)")
    parser.add_argument("-ioe_pattern", type=str, default=r"output_prefix_.*_strict.ioe",
                         help="Filename regex used to find .ioe files when -input is a "
                              "directory (only used with -format ioe)")
    parser.add_argument("-examples_per_event", type=int, default=5,
                         help="Per event type, keep only this many example clusters with "
                              "the fewest transcripts (0 keeps every cluster); only used "
                              "when -input is a directory with -format ioe")
    args = parser.parse_args()

    if not os.path.exists(args.dochap):
        parser.error(f"DoChaP db not found at {args.dochap}")
    if not os.path.exists(args.input):
        parser.error(f"Input path not found at {args.input}")

    return args


def main():
    args = parse_args()
    con = sqlite3.connect(args.dochap)
    try:
        if args.format == "hadas":
            print_genes = [g.strip() for g in args.gene_ids.split(',')] if args.gene_ids else None
            alternative_splicing.analyze_hadas_input(con, args.input, args.output_csv, print_genes=print_genes)
        elif os.path.isdir(args.input):
            alternative_splicing.analyze_ioe_files(
                con, args.input, args.ioe_pattern, args.output_csv,
                examples_per_event=args.examples_per_event,
            )
        else:
            alternative_splicing.analyze_ioe_file(con, args.input, args.output_csv)
    finally:
        con.close()


if __name__ == "__main__":
    main()
