"""
Frontend CLI for the DOMAS junction/domain analysis pipeline.

Parses and validates command-line arguments, then delegates to the actual
analysis logic in alternative_splicing.py (which wires the chosen input
format into JunctionsAnalysis.analyze_junctions() in junction_analisys.py).
This module owns no analysis logic itself and nothing else should import
from it.
"""
import argparse
import logging
import os
import sqlite3
import sys

import alternative_splicing

_DEFAULT_NUM_WORKERS = os.cpu_count() or 5


def parse_args():
    parser = argparse.ArgumentParser(description="Analyze junctions and detect domain changes across alternative transcripts.")
    parser.add_argument("-format", required=False, choices=["internal", "ioe", "leafcutter", "rmats", "majiq"],
                         help="Input file format: 'internal' for the internal comparative "
                              "splicing Excel file, 'ioe' for a SUPPA .ioe file (or a "
                              "directory of them, see -ioe_pattern), 'leafcutter' for a pair "
                              "of leafcutter_ds output files (see -lc_sig / -lc_effect), "
                              "'rmats' for an rMATS-turbo output directory (the five "
                              "[Event].MATS.JC.txt files, passed via -input), 'majiq' for a "
                              "MAJIQ voila TSV file (passed via -input)")
    parser.add_argument("-input", required=False, default=None, type=str,
                         help="Path to the input for -format internal/ioe/rmats/majiq (Excel for "
                              "internal; a single .ioe file or a directory of .ioe files for ioe; "
                              "an rMATS-turbo output directory for rmats; a voila TSV file for "
                              "majiq). Not used with -format leafcutter.")
    parser.add_argument("-lc_sig", required=False, default=None, type=str,
                         help="Path to the leafcutter_ds_cluster_significance.txt file "
                              "(only used with -format leafcutter)")
    parser.add_argument("-lc_effect", required=False, default=None, type=str,
                         help="Path to the leafcutter_ds_effect_sizes.txt file "
                              "(only used with -format leafcutter)")
    parser.add_argument("-specie", required=False, choices=sorted(alternative_splicing.utils.SPECIE_DB_NAME),
                         help="Species the input was produced from: human, mouse, rat, "
                              "zebrafish or frog. Required for every format except internal, "
                              "which is a human/mouse comparison and carries the species per "
                              "row. rMATS, MAJIQ and SUPPA files contain no species field, so "
                              "it cannot be read from them. DOMAS aborts if the gene ids turn "
                              "out to belong to a different species.")
    parser.add_argument("-dochap", required=False, type=str, help="Path to the DoChaP sqlite db")
    parser.add_argument("-output_csv", type=str, default="junctions_analysis.csv", help="Path to the output csv")
    parser.add_argument("-gene_ids", type=str, default=None,
                         help="Comma-separated list of gene symbols to generate PDFs for "
                              "(only honored with -format internal together with -pdf; "
                              "default is all genes)")
    parser.add_argument("-ioe_pattern", type=str, default=r"output_prefix_.*_strict.ioe",
                         help="Filename regex used to find .ioe files when -input is a "
                              "directory (only used with -format ioe)")
    parser.add_argument("-examples_per_event", type=int, default=0,
                         help="Per event type, keep only this many example clusters with "
                              "the fewest transcripts; only used when -input is a directory "
                              "with -format ioe (SUPPA). 0 (default) keeps every cluster; a "
                              "subset is a debugging convenience, not something to take "
                              "silently from a SUPPA directory.")
    parser.add_argument("-num_workers", type=int, default=_DEFAULT_NUM_WORKERS,
                         help=f"Number of parallel worker processes (default: this machine's "
                              f"CPU count, {_DEFAULT_NUM_WORKERS})")
    parser.add_argument("-no_representative_domains", action="store_true",
                         help="Use DomainEvent/DomainType only. By default domains come "
                              "from the RepresentativeDomains table where available, "
                              "falling back to DomainEvent/DomainType per protein.")
    parser.add_argument("-pdf", action="store_true",
                         help="Generate a per-gene PDF (only honored with -format internal; "
                              "the ioe path never generates PDFs regardless of this flag). "
                              "Off by default: a full-scale internal run would otherwise "
                              "produce a PDF per gene.")
    parser.add_argument("-no_stats", action="store_true",
                         help="Skip the results_stats.py report. By default the report "
                              "(event distribution, domain frequency, etc.) is generated "
                              "for the produced -output_csv after the run.")
    parser.add_argument("-stats_out_dir", type=str, default=None,
                         help="Directory for the stats report (default: same "
                              "directory as -output_csv)")
    parser.add_argument("-max_clusters", type=int, default=0,
                         help="If > 0, analyze only the first N clusters (sorted). Caps the "
                              "amount of work; used by the web GUI to process the first 100 "
                              "clusters. 0 (default) means no limit.")
    parser.add_argument("-keep_non_comparable", action="store_true",
                         help="Keep rows for non-comparable transcripts (e.g. those with a "
                              "gene_not_in_db / feature_not_mapped / no_unique_features "
                              "event). By default the output CSV holds only transcripts "
                              "that were actually compared to the canonical transcript.")
    parser.add_argument("-write_all_comparable", action="store_true",
                         help="Compare every comparable transcript to the canonical one and "
                              "keep a row for each, adding the is_most_like_canonical and "
                              "is_longest_cds columns that say which one the selection rule "
                              "picks. By default only that transcript is compared - the one "
                              "most like the canonical, or the longest-CDS one where none "
                              "qualifies - and the two columns are omitted.")

    args = parser.parse_args()

    if not args.format:
        parser.error("-format is required")

    # internal states a species per row (it is a human/mouse comparison), so one
    # -specie cannot describe it. Every other format must state one: three carry
    # no species field, and a GeneID-keyed gene names none in its id either.
    if args.format == "internal":
        if args.specie:
            parser.error("-specie does not apply to -format internal: that input is a "
                         "human/mouse comparison and carries the species per row.")
    elif not args.specie:
        parser.error(f"-specie is required for -format {args.format} "
                     f"(one of: {', '.join(sorted(alternative_splicing.utils.SPECIE_DB_NAME))})")
    if not args.dochap:
        parser.error("-dochap is required")
    if not os.path.exists(args.dochap):
        parser.error(f"DoChaP db not found at {args.dochap}")

    # -input (internal/ioe) and the leafcutter -lc_sig/-lc_effect pair name
    # different formats and cannot be combined.
    input_provided = args.input is not None
    leafcutter_provided = args.lc_sig is not None or args.lc_effect is not None
    if input_provided and leafcutter_provided:
        sys.stderr.write(
            "Error: mixed input formats - provide either -input (for -format internal/ioe/rmats/majiq) "
            "OR -lc_sig/-lc_effect (for -format leafcutter), not both.\n")
        sys.exit(-1)

    if args.format == "leafcutter":
        if input_provided:
            sys.stderr.write(
                "Error: mixed input formats - -format leafcutter uses -lc_sig/-lc_effect, "
                "not -input.\n")
            sys.exit(-1)
        if args.lc_sig is None or args.lc_effect is None:
            parser.error("-format leafcutter requires both -lc_sig and -lc_effect")
        if not os.path.exists(args.lc_sig):
            parser.error(f"leafcutter significance file not found at {args.lc_sig}")
        if not os.path.exists(args.lc_effect):
            parser.error(f"leafcutter effect-sizes file not found at {args.lc_effect}")
    else:
        if leafcutter_provided:
            sys.stderr.write(
                f"Error: mixed input formats - -lc_sig/-lc_effect are only for "
                f"-format leafcutter, not -format {args.format}.\n")
            sys.exit(-1)
        if not input_provided:
            parser.error(f"-format {args.format} requires -input")
        if not os.path.exists(args.input):
            parser.error(f"Input path not found at {args.input}")

    return args


def main():
    # JunctionsAnalysis logs progress every 10,000 clusters at INFO. Nothing
    # else in this call chain configures a handler, and an unconfigured logger
    # drops everything below WARNING, so those records need this to appear.
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    args = parse_args()

    con = sqlite3.connect(args.dochap)
    try:
        if args.format == "leafcutter":
            alternative_splicing.analyze_leafcutter_input(
                con, significance_file=args.lc_sig, effect_sizes_file=args.lc_effect,
                output_csv=args.output_csv, specie=args.specie, num_workers=args.num_workers,
                use_representative_domains=not args.no_representative_domains, max_clusters=args.max_clusters,
                filter_non_comparable=not args.keep_non_comparable,
                write_all_comparable=args.write_all_comparable,
            )
        elif args.format == "rmats":
            alternative_splicing.analyze_rmats_input(
                con, rmats_dir=args.input, output_csv=args.output_csv, specie=args.specie, num_workers=args.num_workers,
                use_representative_domains=not args.no_representative_domains, max_clusters=args.max_clusters,
                filter_non_comparable=not args.keep_non_comparable,
                write_all_comparable=args.write_all_comparable,
            )
        elif args.format == "majiq":
            alternative_splicing.analyze_voila_input(
                con, voila_tsv=args.input, output_csv=args.output_csv, specie=args.specie, num_workers=args.num_workers,
                use_representative_domains=not args.no_representative_domains, max_clusters=args.max_clusters,
                filter_non_comparable=not args.keep_non_comparable,
                write_all_comparable=args.write_all_comparable,
            )
        elif args.format == "internal":
            print_genes = [g.strip() for g in args.gene_ids.split(',')] if args.gene_ids else None
            alternative_splicing.analyze_hadas_input(
                con, input_file=args.input, output_csv=args.output_csv,
                print_genes=print_genes, num_workers=args.num_workers,
                use_representative_domains=not args.no_representative_domains, create_pdf=args.pdf,
                max_clusters=args.max_clusters, filter_non_comparable=not args.keep_non_comparable,
                write_all_comparable=args.write_all_comparable,
            )
        elif os.path.isdir(args.input):
            alternative_splicing.analyze_ioe_files(
                con, input_path=args.input, pattern=args.ioe_pattern, output_csv=args.output_csv, specie=args.specie,
                examples_per_event=args.examples_per_event, num_workers=args.num_workers,
                use_representative_domains=not args.no_representative_domains, max_clusters=args.max_clusters,
                filter_non_comparable=not args.keep_non_comparable,
                write_all_comparable=args.write_all_comparable,
            )
        else:
            alternative_splicing.analyze_ioe_file(
                con, ioe_file=args.input, output_csv=args.output_csv, specie=args.specie, num_workers=args.num_workers,
                use_representative_domains=not args.no_representative_domains, max_clusters=args.max_clusters,
                filter_non_comparable=not args.keep_non_comparable,
                write_all_comparable=args.write_all_comparable,
            )
    finally:
        con.close()

    if not args.no_stats:
        import results_stats
        stats_out_dir = args.stats_out_dir or os.path.dirname(os.path.abspath(args.output_csv))
        if args.format == "internal":
            results_stats.generate_report(hadas_file=args.output_csv, out_dir=stats_out_dir)
        else:
            results_stats.generate_report(ioe_file=args.output_csv, out_dir=stats_out_dir)


if __name__ == "__main__":
    main()
