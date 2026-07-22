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
    parser.add_argument("-format", required=True, choices=["hadas", "ioe", "leafcutter", "rmats", "majiq"],
                         help="Input file format: 'hadas' for a hadas-style comparative "
                              "splicing Excel file, 'ioe' for a SUPPA .ioe file (or a "
                              "directory of them, see -ioe_pattern), 'leafcutter' for a pair "
                              "of leafcutter_ds output files (see -lc_sig / -lc_effect), "
                              "'rmats' for an rMATS-turbo output directory (the five "
                              "[Event].MATS.JC.txt files, passed via -input), 'majiq' for a "
                              "MAJIQ voila TSV file (passed via -input)")
    parser.add_argument("-input", required=False, default=None, type=str,
                         help="Path to the input for -format hadas/ioe/rmats/majiq (Excel for "
                              "hadas; a single .ioe file or a directory of .ioe files for ioe; "
                              "an rMATS-turbo output directory for rmats; a voila TSV file for "
                              "majiq). Not used with -format leafcutter.")
    parser.add_argument("-lc_sig", required=False, default=None, type=str,
                         help="Path to the leafcutter_ds_cluster_significance.txt file "
                              "(only used with -format leafcutter)")
    parser.add_argument("-lc_effect", required=False, default=None, type=str,
                         help="Path to the leafcutter_ds_effect_sizes.txt file "
                              "(only used with -format leafcutter)")
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
    parser.add_argument("-num_workers", type=int, default=_DEFAULT_NUM_WORKERS,
                         help=f"Number of parallel worker processes (default: this machine's "
                              f"CPU count, {_DEFAULT_NUM_WORKERS})")
    parser.add_argument("-use_representative_domains", action="store_true",
                         help="Pull domains from the RepresentativeDomains table where "
                              "available, falling back to DomainEvent/DomainType per "
                              "protein (default: DomainEvent/DomainType only)")
    parser.add_argument("-no_pdf", action="store_true",
                         help="Skip PDF generation (only honored with -format hadas; the "
                              "ioe path never generates PDFs regardless of this flag). "
                              "Useful for a full-scale hadas run where a PDF per gene "
                              "would otherwise always be produced.")
    parser.add_argument("-run_stats", action="store_true",
                         help="After the run, generate the results_stats.py report "
                              "(event distribution, domain frequency, etc.) for the "
                              "produced -output_csv.")
    parser.add_argument("-stats_out_dir", type=str, default=None,
                         help="Directory for the -run_stats report (default: same "
                              "directory as -output_csv)")
    parser.add_argument("-max_clusters", type=int, default=0,
                         help="If > 0, analyze only the first N clusters (sorted). Caps the "
                              "amount of work; used by the web GUI to process the first 100 "
                              "clusters. 0 (default) means no limit.")
    parser.add_argument("-filter_non_comparable", action="store_true",
                         help="Keep only transcripts that were actually compared to the "
                              "canonical transcript in the output CSV; drop rows for "
                              "non-comparable (not chosen) transcripts, e.g. those with a "
                              "gene_not_in_db / junction_not_mapped / no_unique_junctions event.")
    args = parser.parse_args()

    if not os.path.exists(args.dochap):
        parser.error(f"DoChaP db not found at {args.dochap}")

    # Reject mixing inputs from different formats: -input (hadas/ioe) must not be
    # combined with the leafcutter -lc_sig/-lc_effect pair. Exit with -1 per spec.
    input_provided = args.input is not None
    leafcutter_provided = args.lc_sig is not None or args.lc_effect is not None
    if input_provided and leafcutter_provided:
        sys.stderr.write(
            "Error: mixed input formats - provide either -input (for -format hadas/ioe/rmats/majiq) "
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
    # junction_analisys.py's JunctionsAnalysis already logs progress every
    # 10,000 clusters (with an ETA) via self.logger.info(...) - but nothing
    # in the actual call chain from this CLI ever configured a logging
    # handler (alternative_splicing.py only does that inside its own
    # __main__ guard, which doesn't run when it's imported as a module, as
    # it is here), so every one of those log records was silently dropped:
    # logging with no handler configured produces no output at all below
    # WARNING. This is what actually enables that existing progress logging
    # for a real domas.py run, not new logging of its own.
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    args = parse_args()
    con = sqlite3.connect(args.dochap)
    try:
        if args.format == "leafcutter":
            alternative_splicing.analyze_leafcutter_input(
                con, args.lc_sig, args.lc_effect, args.output_csv, num_workers=args.num_workers,
                use_representative_domains=args.use_representative_domains, max_clusters=args.max_clusters,
                filter_non_comparable=args.filter_non_comparable,
            )
        elif args.format == "rmats":
            alternative_splicing.analyze_rmats_input(
                con, args.input, args.output_csv, num_workers=args.num_workers,
                use_representative_domains=args.use_representative_domains, max_clusters=args.max_clusters,
                filter_non_comparable=args.filter_non_comparable,
            )
        elif args.format == "majiq":
            alternative_splicing.analyze_voila_input(
                con, args.input, args.output_csv, num_workers=args.num_workers,
                use_representative_domains=args.use_representative_domains, max_clusters=args.max_clusters,
                filter_non_comparable=args.filter_non_comparable,
            )
        elif args.format == "hadas":
            print_genes = [g.strip() for g in args.gene_ids.split(',')] if args.gene_ids else None
            alternative_splicing.analyze_hadas_input(
                con, args.input, args.output_csv, print_genes=print_genes, num_workers=args.num_workers,
                use_representative_domains=args.use_representative_domains, create_pdf=not args.no_pdf,
                max_clusters=args.max_clusters, filter_non_comparable=args.filter_non_comparable,
            )
        elif os.path.isdir(args.input):
            alternative_splicing.analyze_ioe_files(
                con, args.input, args.ioe_pattern, args.output_csv,
                examples_per_event=args.examples_per_event, num_workers=args.num_workers,
                use_representative_domains=args.use_representative_domains, max_clusters=args.max_clusters,
                filter_non_comparable=args.filter_non_comparable,
            )
        else:
            alternative_splicing.analyze_ioe_file(
                con, args.input, args.output_csv, num_workers=args.num_workers,
                use_representative_domains=args.use_representative_domains, max_clusters=args.max_clusters,
                filter_non_comparable=args.filter_non_comparable,
            )
    finally:
        con.close()

    if args.run_stats:
        import results_stats
        stats_out_dir = args.stats_out_dir or os.path.dirname(os.path.abspath(args.output_csv))
        if args.format == "hadas":
            results_stats.generate_report(hadas_file=args.output_csv, out_dir=stats_out_dir)
        else:
            results_stats.generate_report(ioe_file=args.output_csv, out_dir=stats_out_dir)


if __name__ == "__main__":
    main()
