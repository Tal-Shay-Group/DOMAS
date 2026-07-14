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

import alternative_splicing

_DEFAULT_NUM_WORKERS = os.cpu_count() or 5


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
    args = parser.parse_args()

    if not os.path.exists(args.dochap):
        parser.error(f"DoChaP db not found at {args.dochap}")
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
        if args.format == "hadas":
            print_genes = [g.strip() for g in args.gene_ids.split(',')] if args.gene_ids else None
            alternative_splicing.analyze_hadas_input(
                con, args.input, args.output_csv, print_genes=print_genes, num_workers=args.num_workers,
                use_representative_domains=args.use_representative_domains, create_pdf=not args.no_pdf,
            )
        elif os.path.isdir(args.input):
            alternative_splicing.analyze_ioe_files(
                con, args.input, args.ioe_pattern, args.output_csv,
                examples_per_event=args.examples_per_event, num_workers=args.num_workers,
                use_representative_domains=args.use_representative_domains,
            )
        else:
            alternative_splicing.analyze_ioe_file(
                con, args.input, args.output_csv, num_workers=args.num_workers,
                use_representative_domains=args.use_representative_domains,
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
