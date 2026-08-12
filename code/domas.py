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
import utils

_DEFAULT_NUM_WORKERS = os.cpu_count() or 5


def parse_args():
    parser = argparse.ArgumentParser(description="Analyze junctions and detect domain changes across alternative transcripts.")
    parser.add_argument("-format", required=False, choices=["internal", "internal2", "ioe", "leafcutter", "rmats", "majiq"],
                         help="Input file format: 'internal' for the internal comparative "
                              "splicing Excel file, 'internal2' for the internal supplementary "
                              "table Excel file (LeafCutter significant results: Cluster, "
                              "Splicing junction, deltaPSI, p.Adjust, Genes), 'ioe' for a SUPPA "
                              ".ioe file (or a directory of them, see -ioe_pattern), 'leafcutter' "
                              "for a pair of leafcutter_ds output files (see -lc_sig / "
                              "-lc_effect), 'rmats' for an rMATS-turbo output directory (the five "
                              "[Event].MATS.JC.txt files, passed via -input), 'majiq' for a "
                              "MAJIQ voila TSV file (passed via -input)")
    parser.add_argument("-input", required=False, default=None, type=str,
                         help="Path to the input for -format internal/internal2/ioe/rmats/majiq "
                              "(Excel for internal and internal2; a single .ioe file or a "
                              "directory of .ioe files for ioe; an rMATS-turbo output directory "
                              "for rmats; a voila TSV file for majiq). Not used with -format "
                              "leafcutter.")
    parser.add_argument("-lc_sig", required=False, default=None, type=str,
                         help="Path to the leafcutter_ds_cluster_significance.txt file "
                              "(only used with -format leafcutter)")
    parser.add_argument("-lc_effect", required=False, default=None, type=str,
                         help="Path to the leafcutter_ds_effect_sizes.txt file "
                              "(only used with -format leafcutter)")
    parser.add_argument("-species", required=False, choices=sorted(alternative_splicing.utils.SPECIE_DB_NAME),
                         help="Species the input was produced from: human, mouse, rat, "
                              "zebrafish or frog. Required for every format except internal, "
                              "which is a human/mouse comparison and carries the species per "
                              "row. rMATS, MAJIQ, SUPPA and internal2 files contain no species "
                              "field, so it cannot be read from them. DOMAS aborts if the gene "
                              "ids turn out to belong to a different species.")
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
    parser.add_argument("-pdf", action="store_true",
                         help="Generate a per-gene PDF (only honored with -format internal; "
                              "the ioe path never generates PDFs regardless of this flag). "
                              "Off by default: a full-scale internal run would otherwise "
                              "produce a PDF per gene.")
    parser.add_argument("-stats", action="store_true",
                         help="Generate the results_stats.py report (event distribution, "
                              "domain frequency, etc.) for the produced -output_csv after "
                              "the run. Off by default: the report is a separate set of "
                              "figures and tables, and takes longer than the analysis on a "
                              "large result.")
    parser.add_argument("-stats_out_dir", type=str, default=None,
                         help="Directory for the -stats report (default: same "
                              "directory as -output_csv)")
    parser.add_argument("-max_clusters", type=int, default=0,
                         help="If > 0, analyze only the first N clusters (sorted). Caps the "
                              "amount of work; used by the web GUI to process the first 100 "
                              "clusters. 0 (default) means no limit.")
    parser.add_argument("-omit_non_comparable", action="store_true",
                         help="Leave out rows for non-comparable transcripts (e.g. those with "
                              "a gene_not_in_db / feature_not_mapped / no_unique_features "
                              "event), keeping only transcripts that were actually compared "
                              "to the canonical one. By default they are written, so the "
                              "output CSV accounts for every cluster in the input - including "
                              "the ones no comparison was possible for.")
    parser.add_argument("-show_only_compared", action="store_true",
                         help="Draw only the canonical transcript and the transcripts whose "
                              "domains were actually compared to it (only honored with -pdf). "
                              "Leaves out the ones the analysis never compared: a transcript "
                              "carrying no feature of the event, none unique to it, one subsumed "
                              "by a larger event, and one its own group passed over in favour of "
                              "another (transcript_not_chosen). By default every transcript of "
                              "the gene is drawn.")
    parser.add_argument("-extra_columns", action="store_true",
                         help="Add three columns to the output CSV: `rank`, naming the "
                              "canonical transcript's exons that the event's junctions join "
                              "(E2_E4, or *_E5 where a junction's splice site is not an exon "
                              "edge of it), and canonical_junction_in_cds / alternative_junction_in_cds, "
                              "saying "
                              "whether those junctions fall inside the coding sequence of the "
                              "canonical and of the compared transcript (yes / partial / no, "
                              "or no_cds for a transcript with no annotated protein). Works "
                              "for every input format. Off by default.")
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
    # -species cannot describe it. Every other format must state one: three carry
    # no species field, and a GeneID-keyed gene names none in its id either.
    if args.format == "internal":
        if args.species:
            parser.error("-species does not apply to -format internal: that input is a "
                         "human/mouse comparison and carries the species per row.")
    elif not args.species:
        parser.error(f"-species is required for -format {args.format} "
                     f"(one of: {', '.join(sorted(alternative_splicing.utils.SPECIE_DB_NAME))})")
    if not args.dochap:
        parser.error("-dochap is required")
    if not os.path.exists(args.dochap):
        parser.error(f"DoChaP db not found at {args.dochap}")

    # -input (internal/internal2/ioe) and the leafcutter -lc_sig/-lc_effect pair
    # name different formats and cannot be combined.
    input_provided = args.input is not None
    leafcutter_provided = args.lc_sig is not None or args.lc_effect is not None
    if input_provided and leafcutter_provided:
        sys.stderr.write(
            "Error: mixed input formats - provide either -input (for -format internal/internal2/ioe/rmats/majiq) "
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
    args = parse_args()

    # -output_csv may name a directory that does not exist yet; the analysis
    # writes straight to that path, so create it before the run rather than
    # failing at the end with the results already computed.
    output_dir = os.path.dirname(os.path.abspath(args.output_csv))
    os.makedirs(output_dir, exist_ok=True)

    # One handler, a file beside the output CSV. Every module logs through
    # logging.getLogger(__name__), so this is the only place a destination is
    # chosen and the run leaves nothing on the console.
    log_file = os.path.join(output_dir, 'domas.log')
    # Emptied once here, then opened for APPENDING - not with filemode='w'. The
    # worker processes open the same file themselves (see _init_worker), and a
    # 'w' handle in the parent keeps its own offset from 0, so every parent
    # record overwrites bytes the workers have already appended: the run's first
    # per-cluster records come back torn in half, or vanish. Appending puts every
    # write at the current end of file, in either process.
    open(log_file, 'w').close()
    logging.basicConfig(
        filename=log_file,
        filemode='a',
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        force=True,
    )
    # Clusters are analysed in worker processes, which are spawned and so start
    # with no logging configuration of their own. Naming the file in the
    # environment - inherited by the children - keeps their records in the same
    # log instead of the console (see junction_analisys._init_worker).
    os.environ['DOMAS_LOG_FILE'] = log_file

    # Progress on the console: milestones, warnings and errors only. Everything
    # at INFO - per-chunk writes, per-cluster reasons - stays in the log file.
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(utils.PROGRESS)
    console.setFormatter(logging.Formatter("%(message)s"))
    logging.getLogger().addHandler(console)

    con = sqlite3.connect(args.dochap)
    try:
        if args.format == "leafcutter":
            alternative_splicing.analyze_leafcutter_input(
                con, significance_file=args.lc_sig, effect_sizes_file=args.lc_effect,
                output_csv=args.output_csv, specie=args.species, num_workers=args.num_workers, max_clusters=args.max_clusters,
                filter_non_comparable=args.omit_non_comparable,
                write_all_comparable=args.write_all_comparable,
                extra_columns=args.extra_columns,
            )
        elif args.format == "rmats":
            alternative_splicing.analyze_rmats_input(
                con, rmats_dir=args.input, output_csv=args.output_csv, specie=args.species, num_workers=args.num_workers, max_clusters=args.max_clusters,
                filter_non_comparable=args.omit_non_comparable,
                write_all_comparable=args.write_all_comparable,
                extra_columns=args.extra_columns,
            )
        elif args.format == "majiq":
            alternative_splicing.analyze_voila_input(
                con, voila_tsv=args.input, output_csv=args.output_csv, specie=args.species, num_workers=args.num_workers, max_clusters=args.max_clusters,
                filter_non_comparable=args.omit_non_comparable,
                write_all_comparable=args.write_all_comparable,
                extra_columns=args.extra_columns,
            )
        elif args.format == "internal2":
            alternative_splicing.analyze_internal2_input(
                con, input_file=args.input, output_csv=args.output_csv, specie=args.species,
                num_workers=args.num_workers, max_clusters=args.max_clusters,
                filter_non_comparable=args.omit_non_comparable,
                write_all_comparable=args.write_all_comparable,
                extra_columns=args.extra_columns,
            )
        elif args.format == "internal":
            print_genes = [g.strip() for g in args.gene_ids.split(',')] if args.gene_ids else None
            alternative_splicing.analyze_hadas_input(
                con, input_file=args.input, output_csv=args.output_csv,
                print_genes=print_genes, num_workers=args.num_workers, create_pdf=args.pdf,
                max_clusters=args.max_clusters, filter_non_comparable=args.omit_non_comparable,
                write_all_comparable=args.write_all_comparable,
                extra_columns=args.extra_columns,
                restrict_pdf_to_comparable=args.show_only_compared,
            )
        elif os.path.isdir(args.input):
            alternative_splicing.analyze_ioe_files(
                con, input_path=args.input, pattern=args.ioe_pattern, output_csv=args.output_csv, specie=args.species,
                examples_per_event=args.examples_per_event, num_workers=args.num_workers, max_clusters=args.max_clusters,
                filter_non_comparable=args.omit_non_comparable,
                write_all_comparable=args.write_all_comparable,
                extra_columns=args.extra_columns,
            )
        else:
            alternative_splicing.analyze_ioe_file(
                con, ioe_file=args.input, output_csv=args.output_csv, specie=args.species, num_workers=args.num_workers, max_clusters=args.max_clusters,
                filter_non_comparable=args.omit_non_comparable,
                write_all_comparable=args.write_all_comparable,
                extra_columns=args.extra_columns,
            )
    finally:
        con.close()

    if args.stats:
        import results_stats
        stats_out_dir = args.stats_out_dir or os.path.dirname(os.path.abspath(args.output_csv))
        if args.format == "internal":
            results_stats.generate_report(hadas_file=args.output_csv, out_dir=stats_out_dir)
        else:
            results_stats.generate_report(ioe_file=args.output_csv, out_dir=stats_out_dir)


if __name__ == "__main__":
    main()
