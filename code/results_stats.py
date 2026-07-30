"""
Junction analysis — statistics for a single results.csv, plus cross-file
comparison for two or more of them (e.g. IOE genome-derived vs Hadas real
reads).

Per-file analyses (analyze_file()):
  1. Event-type distribution          — unanalyzable reasons + analyzable outcomes,
                                         one result per cluster (see _cluster_event_labels);
                                         "mixed" when a cluster reduces to more than one
                                         reason, "domain swap" for the common
                                         added_domain+dropped domain case specifically
  2. Mixed combinations               — which exact reason-sets make up "mixed", and how often
  3. Most frequent domains            — overall and per event type
  4. Domain frequency scatter         — per-domain rank and raw occurrence count,
                                         one species vs the other; both across every
                                         shared domain and restricted to #3's top-N
                                         (skipped with a warning if there isn't
                                         exactly one 'specie' pair present)
  5. Length change — shorter          — absolute Δaa and relative Δ%
  6. Length change — longer           — absolute Δaa and relative Δ%
  7. Species comparison               — human vs mouse, overall and per gene
                                         (skipped with a warning if there's no
                                         'specie' column)
  8. AS splice type vs outcome        — cluster prefix encodes A3/A5/SE/RI/…
                                         (skipped with a warning if no cluster
                                         carries an AS-type prefix)
  9. Domain count change              — copies gained/lost vs canonical
 10. Domain severity spectrum         — continuous % of canonical domain retained

Cross-file analysis (compare_files()):
 11. Event-type distribution comparison — side by side across files
 12. Severity spectrum comparison        — overlaid histograms across files

A column that's expected but absent from a given file (e.g. IOE files have no
'specie' column) triggers a warning at load time, and analyses that depend on
it are skipped rather than failing.

Outputs: with generate_report(pdf_report=True) (see _run_pdf_report), every
chart and table above is rendered straight into one PDF per file
(<label>_report.pdf) - no intermediate PNG/CSV files, except a table too
large for the PDF's row cap (e.g. species_per_gene_comparison, which can
have tens of thousands of rows) still gets a full CSV alongside the PDF
(see _save_table). Without pdf_report=True (the default, and always the
case when calling analyze_file() directly rather than through
generate_report()), each chart/table is instead written as its own
standalone file in OUT_DIR, prefixed with its label:
  event_distribution_<label>.png
  mixed_combinations_<label>.csv
  domain_frequency_<label>.png
  domain_rank_scatter_<label>.png / domain_rank_scatter_top20_<label>.png
  domain_count_scatter_<label>.png / domain_count_scatter_top20_<label>.png
  domain_descriptions_<label>.csv (only with fetch_domain_descriptions=True)
  length_change_shorter_<label>.png
  length_change_longer_<label>.png
  species_overall_<label>.png
  species_per_gene_<label>.png
  species_per_gene_comparison_<label>.csv
  splice_type_vs_outcome_<label>.png
  domain_count_change_<label>.png
  severity_spectrum_<label>.png
  event_distribution_comparison.png
  severity_spectrum_comparison.png
"""

import argparse
import io
import math
import os
import re
import sqlite3
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────────
CODE_DIR     = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(CODE_DIR)
OUT_DIR      = PROJECT_ROOT

IOE_FILE   = os.path.join(PROJECT_ROOT, "ioe_full_results.csv")
HADAS_FILE = os.path.join(PROJECT_ROOT, "hadas_results.csv")

# DoChaP merged DB - source of the RepresentativeDomains.description text used
# for the top-domains description table (matches alternative_splicing.py).
DEFAULT_DOCHAP_DB_PATH = "/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-web/DB_merged.sqlite"

# ── Event-type taxonomy ────────────────────────────────────────────────────────
# Mirrors the full set of event_type values ClusterAnalysisResult.add_event()
# can actually produce (see junction_analisys.py: the literal add_event(...)
# calls for the unanalyzable/skip reasons, classify_domain_change() for the
# analyzed outcomes).
# "gene_not_in_db" = gene's ENSG ID not found in the reference database (e.g., lncRNAs)
# "no_canonical_transcript" = gene found in DB but no canonical transcript
UNANALYZABLE_TYPES = [
    "gene_not_in_db",
    "no_canonical_transcript",
    "only_one_transcript",
    "junction_not_mapped",
    "no_unique_junctions",
    "transcript_doesnt_have_junctions",
    "no_canonical_junctions",
    "no_domains_in_region",
]

ANALYZED_TYPES = [
    "dropped domain",
    "added_domain",
    "shorter",
    "longer",
    "same",
    "split domain",
    "merged domain",
    "reduced_domain_number",
    "increased_domain_number",
]

ALL_EVENT_TYPES = UNANALYZABLE_TYPES + ANALYZED_TYPES

# Columns that aren't always present in a results.csv (depends on input
# format / how it was produced) - analyses that need one of these are
# skipped with a warning instead of failing when it's absent.
OPTIONAL_COLUMNS = ["specie", "cluster"]

# ── Visual style ───────────────────────────────────────────────────────────────
COLORS = {
    "human": "#2E75B6",
    "mouse": "#ED7D31",
    "IOE":   "#5B9BD5",
    "Hadas": "#ED7D31",
}

LENGTH_CHANGE_COLORS = {"shorter": "#FF8C00", "longer": "#FFC000"}


def _warn(message):
    print(f"  WARNING: {message}")


# unanalyzable: muted, distinguishable hues (not a grey ramp - those were hard
# to tell apart). analyzed: distinct, vivid colours. Each list's length is
# asserted against its UNANALYZABLE_TYPES/ANALYZED_TYPES below - a silent
# length mismatch (e.g. forgetting to add a color when a new event type is
# added) would otherwise make zip() truncate and shift every later color by
# one, which happened once already when "gene_not_in_db" was added.
_UNANALYZABLE_COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2", "#937860", "#CCB974", "#64B5CD"]
_ANALYZED_COLORS = ["#C00000", "#2E8B57", "#FF8C00", "#FFC000", "#70AD47", "#7030A0", "#4472C4", "#A6761D", "#E7298A"]
assert len(_UNANALYZABLE_COLORS) == len(UNANALYZABLE_TYPES), \
    f"event_color(): {len(_UNANALYZABLE_COLORS)} unanalyzable colors for {len(UNANALYZABLE_TYPES)} UNANALYZABLE_TYPES"
assert len(_ANALYZED_COLORS) == len(ANALYZED_TYPES), \
    f"event_color(): {len(_ANALYZED_COLORS)} analyzed colors for {len(ANALYZED_TYPES)} ANALYZED_TYPES"


def event_color(event_type):
    palette = dict(zip(ALL_EVENT_TYPES, _UNANALYZABLE_COLORS + _ANALYZED_COLORS))
    palette["domain swap"] = "#8B4513"  # distinct from both added_domain (sea green) and dropped domain (dark red)
    palette["mixed"] = "#555555"        # only used outside the pie (e.g. species_comparison bars) - the pie explodes "mixed" into MIXED_COMBO_SHADES instead
    return palette.get(event_type, "#888888")


# Shade ramp for the pie's exploded "mixed" sub-wedges (see event_distribution
# and _pdf_event_distribution_combined) - a distinct blue-grey hue not used by
# any single-reason wedge, so the whole family still reads as "part of mixed"
# while each shade stays distinguishable. Order means different things in the
# two callers: event_distribution() assigns by rank (darkest = most frequent
# combination in that one pie); _pdf_event_distribution_combined assigns by
# combo *identity* instead, so the same combination gets the same shade in
# every label's panel - see that function for why rank-based shading isn't
# safe to share across panels. Either way, the last shade is reserved for
# the "other mixed" catch-all.
MIXED_COMBO_SHADES = ["#37474F", "#546E7A", "#78909C", "#90A4AE", "#B0BEC5", "#CFD8DC"]


def event_sort_key(event_type):
    try:
        return ALL_EVENT_TYPES.index(event_type)
    except ValueError:
        pass
    if event_type == "domain swap":
        return len(ALL_EVENT_TYPES)
    if event_type == "mixed":
        return len(ALL_EVENT_TYPES) + 1
    return len(ALL_EVENT_TYPES) + 2


# Shorter chart-only labels. The analyzed side always renders under a title
# that already says "domain" (Domain Frequency, Domain Copy-Count Change,
# "Analyzable clusters" next to a domain-change legend, ...), so repeating
# "domain" in every category name just adds width without adding meaning.
# "transcript_doesnt_have_junctions" is shortened on its own merits - it's
# nearly twice as long as any other unanalyzable reason.
# Never used for filtering/CSV output - ANALYZED_TYPES, UNANALYZABLE_TYPES,
# ALL_EVENT_TYPES and the "event_type"/"label" column values stay unchanged.
SHORT_LABELS = {
    "gene_not_in_db": "not in DB",
    "dropped domain": "dropped",
    "added_domain": "added",
    "split domain": "split",
    "merged domain": "merged",
    # Not "fewer/more copies" - that phrasing is already the stacked-bar
    # legend text in domain_count_change() for a different, row-level
    # concept (this is the row's classified outcome type).
    "reduced_domain_number": "fewer domains",
    "increased_domain_number": "more domains",
    "domain swap": "swap",
    "transcript_doesnt_have_junctions": "lacks junctions",
}


def display_label(event_type):
    return SHORT_LABELS.get(event_type, event_type.replace("_", " "))


# Set (via _run_pdf_report) while a PDF report is being built, so save()/
# _save_table() route straight into it instead of writing scattered
# PNG/CSV files - see both functions below.
_ACTIVE_PDF = None

# Set (via _run_pdf_report) to a list while a PDF report is being built, so
# _save_table() collects tables instead of rendering them immediately -
# _run_pdf_report then renders them all together at the very end, after
# every section's charts, instead of each interleaved into its own section.
_DEFERRED_TABLES = None


def save(fig, filename):
    """
    Add `fig` as one or more pages to the currently-open PDF report (see
    _run_pdf_report), or - outside any PDF context, e.g. an ad-hoc call to
    a single analysis function with no report being built - fall back to
    writing a standalone PNG to OUT_DIR.
    """
    if _ACTIVE_PDF is not None:
        _savefig_to_pdf(_ACTIVE_PDF, fig, filename)
        return
    path = os.path.join(OUT_DIR, filename)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {path}")


def _save_table(df, filename, title=None, max_rows=300):
    """
    Add `df` as one or more table pages to the currently-open PDF report
    (see save()), or - outside any PDF context - write it as a standalone
    CSV to OUT_DIR, same fallback as save(). Within a PDF context, a table
    with more than max_rows still also gets a full CSV on disk (e.g.
    species_per_gene_comparison, which can have tens of thousands of rows) -
    the PDF only ever shows the first max_rows, so anything bigger needs
    somewhere to keep the rest.

    `title` is the human-readable page heading (should name both the table
    and its scope, e.g. "Mixed Combinations — Hadas Human") - falls back to
    the filename if not given. When _DEFERRED_TABLES is active (see
    _run_pdf_report), the table is collected rather than rendered
    immediately, so every table in the report ends up together at the end
    instead of interleaved into each section.
    """
    title = title or filename
    if _ACTIVE_PDF is not None:
        if _DEFERRED_TABLES is not None:
            _DEFERRED_TABLES.append((df, title))
        else:
            _pdf_csv_pages_from_df(_ACTIVE_PDF, df, title, max_rows=max_rows)
        if len(df) > max_rows:
            path = os.path.join(OUT_DIR, filename)
            df.to_csv(path, index=False)
            print(f"  → {path}  (full {len(df)} row(s) - the PDF only shows the first {max_rows})")
        return
    path = os.path.join(OUT_DIR, filename)
    df.to_csv(path, index=False)
    print(f"  → {path}")


# ── Data loading ───────────────────────────────────────────────────────────────

RESULTS_CSV_DTYPES = {
    "gene_symbol": "category",
    "specie": "category",
    "event_type": "category",
    "domain_name": "category",
    "cluster": "category",
    "transcript_id": "category",
    "c_domain_length": "float32",
    "t_domain_length": "float32",
    "c_domains_number": "float32",
    "t_domains_number": "float32",
}
# is_longest_cds/is_most_like_canonical are left to pandas' default (object,
# holding True/False/NaN) rather than forced to a dtype here - they're plain
# numpy-backed object arrays either way, so the multi-chunk concatenation
# path below (np.concatenate on .to_numpy()) handles them with no special
# casing, unlike pandas' nullable "boolean" extension dtype would need.

# Columns loaded but never referenced by any analysis below (verified by
# grep across this file) — safe to skip entirely to save memory on huge
# result files.
_UNUSED_COLS = {
    "canonical_transcript_id",
}


def _cluster_prefix(c):
    """Reduce a cluster id to its AS-type prefix (e.g. 'A3_ENSG...' -> 'A3')."""
    if not c or c.startswith("chr"):
        return ""
    return c.split("_")[0]


def _read_results_csv(path, chunk_rows=1_000_000):
    """Memory-efficient reader for a results.csv of any size.

    Reads in chunks with compact dtypes (category/float32/nullable-boolean)
    and drops columns that are never used downstream. The full "cluster" id
    is kept (as a category, not collapsed to its AS-type prefix) since
    select_representative_transcript() needs real per-cluster identity, not
    just its AS-type category - at full IOE scale (21M rows, ~1.1M distinct
    clusters) this costs well under 1GB, not the "materialise tens of
    millions of mostly-unique strings" blowup a plain object/string column
    would. For a small file this just runs as a single chunk - no
    meaningful overhead either way.
    """
    header_cols = list(pd.read_csv(path, nrows=0).columns)
    usecols = [c for c in header_cols if c not in _UNUSED_COLS]
    dtype = {k: v for k, v in RESULTS_CSV_DTYPES.items() if k in usecols}

    chunks = list(pd.read_csv(path, usecols=usecols, dtype=dtype, chunksize=chunk_rows))
    if len(chunks) == 1:
        return chunks[0]

    # NOTE: plain pd.concat(chunks) silently upcasts categorical columns back
    # to plain object dtype whenever the per-chunk categories differ (which
    # they always will, since each chunk only sees its own slice of values).
    # For a many-million-row file that turns our compact category/float32
    # columns back into huge Python-object arrays and is what was causing
    # this loader to run out of memory. Combine each column explicitly
    # instead so categoricals stay categoricals.
    cat_cols = [c for c in chunks[0].columns if isinstance(chunks[0][c].dtype, pd.CategoricalDtype)]
    other_cols = [c for c in chunks[0].columns if c not in cat_cols]

    data = {}
    for col in cat_cols:
        data[col] = pd.api.types.union_categoricals(
            [chunk[col] for chunk in chunks], sort_categories=False
        )
    for col in other_cols:
        data[col] = np.concatenate([chunk[col].to_numpy() for chunk in chunks])

    return pd.DataFrame(data, columns=list(chunks[0].columns))


def normalize_event_types(df):
    """
    classify_domain_change() suffixes its length-comparison outcome with
    "_domains" when comparing equal-count multi-domain groups (e.g.
    "longer_domains") instead of a single domain pair ("longer") - same
    underlying outcome, so normalise both forms to one label.

    Idempotent (safe to call on an already-normalised df), and independent
    of _read_results_csv()'s cluster-id collapsing - call this on any
    results.csv DataFrame, however it was loaded, before filtering by
    ANALYZED_TYPES/UNANALYZABLE_TYPES. Without it, "_domains"-suffixed rows
    silently fall outside both lists and get miscounted as neither analyzed
    nor unanalyzable.
    """
    df["event_type"] = df["event_type"].replace({
        "same_domains": "same",
        "longer_domains": "longer",
        "shorter_domains": "shorter",
    })
    return df


def _load_single(path, label):
    """Load and normalise one results.csv. Warns (and fills a blank column)
    for any of OPTIONAL_COLUMNS that isn't present, rather than failing."""
    df = _read_results_csv(path)
    df = normalize_event_types(df)

    # Every event_type the algorithm can produce must be in ALL_EVENT_TYPES
    # (see its definition above) - an unknown one means the taxonomy here has
    # drifted out of sync with junction_analisys.py and every count/percentage
    # below would silently be wrong (as happened before added_domain etc.
    # were added), so fail loudly instead of quietly under-reporting.
    unknown = set(df["event_type"].unique()) - set(ALL_EVENT_TYPES)
    if unknown:
        counts = df["event_type"].value_counts()
        details = ", ".join(f"{et!r} ({counts[et]} rows)" for et in sorted(unknown))
        raise ValueError(
            f"Unknown event_type value(s) in {label} ({path}): {details}. "
            f"Add them to UNANALYZABLE_TYPES or ANALYZED_TYPES in results_stats.py."
        )

    for col in OPTIONAL_COLUMNS:
        if col not in df.columns:
            _warn(f"'{col}' column not found in {label} ({path}) — analyses that need it will be skipped.")
            df[col] = ""
        elif df[col].isna().any():
            # A blank CSV field parses as NaN, not "" - and groupby() drops
            # NaN keys by default, which would silently vanish every row for
            # a format with no specie data (e.g. IOE) the moment anything
            # groups by specie together with another column.
            if isinstance(df[col].dtype, pd.CategoricalDtype) and "" not in df[col].cat.categories:
                df[col] = df[col].cat.add_categories([""])
            df[col] = df[col].fillna("")

    # Extract AS splice type from cluster prefix (e.g. "A3_ENSG…" → "A3").
    # Rows whose cluster id doesn't carry one (e.g. Hadas' "chr11:clu_6611",
    # or a missing cluster column) get "" and are excluded downstream.
    df["as_type"] = df["cluster"].apply(_cluster_prefix)

    return df


def select_representative_transcript(df, on_ambiguous="raise"):
    """
    Pick a single transcript per (specie, cluster) to keep for downstream
    analysis, so a cluster with several comparable transcripts isn't counted
    once per transcript.

    Only rows with an analyzed event_type are considered (see
    ANALYZED_TYPES) - unanalyzable rows don't represent an actual transcript
    comparison and are dropped before selection. Priority per cluster:
      1. the transcript tagged is_most_like_canonical
      2. if none is tagged, the transcript tagged is_longest_cds
      3. if neither tag is set, the cluster's sole remaining transcript - if
         more than one distinct transcript remains at this point, the
         cluster is ambiguous and is handled per `on_ambiguous`:
           "raise" (default): raise ValueError, naming a few examples.
           "drop": exclude the cluster's rows from the result entirely, and
             print a warning with the count dropped.

    analyze_junction() always tags exactly one transcript is_longest_cds for
    any cluster with >=1 comparable transcript (is_most_like_canonical goes
    unset when no transcript qualifies that rule, which is why step 2 exists),
    so the ambiguous case isn't just bad/missing input: it happens for real
    when the tagged transcript's own
    comparison landed on an unanalyzable outcome (typically
    no_domains_in_region) while other, untagged transcripts produced
    analyzed outcomes - the tag "wins" structurally but leaves nothing in
    the analyzed subset, and there's no principled way to prefer one of the
    remaining untagged transcripts over another.

    `df` must carry the real, full cluster id in its 'cluster' column - not
    the AS-type-prefix-collapsed one older versions of _read_results_csv()
    produced, which would group every cluster of the same AS type together
    instead of keeping each cluster distinct. normalize_event_types() is
    applied here regardless of how `df` was loaded, so "_domains"-suffixed
    event types are handled either way.

    Returns the subset of `df` belonging to each cluster's chosen transcript
    (every original column preserved).
    """
    if on_ambiguous not in ("raise", "drop"):
        raise ValueError(f"on_ambiguous must be 'raise' or 'drop', got {on_ambiguous!r}")

    df = normalize_event_types(df.copy())
    analyzed_df = df[df["event_type"].isin(ANALYZED_TYPES)]
    if analyzed_df.empty:
        return analyzed_df

    group_cols = ["specie", "cluster"] if "specie" in analyzed_df.columns else ["cluster"]

    most_like_pick = (
        analyzed_df[analyzed_df["is_most_like_canonical"] == True]
        .groupby(group_cols, observed=True, dropna=False)["transcript_id"].first()
    )
    longest_pick = (
        analyzed_df[analyzed_df["is_longest_cds"] == True]
        .groupby(group_cols, observed=True, dropna=False)["transcript_id"].first()
    )
    nunique = analyzed_df.groupby(group_cols, observed=True, dropna=False)["transcript_id"].nunique()

    chosen = most_like_pick.combine_first(longest_pick)

    unresolved = nunique.index.difference(chosen.index)
    if len(unresolved) > 0:
        ambiguous = nunique.loc[unresolved]
        ambiguous = ambiguous[ambiguous != 1]
        if len(ambiguous) > 0:
            if on_ambiguous == "raise":
                examples = ambiguous.head(5).to_dict()
                raise ValueError(
                    f"Cannot uniquely select a transcript for {len(ambiguous)} cluster(s): no row "
                    f"tagged is_most_like_canonical/is_longest_cds, and more than one distinct "
                    f"transcript remains. Examples ({{cluster: n_transcripts}}): {examples}"
                )
            _warn(f"select_representative_transcript: dropping {len(ambiguous)} ambiguous cluster(s) "
                  f"(no is_most_like_canonical/is_longest_cds tag, more than one transcript remains).")
        sole_pick = analyzed_df.groupby(group_cols, observed=True, dropna=False)["transcript_id"].first()
        resolvable = unresolved.difference(ambiguous.index)
        chosen = chosen.combine_first(sole_pick.loc[resolvable])

    chosen_df = chosen.rename("_chosen_transcript_id").reset_index()
    merged = analyzed_df.merge(chosen_df, on=group_cols, how="left")
    return merged[merged["transcript_id"] == merged["_chosen_transcript_id"]].drop(columns="_chosen_transcript_id")


# ══════════════════════════════════════════════════════════════════════════════
# 1. EVENT-TYPE DISTRIBUTION
# ══════════════════════════════════════════════════════════════════════════════

DOMAIN_SWAP_COMBO = frozenset({"added_domain", "dropped domain"})


def _collapse_to_cluster_label(sub_df, group_cols):
    """
    One label per group in sub_df: its sole distinct event_type; "domain
    swap" if the group's rows are exactly {added_domain, dropped domain} -
    common enough, and semantically distinct enough (one domain gained, one
    lost, in the same comparison) to name explicitly instead of burying in
    the generic bucket (see mixed_combinations() for what else ends up
    there); or "mixed" for any other multi-value case. Returns a Series
    indexed by group_cols.
    """
    grouped = sub_df.groupby(group_cols, observed=True, dropna=False)["event_type"]
    event_sets = grouped.apply(lambda s: frozenset(s))
    # .astype(object) first - event_type is a category dtype, which can't
    # hold "mixed"/"domain swap" via .where() unless already a category.
    sole = grouped.first().astype(object)

    label = sole.where(event_sets.apply(len) == 1, "mixed")
    label = label.where(event_sets != DOMAIN_SWAP_COMBO, "domain swap")
    return label


def _cluster_event_labels(df):
    """
    One row per (specie, cluster): whether it's "unanalyzable" or
    "analyzable", and its label - the sole event_type behind that, or
    "mixed" if more than one applies. Multiple transcript-comparison rows
    per cluster (e.g. no_unique_junctions can fire once per non-matching
    candidate transcript) would otherwise inflate a per-row count in a way
    that mixes "one row per cluster" event types (e.g.
    no_canonical_transcript) with "one row per candidate transcript" ones in
    the same tally.

    Unanalyzable clusters (no analyzed-event row at all): label from all of
    the cluster's (necessarily unanalyzable) rows.

    Analyzable clusters: label from just the representative transcript's
    rows (see select_representative_transcript) - a single transcript can
    still match several domain groups, each independently classified, hence
    still needs its own "mixed" case. A cluster select_representative_transcript
    can't resolve (ambiguous - see its docstring) is also labeled "mixed",
    since by definition more than one candidate transcript with a
    (potentially different) analyzed outcome remains.

    Returns a DataFrame with columns group_cols + ['cluster_status', 'label',
    'gene_symbol'] ('gene_symbol' included whenever `df` has one - constant
    within a cluster, so callers needing a per-gene breakdown - e.g.
    species_comparison() - don't need to re-derive it).
    """
    group_cols = ["specie", "cluster"] if "specie" in df.columns else ["cluster"]
    is_analyzed = df["event_type"].isin(ANALYZED_TYPES)

    has_analyzed = (
        df.assign(_is_analyzed=is_analyzed)
        .groupby(group_cols, observed=True, dropna=False)["_is_analyzed"].any()
    )
    analyzable_keys = has_analyzed[has_analyzed].index
    unanalyzable_keys = has_analyzed[~has_analyzed].index

    un_labels = _collapse_to_cluster_label(df[~is_analyzed], group_cols).reindex(unanalyzable_keys)

    representative_df = select_representative_transcript(df, on_ambiguous="drop")
    an_labels = _collapse_to_cluster_label(representative_df, group_cols).reindex(analyzable_keys, fill_value="mixed")

    result = pd.concat([
        un_labels.rename("label").to_frame().assign(cluster_status="unanalyzable"),
        an_labels.rename("label").to_frame().assign(cluster_status="analyzable"),
    ]).reset_index()

    if "gene_symbol" in df.columns:
        gene_lookup = df.groupby(group_cols, observed=True, dropna=False)["gene_symbol"].first()
        result = result.merge(gene_lookup.reset_index(), on=group_cols, how="left")

    return result


def event_distribution(df, label):
    """
    Plot and print the split between unanalyzable and analyzable clusters -
    one result per (specie, cluster), not per row (see
    _cluster_event_labels) - with a percentage breakdown of specific reasons
    within each group. The printed breakdown still reports a single "mixed"
    bucket for clusters that don't reduce to one reason, but the pie explodes
    that bucket into its actual top combinations (see _mixed_combo_sets),
    since "mixed" alone can be the largest slice and hides what's going on.
    """
    cluster_labels = _cluster_event_labels(df)
    unanalyzable_df = cluster_labels[cluster_labels["cluster_status"] == "unanalyzable"]
    analyzable_df = cluster_labels[cluster_labels["cluster_status"] == "analyzable"]

    n_total    = len(cluster_labels)
    n_unanalyz = len(unanalyzable_df)
    n_analyzed = len(analyzable_df)

    print(f"\n{'='*70}")
    print(f"EVENT DISTRIBUTION — {label}  (n={n_total} clusters)")
    print(f"{'='*70}")
    print(f"  Unanalyzable : {n_unanalyz}  ({100*n_unanalyz/n_total:.1f}%)")
    for et in UNANALYZABLE_TYPES + ["mixed"]:
        count = (unanalyzable_df["label"] == et).sum()
        if count:
            print(f"    {et:<42s} {count:3d}  ({100*count/n_unanalyz:.1f}% of unanalyzable)")
    print(f"  Analyzable   : {n_analyzed}  ({100*n_analyzed/n_total:.1f}%)")
    for et in ANALYZED_TYPES + ["domain swap", "mixed"]:
        count = (analyzable_df["label"] == et).sum()
        if count:
            print(f"    {et:<42s} {count:3d}  ({100*count/n_analyzed:.1f}% of analyzable)")

    # "mixed" can be the biggest slice of either pie (e.g. it's ~90% of
    # unanalyzable clusters in real data) - a flat "mixed" wedge would hide
    # what's actually going on, so break it into its real combinations
    # (same data as mixed_combinations()) and give each its own wedge.
    un_sets, an_sets = _mixed_combo_sets(df)

    # ── plot: two pie charts side by side ──
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Event-Type Distribution — {label}", fontsize=13, fontweight="bold")

    for ax, group_df, group_types, title, combo_sets in [
        (axes[0], unanalyzable_df, UNANALYZABLE_TYPES, "Unanalyzable clusters", un_sets),
        (axes[1], analyzable_df,   ANALYZED_TYPES + ["domain swap"], "Analyzable clusters", an_sets),
    ]:
        counts = group_df["label"].value_counts()
        # keep canonical order, drop types with zero count
        labels = [et for et in group_types if et in counts.index]
        sizes  = [counts[et] for et in labels]
        colors = [event_color(et) for et in labels]

        mixed_labels, mixed_sizes, mixed_colors = _mixed_wedges(combo_sets, len(group_df))
        labels, sizes, colors = (
            [display_label(et) for et in labels] + mixed_labels,
            sizes + mixed_sizes,
            colors + mixed_colors,
        )

        wedges = _draw_event_type_pie(ax, labels, sizes, colors)
        if wedges:
            # Long "reason + reason" combo labels in ncol=3 can grow wider
            # than this subplot and spill into its neighbour - drop to a
            # single column when that's likely.
            ncol = 1 if any(len(l) > 22 for l in labels) else 3
            ax.legend(
                wedges, labels,
                loc="upper center", bbox_to_anchor=(0.5, -0.02),
                fontsize=7.5, ncol=ncol, frameon=False,
            )

        total = sum(sizes)
        ax.set_title(f"{title}\n(n={total})", fontsize=10)

    save(fig, f"event_distribution_{label.lower().replace(' ', '_')}.png")


def _mixed_combo_breakdown(sets, total, threshold_pct=5, max_items=None):
    """
    Break the "mixed" (more-than-one-reason) subset of `sets` (a per-cluster
    Series of frozensets, as returned by _mixed_combo_sets) into its actual
    combinations - each is kept as its own item as long as it's individually
    >= threshold_pct of `total` (the whole pie's cluster count, not just the
    mixed subset, so the returned counts sum correctly with the rest of the
    pie's wedges - a single "mixed" wedge covering 90%+ of a pie isn't
    useful, but most of its *individual* combinations usually are big enough
    on their own to be worth their own slice). Anything under threshold_pct
    folds into one "other mixed (N combos)" residual. max_items caps how
    many combos can be their own item even if they clear the threshold - a
    safety net for the rare case there are more qualifying combos than
    there are colors to give them; the excess folds into the residual too.

    Returns an ordered list of (combo_str, count) pairs, descending by
    count, with the "other mixed" residual (if any) last.
    """
    mixed = sets[(sets.apply(len) > 1) & (sets != DOMAIN_SWAP_COMBO)]
    if mixed.empty:
        return []
    combo_counts = (
        mixed.apply(lambda s: " + ".join(display_label(et) for et in sorted(s, key=event_sort_key)))
        .value_counts()
    )  # already sorted descending
    threshold_count = threshold_pct / 100 * total
    qualifying = combo_counts[combo_counts >= threshold_count]
    if max_items is not None and len(qualifying) > max_items:
        qualifying = qualifying.iloc[:max_items]
    residual = combo_counts.drop(qualifying.index)
    items = list(qualifying.items())
    if len(residual) > 0:
        items.append((f"other mixed ({len(residual)} combos)", int(residual.sum())))
    return items


def _mixed_wedges(sets, total, threshold_pct=5):
    """
    (wedge_labels, wedge_sizes, wedge_colors) for the "mixed" slice of one
    event-type pie, broken into its individual combinations (see
    _mixed_combo_breakdown for the >=threshold_pct selection rule). Shades
    from MIXED_COMBO_SHADES assigned by frequency rank *within this call*,
    so the same shade means a different combo in a different call (e.g.
    Human vs Mouse). Fine for a single pie's own legend
    (event_distribution()); not safe to reuse across panels in one shared
    legend - see _pdf_event_distribution_combined, which assigns colors by
    combo identity instead for that reason.
    """
    items = _mixed_combo_breakdown(sets, total, threshold_pct, max_items=len(MIXED_COMBO_SHADES) - 1)
    if not items:
        return [], [], []
    wedge_labels = [combo for combo, _ in items]
    wedge_sizes  = [count for _, count in items]
    wedge_colors = [MIXED_COMBO_SHADES[min(i, len(MIXED_COMBO_SHADES) - 1)] for i in range(len(items))]
    if wedge_labels[-1].startswith("other mixed"):
        wedge_colors[-1] = MIXED_COMBO_SHADES[-1]
    return wedge_labels, wedge_sizes, wedge_colors


def _draw_event_type_pie(ax, labels, sizes, colors, autopct_threshold=3):
    """Draw one event-type pie's wedges + in-slice %% (no legend) on `ax`. Returns the wedges, or [] if there's no data."""
    if not sizes:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return []
    wedges, _, autotexts = ax.pie(
        sizes, colors=colors,
        autopct=lambda p: f"{p:.1f}%" if p > autopct_threshold else "",
        startangle=90, pctdistance=0.75,
    )
    for t in autotexts:
        t.set_fontsize(7)
    return wedges


def _event_type_pie_data(group_df, group_types):
    """
    (display_labels, sizes, colors) for one event-type pie panel, "mixed"
    left as one plain wedge (not exploded into combinations, unlike
    _mixed_wedges) - every panel then shares the exact same category ->
    color mapping (event_color() is a pure function of the category name),
    so a single legend can validly cover several panels at once. Used by
    the combined multi-species event-distribution page; event_distribution()
    itself still explodes "mixed" for its own single pie, where a per-combo
    color has an unambiguous meaning.
    """
    counts = group_df["label"].value_counts()
    labels = [et for et in group_types if et in counts.index]
    sizes = [counts[et] for et in labels]
    colors = [event_color(et) for et in labels]
    return [display_label(et) for et in labels], sizes, colors


def _mixed_combo_sets(df):
    """
    Per-cluster (specie/cluster) frozenset of the event types actually
    present - for unanalyzable clusters, from all their rows; for analyzable
    clusters, from just the representative transcript's rows (see
    select_representative_transcript). Returns (un_sets, an_sets), each a
    Series covering *every* cluster of that status, not just "mixed" ones -
    filter to `.apply(len) > 1` (and, for an_sets, `!= DOMAIN_SWAP_COMBO`) to
    get only the "mixed" ones, exactly like _collapse_to_cluster_label does.
    Shared by mixed_combinations() (the full breakdown table) and
    event_distribution() (which explodes the pie's "mixed" wedge by combo).
    """
    group_cols = ["specie", "cluster"] if "specie" in df.columns else ["cluster"]
    is_analyzed = df["event_type"].isin(ANALYZED_TYPES)

    has_analyzed = (
        df.assign(_is_analyzed=is_analyzed)
        .groupby(group_cols, observed=True, dropna=False)["_is_analyzed"].any()
    )
    unanalyzable_keys = has_analyzed[~has_analyzed].index
    analyzable_keys = has_analyzed[has_analyzed].index

    un_sets = (
        df[~is_analyzed].groupby(group_cols, observed=True, dropna=False)["event_type"]
        .apply(lambda s: frozenset(s)).reindex(unanalyzable_keys)
    )
    representative_df = select_representative_transcript(df, on_ambiguous="drop")
    an_sets = (
        representative_df.groupby(group_cols, observed=True, dropna=False)["event_type"]
        .apply(lambda s: frozenset(s)).reindex(analyzable_keys).dropna()
    )
    return un_sets, an_sets


def mixed_combinations(df, label, top_n=15):
    """
    "mixed" (see event_distribution/_cluster_event_labels) just says "more
    than one reason" - this breaks down what those reasons actually are, by
    counting how often each exact combination occurs, for unanalyzable
    clusters (all their rows) and analyzable clusters (their representative
    transcript's rows) separately. Prints the top `top_n` combinations for
    each and saves the full breakdown to CSV.

    "domain swap" (added_domain + dropped domain) is excluded from the
    analyzable side here - it's promoted to its own label rather than
    generic "mixed" (see _collapse_to_cluster_label), so it wouldn't be
    counted as "mixed" by event_distribution either; the combination is
    still exactly one thing, not multiple competing reasons to enumerate.
    """
    un_sets, an_sets = _mixed_combo_sets(df)

    print(f"\n{'='*70}")
    print(f"MIXED EVENT-TYPE COMBINATIONS — {label}")
    print(f"{'='*70}")

    csv_rows = []
    for cluster_status, sets in [("unanalyzable", un_sets), ("analyzable", an_sets)]:
        mixed = sets[(sets.apply(len) > 1) & (sets != DOMAIN_SWAP_COMBO)]
        n_total = len(sets)
        n_mixed = len(mixed)
        print(f"\n  {cluster_status.capitalize()}: {n_mixed} mixed cluster(s) "
              f"out of {n_total} ({100*n_mixed/n_total:.1f}%)" if n_total else
              f"\n  {cluster_status.capitalize()}: no clusters.")
        if mixed.empty:
            continue

        combo_counts = mixed.value_counts()  # already sorted descending
        print(f"    {'Combination':<65s} {'Count':>8s} {'% of mixed':>10s}")
        for rank, (combo, count) in enumerate(combo_counts.items()):
            combo_str = " + ".join(sorted(combo))
            pct = 100 * count / n_mixed
            csv_rows.append({"cluster_status": cluster_status, "combination": combo_str,
                              "count": count, "pct_of_mixed": round(pct, 1)})
            if rank < top_n:
                print(f"    {combo_str:<65s} {count:>8d} {pct:>9.1f}%")
        if len(combo_counts) > top_n:
            print(f"    ... and {len(combo_counts) - top_n} more combination(s) - see the CSV below.")

    if csv_rows:
        _save_table(pd.DataFrame(csv_rows), f"mixed_combinations_{label.lower().replace(' ', '_')}.csv",
                    title=f"Mixed Combinations — {label}")


# ══════════════════════════════════════════════════════════════════════════════
# 2. MOST FREQUENT DOMAINS
# ══════════════════════════════════════════════════════════════════════════════

def _interpro_lookup(domain_name):
    """
    Map a results.csv domain_name to the (member_db, accession) pair the
    InterPro REST API expects: /entry/<member_db>/<accession> (see
    https://www.ebi.ac.uk/interpro/api/entry/interpro/IPR011992 for the
    response shape). domain_name is an InterPro accession or a member-
    database signature ID, in whatever casing DoChaP happened to store it -
    both "pfam04503" and "PF01633" show up in real data for the same
    member db. Returns None for a format this doesn't recognize.
    """
    name = str(domain_name).strip()
    if re.fullmatch(r"IPR\d+", name):
        return "interpro", name
    if name.startswith("G3DSA:"):
        return "cathgene3d", name
    if re.fullmatch(r"PF\d+", name):
        return "pfam", name
    m = re.fullmatch(r"(?i:pfam)(\d+)", name)
    if m:
        return "pfam", f"PF{m.group(1)}"
    if re.fullmatch(r"SM\d+", name):
        return "smart", name
    m = re.fullmatch(r"(?i:smart)(\d+)", name)
    if m:
        return "smart", f"SM{m.group(1)}"
    if re.fullmatch(r"TIGR\d+", name):
        return "ncbifam", name  # TIGRFAM was folded into NCBIFAM in InterPro
    m = re.fullmatch(r"(?i:tigr)(\d+)", name)
    if m:
        return "ncbifam", f"TIGR{m.group(1)}"
    if re.fullmatch(r"cd\d+", name):
        return "cdd", name
    if re.fullmatch(r"PTHR\d+", name):
        return "panther", name
    if re.fullmatch(r"SSF\d+", name):
        return "ssf", name
    if re.fullmatch(r"PS\d+", name):
        return "profile", name
    return None


def _clean_interpro_text(text):
    """Normalize an InterPro description to plain text.

    Handles both the raw REST payload (HTML tags + inline [[cite:...]] markers)
    and the DoChaP RepresentativeDomains.description text, which has already had
    its citations removed but left behind empty reference brackets like
    "[ , , ]". Strips all of those, tidies the space a removed citation leaves
    before punctuation (e.g. "RNAs ." -> "RNAs."), and collapses whitespace.
    """
    text = re.sub(r"\[\[cite:.*?\]\]", "", text)      # inline citation markers
    text = re.sub(r"<[^>]+>", "", text)                # HTML tags
    text = re.sub(r"\[\s*(?:,\s*)*\]", "", text)       # empty "[ , , ]" citation remnants
    text = re.sub(r"\s+([.,;:])", r"\1", text)          # space left before punctuation
    return re.sub(r"\s+", " ", text).strip()


def _web_description_fallback(missing_domains, result_map, min_interval=0.4):
    """Fill descriptions for `missing_domains` from the InterPro REST API, in place.

    Best-effort and fully guarded: on a server where outbound web access is
    blocked (the production case), this must never raise or hang. A quick
    reachability probe short-circuits a firewalled host in one short timeout
    instead of retrying per domain, and the whole thing is wrapped in try/except
    so any network failure just leaves those domains without a description.

    Returns the number of domains it managed to fill.
    """
    import urllib.error
    import urllib.request

    n_web = 0
    try:
        # Probe once: an HTTPError means the host answered (reachable); a bare
        # URLError means it's unreachable/blocked -> skip the whole fallback.
        try:
            urllib.request.urlopen("https://www.ebi.ac.uk/interpro/api/", timeout=5)
        except urllib.error.HTTPError:
            pass

        from check_db import RateLimitedSession
        # max_retries=0 so a blocked host costs one short timeout per call, not a
        # retry storm; timeout kept modest for the same reason.
        session = RateLimitedSession(min_interval=min_interval, max_retries=0, timeout=10)
        for domain_name in missing_domains:
            lookup = _interpro_lookup(domain_name)
            if lookup is None:
                continue
            member_db, accession = lookup
            url = f"https://www.ebi.ac.uk/interpro/api/entry/{member_db}/{accession}"
            metadata = (session.get_json(url) or {}).get("metadata")
            if metadata is None:
                continue
            blocks = metadata.get("description") or []
            web_desc = _clean_interpro_text(blocks[0]["text"]) if blocks else None
            if web_desc:
                result_map[domain_name][0] = result_map[domain_name][0] or (metadata.get("name") or {}).get("name")
                result_map[domain_name][1] = web_desc
                n_web += 1
    except Exception as exc:  # noqa: BLE001 - network access is best-effort
        print(f"  web fallback unavailable ({type(exc).__name__}: {exc}); "
              f"{len(missing_domains)} domain(s) left without a description.")
    return n_web


def fetch_domain_descriptions(domain_names, out_csv=None, con=None, db_path=None, web_fallback=True):
    """
    Look up each top domain's name + description from the DoChaP
    RepresentativeDomains table (its `domain_name` and `description` columns),
    keyed by `domain_id` - the accession that results.csv stores in its own
    `domain_name` column. Returns a DataFrame [domain_name, name, description],
    one row per unique input domain.

    This replaces the previous InterPro REST lookup as the *primary* source: the
    descriptions now come from the local DB instead of the network. For the ~26%
    of domains whose RepresentativeDomains.description is empty (and any domain
    not in that table at all), an optional InterPro web fallback (web_fallback,
    default on) tries to fill the gap - but it is fully guarded so it never
    raises or hangs where outbound access is blocked (see
    _web_description_fallback). Domains still unresolved get an empty
    description rather than raising.

    Pass an already-open sqlite3 connection as `con`, or a `db_path` to open
    (defaults to DEFAULT_DOCHAP_DB_PATH).
    """
    domain_names = pd.unique(pd.Series(list(domain_names)).dropna())

    own_con = con is None
    if own_con:
        con = sqlite3.connect(db_path or DEFAULT_DOCHAP_DB_PATH)
    try:
        lookup = {}
        if len(domain_names):
            placeholders = ",".join("?" * len(domain_names))
            # domain_name/description are constant per domain_id, so GROUP BY
            # collapses the many per-protein rows down to one per accession.
            query = (
                "SELECT domain_id, domain_name, description "
                "FROM RepresentativeDomains "
                f"WHERE domain_id IN ({placeholders}) GROUP BY domain_id"
            )
            for domain_id, name, description in con.execute(query, list(domain_names)):
                cleaned = _clean_interpro_text(description) if description else None
                lookup[domain_id] = (name, cleaned or None)
    finally:
        if own_con:
            con.close()

    # result_map: domain_name -> [name, description]; mutated in place by the fallback.
    result_map = {d: list(lookup.get(d, (None, None))) for d in domain_names}
    n_db = sum(1 for d in domain_names if result_map[d][1])

    missing = [d for d in domain_names if not result_map[d][1]]
    n_web = _web_description_fallback(missing, result_map) if (web_fallback and missing) else 0

    result = pd.DataFrame(
        [{"domain_name": d, "name": result_map[d][0], "description": result_map[d][1]} for d in domain_names]
    )
    print(f"  descriptions: {n_db} from RepresentativeDomains"
          f"{f', {n_web} via InterPro web fallback' if web_fallback else ''} "
          f"({len(result)} domains total).")
    if out_csv:
        result.to_csv(out_csv, index=False)
        print(f"  → {out_csv}")
    return result


def domain_frequency(df, label, top_n=20):
    """
    Show the top_n most frequent domains overall, then break each down by
    which event types they appear in. Returns the top_domains Series
    (domain_name -> count), or None if there's no domain data - callers
    that want InterPro descriptions for these domains (see
    fetch_domain_descriptions) do that separately, rather than this
    function doing it inline, so that table can be positioned wherever
    makes sense in a report rather than always right after this chart.
    """
    # Only rows where a domain is recorded
    domain_df = df[df["domain_name"].notna() & (df["domain_name"] != "")]

    if domain_df.empty:
        print(f"\nNo domain data found for {label}.")
        return None

    top_domains = domain_df["domain_name"].value_counts().head(top_n)

    print(f"\n{'='*70}")
    print(f"TOP {top_n} DOMAINS — {label}")
    print(f"{'='*70}")
    for domain, count in top_domains.items():
        breakdown = domain_df[domain_df["domain_name"] == domain]["event_type"].value_counts()
        breakdown_str = ", ".join(f"{et}: {c}" for et, c in breakdown.items())
        print(f"  {domain:<30s}  {count:3d}  [{breakdown_str}]")

    # Per-species rank (not capped to top_n) for each of these domains, so a
    # reader can see whether a domain's overall ranking is driven mostly by
    # one species or is fairly common to both - '-' if it doesn't occur for
    # that species at all. Only meaningful with 2+ species actually present
    # (skipped for IOE, which has no specie column, and for an
    # already-species-filtered df like "Hadas Human", where every row is
    # already that one species and a rank-comparison is meaningless).
    species_ranks = None
    if "specie" in domain_df.columns:
        species_present = sorted(s for s in domain_df["specie"].dropna().unique() if s)
        if len(species_present) > 1:
            species_ranks = {}
            for sp in species_present:
                sp_rank_of = {
                    d: i + 1 for i, d in enumerate(domain_df[domain_df["specie"] == sp]["domain_name"]
                                                    .value_counts().index)
                }
                species_ranks[sp] = [sp_rank_of.get(d) for d in top_domains.index]

    # ── plot: horizontal bar chart coloured by event type ──
    left_margin = 0.34 if species_ranks else 0.22
    fig, ax = plt.subplots(figsize=(10, max(4, top_n * 0.5)))
    fig.subplots_adjust(left=left_margin)
    fig.suptitle(f"Top {top_n} Domains — {label}", fontsize=13, fontweight="bold")
    fig.text(0.5, 0.955, "Counts occurrences: one bar segment per matched domain-comparison event"
                          " (colored by outcome, human+mouse combined - see the per-species cluster"
                          " count below).",
             ha="center", fontsize=8, color="#666666", style="italic")

    # Stacked bars: each segment is an event type
    y = np.arange(len(top_domains))
    lefts = np.zeros(len(top_domains))
    event_types_present = sorted(
        domain_df[domain_df["domain_name"].isin(top_domains.index)]["event_type"].unique(),
        key=event_sort_key
    )

    for et in event_types_present:
        widths = [
            (domain_df[(domain_df["domain_name"] == d) & (domain_df["event_type"] == et)].shape[0])
            for d in top_domains.index
        ]
        ax.barh(y, widths, left=lefts, color=event_color(et), label=display_label(et), edgecolor="white")
        lefts += np.array(widths)

    ax.set_yticks(y)
    ax.set_yticklabels(top_domains.index, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Number of occurrences")
    ax.legend(loc="lower right", fontsize=8)

    if species_ranks:
        # x in axes-fraction (so negative values land in the left margin,
        # to the left of x=0 where the plot area / bars start), y in data
        # coordinates (so rows line up with each domain's bar) - the
        # standard matplotlib trick for annotating outside the axes proper.
        # Text color (not a label/legend) carries the human/mouse distinction,
        # matching COLORS elsewhere in this report - a one-letter header
        # above each column spells out which is which.
        trans = ax.get_yaxis_transform()
        col_x = {sp: -0.30 + i * 0.11 for i, sp in enumerate(species_ranks)}
        for sp, x in col_x.items():
            ax.text(x, -1.0, sp[0].upper(), transform=trans, color=COLORS.get(sp, "#333"),
                    fontsize=8, fontweight="bold", ha="center", va="bottom", clip_on=False)
        for row_i in range(len(top_domains)):
            for sp, x in col_x.items():
                rank = species_ranks[sp][row_i]
                ax.text(x, row_i, str(rank) if rank is not None else "-", transform=trans,
                        color=COLORS.get(sp, "#333"), fontsize=8, ha="center", va="center", clip_on=False)

    save(fig, f"domain_frequency_{label.lower().replace(' ', '_')}.png")
    return top_domains


def _species_domain_counts(domain_df, specie):
    """domain_name -> occurrence count for one specie's rows, domains with 0 occurrences excluded.

    value_counts() on a categorical column includes every category at count
    0 (even ones absent from this specie's rows entirely, since domain_name
    carries every category seen across the whole df) - drop those so index
    membership below means "actually occurs for this specie", not "exists as
    a category somewhere in the file".
    """
    counts = domain_df[domain_df["specie"] == specie]["domain_name"].value_counts()
    return counts[counts > 0]


def domain_species_rank_scatter(df, label, top_domains=None):
    """
    Two scatter plots comparing per-domain frequency between species: rank
    (1 = most frequent) and occurrence count, one species vs the other. The
    count axes are normalized to % of that species' total domain
    occurrences rather than raw counts, since the two species' domain row
    counts are rarely comparable (Hadas human routinely has several times
    more domain rows than mouse) - plotting raw counts would put nearly
    every point below the y=x line purely from that volume difference, not
    because the domain is actually more common in one species biologically.

    Only meaningful with exactly one species *pair* present (same
    precondition as domain_frequency()'s per-domain rank columns) - skipped
    with a warning for IOE (no 'specie' column), an already species-filtered
    df (e.g. "Hadas Human", from analyze_file(specie_filter=...), which by
    definition only has one species), or more than 2 species.

    Rendered twice: once across every domain the two species share, and once
    restricted to top_domains (domain_frequency()'s own top-N Series, if
    given) - so the report also has a version lined up with the bar chart
    directly above it.
    """
    domain_df = df[df["domain_name"].notna() & (df["domain_name"] != "")]
    if domain_df.empty:
        return
    if "specie" not in domain_df.columns:
        _warn(f"domain_species_rank_scatter: {label} has no 'specie' column - skipping.")
        return

    species_present = sorted(s for s in domain_df["specie"].dropna().unique() if s)
    if len(species_present) != 2:
        _warn(f"domain_species_rank_scatter: {label} has {len(species_present)} specie(s) "
              f"({species_present}), need exactly 2 - skipping.")
        return
    sp1, sp2 = species_present

    counts1 = _species_domain_counts(domain_df, sp1)
    counts2 = _species_domain_counts(domain_df, sp2)
    if counts1.empty or counts2.empty:
        return
    rank1 = {d: i + 1 for i, d in enumerate(counts1.index)}
    rank2 = {d: i + 1 for i, d in enumerate(counts2.index)}
    # Denominators for the count axes: total domain occurrences for that
    # species (not just the domains shared with the other species), so a
    # domain's normalized value reads as "% of every domain occurrence in
    # this species", comparable across species regardless of how much more
    # domain data one species has overall.
    total1 = int(counts1.sum())
    total2 = int(counts2.sum())
    shared_keys = set(counts1.index) & set(counts2.index)
    color = COLORS.get(label, "#4472C4")
    label_slug = label.lower().replace(' ', '_')

    def _render(domains, suffix, scope_label):
        keys = shared_keys if domains is None else shared_keys & set(domains)
        common = sorted(keys)
        if len(common) < 2:
            print(f"  Skipping domain {sp1}/{sp2} scatter ({scope_label}): "
                  f"only {len(common)} domain(s) present in both species.")
            return

        x_rank = np.array([rank1[d] for d in common])
        y_rank = np.array([rank2[d] for d in common])
        x_pct = np.array([counts1[d] for d in common]) / total1 * 100
        y_pct = np.array([counts2[d] for d in common]) / total2 * 100
        corr_rank = np.corrcoef(x_rank, y_rank)[0, 1]
        # Pearson r is scale-invariant, so this is identical whether computed
        # on the normalized % or the raw counts - only the axes/y=x meaning change.
        corr_count = np.corrcoef(x_pct, y_pct)[0, 1]

        fig, ax = plt.subplots(figsize=(7, 7))
        ax.scatter(x_rank, y_rank, s=18, alpha=0.55, color=color, edgecolor="white", linewidth=0.3)
        max_rank = max(x_rank.max(), y_rank.max())
        ax.plot([1, max_rank], [1, max_rank], linestyle="--", color="#999999", linewidth=1, label="y = x (equal rank)")
        ax.set_xlabel(f"Domain rank in {sp1} (1 = most frequent)")
        ax.set_ylabel(f"Domain rank in {sp2} (1 = most frequent)")
        ax.set_title(f"Domain frequency rank: {sp1} vs {sp2} — {label}\n"
                     f"({scope_label}, {len(common)} domains, Pearson r = {corr_rank:.2f})",
                     fontsize=11, fontweight="bold")
        ax.legend(loc="upper left", fontsize=8)
        ax.set_xlim(0, max_rank * 1.02)
        ax.set_ylim(0, max_rank * 1.02)
        fig.tight_layout()
        save(fig, f"domain_rank_scatter{suffix}_{label_slug}.png")

        fig, ax = plt.subplots(figsize=(7, 7))
        ax.scatter(x_pct, y_pct, s=18, alpha=0.55, color=color, edgecolor="white", linewidth=0.3)
        max_pct = max(x_pct.max(), y_pct.max(), 1e-9)
        ax.plot([0, max_pct], [0, max_pct], linestyle="--", color="#999999", linewidth=1,
                label="y = x (equal relative frequency)")
        ax.set_xlabel(f"% of all domain occurrences in {sp1}  (n={total1:,})")
        ax.set_ylabel(f"% of all domain occurrences in {sp2}  (n={total2:,})")
        ax.set_title(f"Domain occurrence frequency: {sp1} vs {sp2} — {label}\n"
                     f"({scope_label}, {len(common)} domains, Pearson r = {corr_count:.2f})",
                     fontsize=11, fontweight="bold")
        ax.legend(loc="upper left", fontsize=8)
        ax.set_xlim(0, max_pct * 1.05)
        ax.set_ylim(0, max_pct * 1.05)
        fig.tight_layout()
        save(fig, f"domain_count_scatter{suffix}_{label_slug}.png")

    _render(None, "", "all shared domains")
    if top_domains is not None and len(top_domains):
        _render(set(top_domains.index), "_top20", f"top {len(top_domains)} domains")


def domain_cluster_species_counts(df, top_domains, label):
    """
    For domain_frequency()'s same top_domains, how many distinct clusters
    each one appears in for human vs mouse - a different question from
    domain_frequency's own bar lengths, which count rows (one per matched
    domain-comparison event on the representative transcript). A domain
    could rack up a lot of *rows* while actually being concentrated in just
    a handful of clusters (e.g. if a cluster's comparison happens to yield
    the same domain name more than once), or be spread across many - this
    shows cluster spread specifically, side by side by species.

    Skipped (no chart) without a 'specie' column or with fewer than 2
    species present - e.g. IOE, or an already-species-filtered df like
    "Hadas Human", where a per-species comparison is meaningless.
    """
    if top_domains is None or top_domains.empty:
        return
    if "specie" not in df.columns:
        return
    species_present = sorted(s for s in df["specie"].dropna().unique() if s)
    if len(species_present) < 2:
        return

    domain_df = df[df["domain_name"].isin(top_domains.index)]
    cluster_counts = {
        sp: domain_df[domain_df["specie"] == sp].groupby("domain_name")["cluster"]
            .nunique().reindex(top_domains.index, fill_value=0)
        for sp in species_present
    }

    print(f"\n{'='*70}")
    print(f"TOP {len(top_domains)} DOMAINS — CLUSTER COUNT BY SPECIES — {label}")
    print(f"{'='*70}")
    header = "".join(f"{sp.capitalize():>10s}" for sp in species_present)
    print(f"  {'Domain':<15s} {header}")
    for d in top_domains.index:
        row = "".join(f"{cluster_counts[sp][d]:>10d}" for sp in species_present)
        print(f"  {d:<15s} {row}")

    # ── plot: paired horizontal bars, one pair per domain ──
    fig, ax = plt.subplots(figsize=(10, max(4, len(top_domains) * 0.5)))
    fig.suptitle(f"Top {len(top_domains)} Domains — Cluster Count by Species — {label}",
                 fontsize=13, fontweight="bold")
    fig.text(0.5, 0.955, "Counts distinct clusters per domain, split by species (the chart above"
                          " combines both species and counts occurrences, not clusters).",
             ha="center", fontsize=8, color="#666666", style="italic")

    y = np.arange(len(top_domains))
    bar_h = 0.7 / len(species_present)
    for i, sp in enumerate(species_present):
        offset = (i - (len(species_present) - 1) / 2) * bar_h
        ax.barh(y + offset, cluster_counts[sp].values, height=bar_h,
                color=COLORS.get(sp, "#999"), label=sp.capitalize())

    ax.set_yticks(y)
    ax.set_yticklabels(top_domains.index, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Number of clusters")
    ax.legend(loc="lower right", fontsize=8)

    save(fig, f"domain_cluster_species_counts_{label.lower().replace(' ', '_')}.png")


def _fetch_and_save_domain_descriptions(top_domains, label):
    """
    Small wrapper so _run_analyses() can look up + save descriptions for
    domain_frequency()'s top_domains without directly calling the
    module-level fetch_domain_descriptions() from inside a scope where a
    same-named boolean parameter (fetch_domain_descriptions=True/False)
    would shadow it.
    """
    result = fetch_domain_descriptions(top_domains.index)
    _save_table(result, f"domain_descriptions_{label.lower().replace(' ', '_')}.csv",
                title=f"Top Domains — {label}")


# ══════════════════════════════════════════════════════════════════════════════
# 3/4. LENGTH CHANGE DISTRIBUTION — shorter and longer analyzed separately
# ══════════════════════════════════════════════════════════════════════════════

def _length_change_data(df, event_type):
    """
    (abs_change, rel_change) Series - magnitude-only domain length change
    for one event type ("shorter"/"longer"), c_domain_length > 0 rows only
    (needed to compute a relative %%). Both empty if there's no matching
    data. Shared by _length_change() (single-file report) and
    _pdf_length_change_combined() (multi-species overview page).
    """
    change_df = df[
        (df["event_type"] == event_type)
        & df["c_domain_length"].notna()
        & df["t_domain_length"].notna()
    ]
    change_df = change_df[change_df["c_domain_length"] > 0]
    if change_df.empty:
        return pd.Series(dtype=float), pd.Series(dtype=float)
    abs_change = (change_df["c_domain_length"] - change_df["t_domain_length"]).abs()
    rel_change = 100 * abs_change / change_df["c_domain_length"]
    return abs_change, rel_change


def _draw_length_change_hist(ax, vals, color, xlabel, title):
    """
    Draw one length-change histogram (log-spaced bins - finer resolution at
    small changes, where most events sit, coarser at the long sparse tail)
    on `ax`. `vals` are magnitudes (always >= 0); a real length change can't
    be exactly 0, and log bins can't include 0 anyway, so those are dropped.
    """
    vals = vals.dropna()
    vals = vals[vals > 0]
    if vals.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        return
    if vals.min() == vals.max():
        # A single unique value (often just one data point) can't seed a
        # log-spaced bin range - geomspace(x, x, n) is constant, not
        # monotonically increasing, which np.histogram/ax.hist rejects outright.
        ax.text(0.5, 0.5, f"All {len(vals)} value(s) = {vals.iloc[0]:.1f}",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_xlabel(xlabel)
        ax.set_title(title)
        return
    bins = np.geomspace(vals.min(), vals.max(), 30)
    ax.hist(vals.values, bins=bins, alpha=0.85, color=color, edgecolor="white")
    ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.set_title(title)


def _length_change(df, label, event_type):
    """
    Plot and print the distribution of domain length change for one event
    type ("shorter" or "longer") - both as magnitudes (always >= 0), not a
    signed c-t delta: the direction is already fixed by which of these two
    this call is for, so a sign would just relabel a positive number as
    negative on the "longer" plot without adding information.
    """
    abs_change, rel_change = _length_change_data(df, event_type)
    if abs_change.empty:
        print(f"\nNo {event_type} rows with domain lengths found for {label}.")
        return

    print(f"\n{'='*70}")
    print(f"LENGTH CHANGE DISTRIBUTION — {event_type.upper()} — {label}  (n={len(abs_change)})")
    print(f"{'='*70}")
    for vals, name, unit in [(abs_change, "Absolute Δaa", "aa"), (rel_change, "Relative Δ%", "%")]:
        print(f"  {name:<14s}  mean: {vals.mean():.1f}{unit}  "
              f"median: {vals.median():.1f}{unit}  range: [{vals.min():.0f}, {vals.max():.0f}]{unit}")

    color = LENGTH_CHANGE_COLORS[event_type]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f"Domain Length Change — {event_type.capitalize()} — {label}",
                 fontsize=13, fontweight="bold")
    _draw_length_change_hist(axes[0], abs_change, color, "Δ amino acids (magnitude)", "Absolute change (aa)")
    _draw_length_change_hist(axes[1], rel_change, color, "Δ % of canonical domain length", "Relative change (%)")

    save(fig, f"length_change_{event_type}_{label.lower().replace(' ', '_')}.png")


def _pdf_length_change_combined(pdf, label_dfs, event_type, page_title):
    """
    One page: a 2 x len(label_dfs) grid of length-change histograms for one
    event type ("shorter"/"longer") - row 1 = absolute Δaa, row 2 =
    relative Δ%%, one column per label - same aligned-grid rationale as
    _pdf_event_distribution_combined. No shared legend needed (plain
    histograms, not categorical).

    Runs on select_representative_transcript()'s output (one transcript per
    cluster), same as every other analysis of analyzed-outcome distributions
    (domain_frequency, domain_count_change, severity_spectrum,
    splice_type_vs_outcome, and length_change_shorter/longer itself when
    not combined) - a cluster with several comparable transcripts otherwise
    gets one histogram entry per transcript instead of one for the cluster,
    overweighting exactly the clusters that happen to have the most
    candidate transcripts.
    """
    n = len(label_dfs)
    fig, axes = plt.subplots(2, n, figsize=(4.3 * n, 8.5), squeeze=False)
    fig.suptitle(page_title, fontsize=14, fontweight="bold")
    color = LENGTH_CHANGE_COLORS[event_type]

    for col_idx, (label, df) in enumerate(label_dfs.items()):
        representative_df = select_representative_transcript(df, on_ambiguous="drop")
        abs_change, rel_change = _length_change_data(representative_df, event_type)
        _draw_length_change_hist(axes[0, col_idx], abs_change, color,
                                  "Δ amino acids (magnitude)", f"{label} — Absolute (n={len(abs_change)})")
        _draw_length_change_hist(axes[1, col_idx], rel_change, color,
                                  "Δ % of canonical length", "Relative (%)")

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    pdf.savefig(fig, dpi=200)
    plt.close(fig)


def length_change_shorter(df, label):
    _length_change(df, label, "shorter")


def length_change_longer(df, label):
    _length_change(df, label, "longer")


# ══════════════════════════════════════════════════════════════════════════════
# 5. SPECIES COMPARISON (OVERALL + PER GENE)
# ══════════════════════════════════════════════════════════════════════════════

def species_comparison(df, label):
    """
    Compare species within a single file (e.g. human vs mouse in Hadas):
      - Overall event-type distribution
      - Per-gene event-type breakdown → CSV + plot
    Skipped (with a warning) if there's no real species data - e.g. an IOE
    file, which has no specie column at all.

    Both use one result per cluster, not per row (see _cluster_event_labels)
    - the same "cluster-level vs per-transcript-comparison" mixing problem
    event_distribution() has applies here too, split by species.
    """
    species_list = sorted(s for s in df["specie"].dropna().unique() if s)
    if not species_list:
        _warn(f"skipping species comparison for {label} — no species data.")
        return

    cluster_labels = _cluster_event_labels(df)
    cluster_labels = cluster_labels[cluster_labels["specie"].isin(species_list)]

    # ── overall ──
    print(f"\n{'='*70}")
    print(f"{label} — SPECIES COMPARISON (OVERALL)")
    print(f"{'='*70}")

    for sp in species_list:
        sp_df      = cluster_labels[cluster_labels["specie"] == sp]
        un_df      = sp_df[sp_df["cluster_status"] == "unanalyzable"]
        an_df      = sp_df[sp_df["cluster_status"] == "analyzable"]
        n_total    = len(sp_df)
        n_unanalyz = len(un_df)
        n_analyzed = len(an_df)

        print(f"\n  {sp.upper()}  (n={n_total} clusters)")
        print(f"    Unanalyzable : {n_unanalyz}  ({100*n_unanalyz/n_total:.1f}%)")
        for et in UNANALYZABLE_TYPES + ["mixed"]:
            count = (un_df["label"] == et).sum()
            if count:
                print(f"      {et:<40s} {count}  ({100*count/n_unanalyz:.1f}%)")
        print(f"    Analyzable   : {n_analyzed}  ({100*n_analyzed/n_total:.1f}%)")
        for et in ANALYZED_TYPES + ["domain swap", "mixed"]:
            count = (an_df["label"] == et).sum()
            if count:
                print(f"      {et:<40s} {count}  ({100*count/n_analyzed:.1f}%)")

    # overall plot
    panels = [
        ("unanalyzable", UNANALYZABLE_TYPES + ["mixed"], "Unanalyzable clusters"),
        ("analyzable",   ANALYZED_TYPES + ["domain swap", "mixed"],     "Analyzable clusters"),
    ]
    fig_width = sum(max(6, 1.0 * len(group_types)) for _, group_types, _ in panels)
    fig, axes = plt.subplots(1, 2, figsize=(fig_width, 5.5))
    fig.suptitle(f"{label} — Event-Type Distribution by Species", fontsize=13, fontweight="bold")

    for ax, (cluster_status, group_types, title) in zip(axes, panels):
        status_df = cluster_labels[cluster_labels["cluster_status"] == cluster_status]
        x = np.arange(len(group_types))
        w = 0.35

        for i, sp in enumerate(species_list):
            sp_df      = status_df[status_df["specie"] == sp]
            n_in_group = len(sp_df)
            pcts = [
                100 * (sp_df["label"] == et).sum() / n_in_group if n_in_group else 0
                for et in group_types
            ]
            ax.bar(x + i * w, pcts, w, label=sp.capitalize(),
                   color=COLORS.get(sp, "#999"), alpha=0.85)

        ax.set_xticks(x + w / 2)
        ax.set_xticklabels([display_label(et) for et in group_types], fontsize=9, rotation=30, ha="right")
        ax.set_ylabel("% within group")
        ax.set_title(title)
        ax.legend()
        ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))

    fig.tight_layout()
    save(fig, f"species_overall_{label.lower().replace(' ', '_')}.png")

    # ── per gene ──
    if "gene_symbol" not in cluster_labels.columns:
        _warn(f"skipping per-gene species breakdown for {label} — no gene_symbol column.")
        return

    genes = sorted(cluster_labels["gene_symbol"].dropna().unique())

    # Build a tidy pivot: rows = (gene, label), columns = species
    pivot = (
        cluster_labels.groupby(["gene_symbol", "label", "specie"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    # Add % columns for each species
    for sp in species_list:
        gene_totals = cluster_labels.groupby("gene_symbol")["specie"].apply(
            lambda s: (s == sp).sum()
        ).rename(f"{sp}_total")
        pivot = pivot.merge(gene_totals, on="gene_symbol", how="left")
        pivot[f"{sp}_pct"] = (
            pivot.get(sp, 0) / pivot[f"{sp}_total"].replace(0, np.nan) * 100
        ).round(1).fillna("").astype(str).apply(lambda v: f"{v}%" if v else "")

    # Rename raw count columns to <species>_count
    count_cols = {sp: f"{sp}_count" for sp in species_list if sp in pivot.columns}
    pivot = pivot.rename(columns=count_cols)
    # Drop total helper columns
    pivot = pivot.drop(columns=[f"{sp}_total" for sp in species_list], errors="ignore")
    # Sort rows by gene and canonical event order
    pivot["_sort"] = pivot["label"].map(event_sort_key)
    pivot = pivot.sort_values(["gene_symbol", "_sort"]).drop(columns="_sort")

    _save_table(pivot, f"species_per_gene_comparison_{label.lower().replace(' ', '_')}.csv",
                title=f"Per-Gene Comparison — {label}")

    # per-gene plot - one page per fixed-size chunk of genes, each its own
    # properly-sized figure (not one giant figure needing pixel-based
    # splitting afterward, which could cut a panel in half depending on
    # exactly where a page boundary happened to fall in raw pixels, and
    # didn't guarantee the same panel count on every page). Cap to the top
    # genes by cluster count; the CSV above always has every gene.
    max_genes_plotted = 20
    plotted_genes = genes
    if len(genes) > max_genes_plotted:
        plotted_genes = sorted(cluster_labels["gene_symbol"].value_counts().head(max_genes_plotted).index)
        print(f"  NOTE: {len(genes)} genes present — plotting only the top {max_genes_plotted} "
              f"by cluster count (full breakdown in the CSV above).")

    genes_per_page = 4
    gene_chunks = [plotted_genes[i:i + genes_per_page] for i in range(0, len(plotted_genes), genes_per_page)]

    for page_num, gene_chunk in enumerate(gene_chunks, start=1):
        fig, axes = plt.subplots(1, len(gene_chunk), figsize=(2.75 * len(gene_chunk), 5), squeeze=False)
        page_note = f" (page {page_num}/{len(gene_chunks)})" if len(gene_chunks) > 1 else ""
        fig.suptitle(f"{label} — Per-Gene Event Types by Species{page_note}", fontsize=13, fontweight="bold")

        for ax, gene in zip(axes[0], gene_chunk):
            gene_df = cluster_labels[cluster_labels["gene_symbol"] == gene]
            labels_present = sorted(gene_df["label"].unique(), key=event_sort_key)
            x = np.arange(len(labels_present))
            w = 0.35

            for i, sp in enumerate(species_list):
                sp_gene_df = gene_df[gene_df["specie"] == sp]
                gene_total = len(sp_gene_df)
                pcts = [
                    100 * (sp_gene_df["label"] == et).sum() / gene_total if gene_total else 0
                    for et in labels_present
                ]
                ax.bar(x + i * w, pcts, w, label=sp.capitalize(),
                       color=COLORS.get(sp, "#999"), alpha=0.85)

            ax.set_title(gene, fontweight="bold")
            ax.set_xticks(x + w / 2)
            ax.set_xticklabels([display_label(et) for et in labels_present], fontsize=8, rotation=30, ha="right")
            ax.set_ylabel("% of gene total")
            ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
            ax.legend(fontsize=8)

        fig.tight_layout(rect=[0, 0, 1, 0.93])
        save(fig, f"species_per_gene_{label.lower().replace(' ', '_')}_page{page_num}.png")


# ══════════════════════════════════════════════════════════════════════════════
# 6. AS SPLICE TYPE vs DOMAIN OUTCOME
# ══════════════════════════════════════════════════════════════════════════════

def splice_type_vs_outcome(df, label):
    """
    Use the AS event type encoded in the cluster name prefix (A3, A5, SE, RI,
    AF, AL, MX) to ask: does the type of alternative splicing predict the
    domain outcome? Skipped (with a warning) if no cluster in this file
    carries an AS-type prefix (e.g. Hadas, whose cluster ids look like
    "chr11:clu_6611").
    """
    as_df = df[df["as_type"] != ""]
    if as_df.empty:
        _warn(f"skipping AS splice type vs outcome for {label} — no AS-type information in cluster ids.")
        return

    as_types      = sorted(as_df["as_type"].unique())
    outcome_types = [et for et in ANALYZED_TYPES if et in as_df["event_type"].values]

    print(f"\n{'='*70}")
    print(f"AS SPLICE TYPE vs DOMAIN OUTCOME — {label}")
    print(f"{'='*70}")

    for as_type in as_types:
        sub = as_df[as_df["as_type"] == as_type]
        n_analyzed = sub["event_type"].isin(ANALYZED_TYPES).sum()
        outcome_counts = sub[sub["event_type"].isin(ANALYZED_TYPES)]["event_type"].value_counts()
        outcome_str = ", ".join(f"{et}: {outcome_counts.get(et,0)}" for et in outcome_types)
        print(f"  {as_type:<6s}  n={len(sub):3d}  analyzed={n_analyzed}  [{outcome_str}]")

    # ── heatmap: AS type × outcome, cells = % of AS-type rows that are that outcome ──
    matrix = np.zeros((len(as_types), len(outcome_types)))
    for i, as_type in enumerate(as_types):
        sub     = as_df[as_df["as_type"] == as_type]
        n_total = len(sub)
        for j, et in enumerate(outcome_types):
            matrix[i, j] = 100 * (sub["event_type"] == et).sum() / n_total if n_total else 0

    fig, ax = plt.subplots(figsize=(8, max(3, len(as_types) * 0.7)))
    fig.suptitle(f"AS Splice Type vs Domain Outcome (% of AS-type rows) — {label}",
                 fontsize=12, fontweight="bold")

    im = ax.imshow(matrix, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(outcome_types)))
    ax.set_xticklabels([display_label(et) for et in outcome_types], rotation=25, ha="right", fontsize=9)
    ax.set_yticks(range(len(as_types)))
    ax.set_yticklabels(as_types)

    # Annotate each cell with its value
    for i in range(len(as_types)):
        for j in range(len(outcome_types)):
            val = matrix[i, j]
            if val > 0:
                ax.text(j, i, f"{val:.0f}%", ha="center", va="center",
                        fontsize=8, color="black" if val < 50 else "white")

    plt.colorbar(im, ax=ax, label="% of AS-type rows")
    save(fig, f"splice_type_vs_outcome_{label.lower().replace(' ', '_')}.png")


# ══════════════════════════════════════════════════════════════════════════════
# 7. DOMAIN COUNT CHANGE
# ══════════════════════════════════════════════════════════════════════════════

def domain_count_change(df, label):
    """
    For analyzed events, compare the number of domain copies in the canonical
    transcript vs the alternative transcript.
    Categories: fewer copies, same number, more copies.
    """
    analyzed_df = df[
        df["event_type"].isin(ANALYZED_TYPES)
        & df["c_domains_number"].notna()
        & df["t_domains_number"].notna()
    ].copy()

    if analyzed_df.empty:
        print(f"\nNo domain count data found for {label}.")
        return

    analyzed_df["count_change"] = analyzed_df["t_domains_number"] - analyzed_df["c_domains_number"]

    def categorise(delta):
        if delta < 0:
            return "fewer copies"
        elif delta == 0:
            return "same count"
        else:
            return "more copies"

    analyzed_df["count_category"] = analyzed_df["count_change"].apply(categorise)
    category_order = ["fewer copies", "same count", "more copies"]

    print(f"\n{'='*70}")
    print(f"DOMAIN COUNT CHANGE — {label}")
    print(f"{'='*70}")
    for cat in category_order:
        count = (analyzed_df["count_category"] == cat).sum()
        pct   = 100 * count / len(analyzed_df)
        print(f"  {cat:<20s} {count:3d}  ({pct:.1f}%)")

    # ── plot: stacked bars per event type ──
    event_types_present = sorted(analyzed_df["event_type"].unique(), key=event_sort_key)
    fig, ax = plt.subplots(figsize=(max(9, 1.1 * len(event_types_present)), 4.5))
    fig.suptitle(f"Domain Copy-Count Change per Event Type — {label}",
                 fontsize=12, fontweight="bold")

    x      = np.arange(len(event_types_present))
    colors = {"fewer copies": "#C00000", "same count": "#70AD47", "more copies": "#2E75B6"}
    bottom = np.zeros(len(event_types_present))

    for cat in category_order:
        heights = [
            (analyzed_df[analyzed_df["event_type"] == et]["count_category"] == cat).sum()
            for et in event_types_present
        ]
        ax.bar(x, heights, bottom=bottom, label=cat, color=colors[cat], edgecolor="white")
        bottom += np.array(heights)

    ax.set_xticks(x)
    ax.set_xticklabels([display_label(et) for et in event_types_present], fontsize=9, rotation=30, ha="right")
    ax.set_ylabel("Number of rows")
    ax.legend()
    fig.tight_layout()

    save(fig, f"domain_count_change_{label.lower().replace(' ', '_')}.png")


# ══════════════════════════════════════════════════════════════════════════════
# 8. DOMAIN SEVERITY SPECTRUM (single file)
# ══════════════════════════════════════════════════════════════════════════════

def _severity_values(df):
    """
    Vectorised computation of, for every domain-impacting row, the %% of the
    canonical domain's length retained in the alternative transcript:
      0%   = domain completely dropped
      100% = domain unchanged
      >100% = domain got longer
    Shorter/longer/same rows require c_domain_length > 0. Split-domain rows
    are excluded (the concept of % retained is ambiguous there).

    Vectorised rather than a row-by-row loop: df.iterrows() would force
    pandas to materialise the entire frame (all columns, incl. our compact
    category dtypes) as one big object-dtype array via `.values`, which is
    what was blowing up memory on a multi-million-row IOE file.
    """
    et = df["event_type"]
    dropped_vals = np.zeros(int((et == "dropped domain").sum()))

    c = df.get("c_domain_length")
    t = df.get("t_domain_length")
    if c is not None and t is not None:
        valid_mask = et.isin(["shorter", "longer", "same"]) & c.notna() & t.notna() & (c > 0)
        other_vals = (100 * t[valid_mask] / c[valid_mask]).to_numpy(dtype=float)
    else:
        other_vals = np.array([])

    return np.concatenate([dropped_vals, other_vals])


def severity_spectrum(df, label):
    """Print severity-spectrum stats and plot a single-file histogram."""
    vals = _severity_values(df)

    print(f"\n{'='*70}")
    print(f"DOMAIN SEVERITY SPECTRUM — {label}  (%% of canonical domain length retained)")
    print(f"{'='*70}")
    print("  0% = completely dropped  |  100% = unchanged  |  >100% = longer")
    if vals.size == 0:
        print("  No domain-length data available.")
        return

    print(f"  n={len(vals)}")
    print(f"    mean: {vals.mean():.1f}%   median: {np.median(vals):.1f}%   "
          f"range: [{vals.min():.0f}%, {vals.max():.0f}%]")
    print(f"    completely dropped (0%): {(vals == 0).sum()}  "
          f"| unchanged (100%): {(vals == 100).sum()}  "
          f"| longer (>100%): {(vals > 100).sum()}")

    fig, ax = plt.subplots(figsize=(9, 4))
    fig.suptitle(f"Domain Severity Spectrum — {label}", fontsize=12, fontweight="bold")

    bins = np.linspace(0, max(vals.max(), 110), 25)
    ax.hist(vals, bins=bins, alpha=0.8, color=COLORS.get(label, "#5B9BD5"), edgecolor="white")
    ax.axvline(100, color="black", linewidth=0.8, linestyle="--", label="Unchanged (100%)")
    ax.axvline(0,   color="#C00000", linewidth=0.8, linestyle="--", label="Dropped (0%)")
    ax.set_xlabel("% of canonical domain length retained in alternative transcript")
    ax.set_ylabel("Count")
    ax.legend()

    save(fig, f"severity_spectrum_{label.lower().replace(' ', '_')}.png")


# ══════════════════════════════════════════════════════════════════════════════
# Per-file entry point
# ══════════════════════════════════════════════════════════════════════════════

def analyze_file(path, label=None, fetch_domain_descriptions=False, specie_filter=None):
    """
    Load one results.csv and run every single-file analysis on it, printing
    a report and saving plots to OUT_DIR. Missing optional columns (e.g. no
    'specie' column) produce a warning and the analyses that need them are
    skipped, rather than failing.

    fetch_domain_descriptions=True looks up the top domains on InterPro and
    writes domain_descriptions_<label>.csv (see fetch_domain_descriptions) -
    off by default since it hits the network.

    specie_filter, e.g. "human"/"mouse", restricts to just that specie's
    rows before running any analysis, and the specie name is appended to
    `label` (so outputs get their own "..._<label>_human.png"-style
    filenames instead of overwriting the unfiltered run's). Requires a
    'specie' column - raises ValueError if the file doesn't have one.

    Analyses that characterise the distribution of analyzed outcomes
    (domain_frequency, length_change_shorter/longer, domain_count_change,
    severity_spectrum, splice_type_vs_outcome) run on the output of
    select_representative_transcript() - one transcript per cluster - so a
    cluster with several comparable transcripts doesn't get weighted by how
    many alternatives happened to be compared. Ambiguous clusters (see
    select_representative_transcript's docstring) are dropped from that view
    with a warning rather than failing the whole report. event_distribution,
    mixed_combinations and species_comparison instead reduce to one result
    per cluster themselves (see _cluster_event_labels) - unanalyzable
    clusters from all their rows, analyzable ones from the representative
    transcript - so "one row per cluster" event types (e.g.
    no_canonical_transcript) aren't mixed with "one row per candidate
    transcript" ones (e.g. no_unique_junctions) in the same tally.

    Returns the loaded, normalised DataFrame, so callers can pass several of
    these into compare_files() without reloading.
    """
    df, label = _load_and_prepare(path, label, specie_filter)
    _run_analyses(df, label, fetch_domain_descriptions=fetch_domain_descriptions)
    return df


def _load_and_prepare(path, label, specie_filter=None):
    """Load + (optionally) species-filter one results.csv - the first half of analyze_file(), split out so _run_pdf_report() can load every label's df up front, before running any analyses."""
    label = label or os.path.splitext(os.path.basename(path))[0]
    df = _load_single(path, label)
    if specie_filter is not None:
        if "specie" not in df.columns:
            raise ValueError(f"specie_filter={specie_filter!r} given but {label} has no 'specie' column.")
        df = df[df["specie"] == specie_filter].copy()
        label = f"{label} {specie_filter.capitalize()}"
    return df, label


def _run_analyses(df, label, fetch_domain_descriptions=False, skip_sections=None):
    """
    Run every single-file analysis on an already-loaded df - the second
    half of analyze_file(). skip_sections (a set of names: "event_distribution",
    "length_change_shorter", "length_change_longer") lets a caller that's
    already rendering those separately (see _run_pdf_report's combined
    multi-species pages) skip the redundant per-label copy.
    """
    skip_sections = skip_sections or ()

    print("\n" + "#" * 70)
    print(f"# ANALYZING {label}")
    print("#" * 70)

    representative_df = select_representative_transcript(df, on_ambiguous="drop")
    n_clusters = len(representative_df.groupby(
        ["specie", "cluster"] if "specie" in representative_df.columns else ["cluster"], observed=True
    ))
    print(f"  Representative-transcript selection: {len(representative_df)} analyzed rows "
          f"across {n_clusters} clusters kept for domain/length/count/severity analyses below.")

    if "event_distribution" not in skip_sections:
        event_distribution(df, label)
    mixed_combinations(df, label)  # unanalyzable clusters, then analyzable clusters
    top_domains = domain_frequency(representative_df, label)  # top-20-domains graph
    domain_species_rank_scatter(representative_df, label, top_domains)  # rank/count scatter, all + top-20
    if top_domains is not None:
        domain_cluster_species_counts(representative_df, top_domains, label)  # ...cluster count by species
    if "length_change_longer" not in skip_sections:
        length_change_longer(representative_df, label)
    if "length_change_shorter" not in skip_sections:
        length_change_shorter(representative_df, label)
    species_comparison(df, label)  # per-gene analysis
    if fetch_domain_descriptions and top_domains is not None:  # top-domains table
        _fetch_and_save_domain_descriptions(top_domains, label)
    splice_type_vs_outcome(representative_df, label)
    domain_count_change(representative_df, label)
    severity_spectrum(representative_df, label)


# ══════════════════════════════════════════════════════════════════════════════
# Cross-file comparison
# ══════════════════════════════════════════════════════════════════════════════

def event_type_comparison(named_dfs):
    """
    Side-by-side comparison of event-type distributions across two or more
    already-loaded files (as returned by analyze_file()). Shows both raw
    counts and percentages so differences in dataset size don't obscure the
    comparison.

    Uses one result per cluster, not per row (see _cluster_event_labels) -
    same fix as event_distribution()/species_comparison(): counting raw
    df["event_type"] rows directly would count a cluster's own event type
    once per alternative transcript compared, instead of once for the
    cluster - overweighting clusters with more candidate transcripts.
    "mixed"/"domain swap" are included as their own bars here too, matching
    event_distribution()'s categories (this function previously only
    iterated the single-reason UNANALYZABLE_TYPES/ANALYZED_TYPES lists,
    which would have silently excluded every mixed/domain-swap cluster from
    both the counts and the % denominator).
    """
    cluster_labels_by_name = {name: _cluster_event_labels(df) for name, df in named_dfs.items()}

    print(f"\n{'='*70}")
    print(f"EVENT-TYPE COMPARISON — {' vs '.join(named_dfs)}")
    print(f"{'='*70}")

    for group_label, status, group_types in [
        ("Unanalyzable", "unanalyzable", UNANALYZABLE_TYPES + ["mixed"]),
        ("Analyzed",     "analyzable",   ANALYZED_TYPES + ["domain swap", "mixed"]),
    ]:
        header = "".join(f"{name:>16s}" for name in named_dfs)
        print(f"\n  {group_label}:")
        print(f"    {'Event type':<42s} {header}")
        print(f"    {'-'*(42 + 16*len(named_dfs))}")
        for et in group_types:
            row = ""
            for name in named_dfs:
                status_labels = cluster_labels_by_name[name]
                status_labels = status_labels[status_labels["cluster_status"] == status]
                n_in_group = len(status_labels)
                count = (status_labels["label"] == et).sum()
                pct   = 100 * count / n_in_group if n_in_group else 0
                row += f"{f'{count} ({pct:.1f}%)':>16s}"
            print(f"    {et:<42s} {row}")

    # ── plot: grouped bars, one panel per group ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f"Event-Type Distribution Comparison — {' vs '.join(named_dfs)}",
                 fontsize=13, fontweight="bold")

    n = len(named_dfs)
    w = 0.8 / n
    for ax, status, group_types, title in [
        (axes[0], "unanalyzable", UNANALYZABLE_TYPES + ["mixed"], "Unanalyzable events"),
        (axes[1], "analyzable",   ANALYZED_TYPES + ["domain swap", "mixed"], "Analyzed outcomes"),
    ]:
        x = np.arange(len(group_types))

        for i, name in enumerate(named_dfs):
            status_labels = cluster_labels_by_name[name]
            status_labels = status_labels[status_labels["cluster_status"] == status]
            n_in_group = len(status_labels)
            pcts = [
                100 * (status_labels["label"] == et).sum() / n_in_group if n_in_group else 0
                for et in group_types
            ]
            ax.bar(x + i * w, pcts, w, label=name, color=COLORS.get(name, "#999"), alpha=0.85)

        ax.set_xticks(x + w * (n - 1) / 2)
        ax.set_xticklabels([display_label(et) for et in group_types], fontsize=8, rotation=30, ha="right")
        ax.set_ylabel("% within group")
        ax.set_title(title)
        ax.legend()
        ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))

    save(fig, "event_distribution_comparison.png")


def severity_spectrum_comparison(named_dfs):
    """Overlay each file's severity-spectrum histogram for a direct comparison.

    Uses one transcript per cluster (see select_representative_transcript),
    matching the single-file severity_spectrum() plot in analyze_file().
    """
    values_by_name = {
        name: _severity_values(select_representative_transcript(df, on_ambiguous="drop"))
        for name, df in named_dfs.items()
    }
    values_by_name = {name: vals for name, vals in values_by_name.items() if vals.size > 0}
    if not values_by_name:
        _warn("skipping severity spectrum comparison — no domain-length data in any file.")
        return

    fig, ax = plt.subplots(figsize=(10, 4))
    fig.suptitle(f"Domain Severity Spectrum Comparison — {' vs '.join(values_by_name)}",
                 fontsize=12, fontweight="bold")

    overall_max = max(vals.max() for vals in values_by_name.values())
    bins = np.linspace(0, max(overall_max, 110), 25)

    for name, vals in values_by_name.items():
        ax.hist(vals, bins=bins, alpha=0.6, label=name, color=COLORS.get(name, "#999"), edgecolor="white")

    ax.axvline(100, color="black", linewidth=0.8, linestyle="--", label="Unchanged (100%)")
    ax.axvline(0,   color="#C00000", linewidth=0.8, linestyle="--", label="Dropped (0%)")
    ax.set_xlabel("% of canonical domain length retained in alternative transcript")
    ax.set_ylabel("Count")
    ax.legend()

    save(fig, "severity_spectrum_comparison.png")


def compare_files(named_dfs):
    """
    Run every cross-file analysis over 2+ already-loaded DataFrames (as
    returned by analyze_file()), e.g. compare_files({"IOE": ioe_df, "Hadas":
    hadas_df}).
    """
    if len(named_dfs) < 2:
        _warn("compare_files() needs at least two files - skipping comparison.")
        return

    print("\n" + "#" * 70)
    print(f"# COMPARING {' vs '.join(named_dfs)}")
    print("#" * 70)

    event_type_comparison(named_dfs)
    severity_spectrum_comparison(named_dfs)


# ══════════════════════════════════════════════════════════════════════════════
# PDF bundling - every PNG/CSV analyze_file() wrote for a label, in one file
# ══════════════════════════════════════════════════════════════════════════════

def _label_slug(label):
    return label.lower().replace(" ", "_")


def _pdf_title_page(pdf, text, fontsize=20):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.5, 0.5, text, ha="center", va="center", fontsize=fontsize, fontweight="bold", wrap=True)
    pdf.savefig(fig)
    plt.close(fig)


def _crop_suptitle_band(img, search_frac=0.15, min_gap=10):
    """
    Crop off a figure-level suptitle band from the top of `img`, if there is
    one - detected as a blank (all-white) horizontal gap of at least
    min_gap rows, within the top `search_frac` of the image, that follows
    some real content above it (so we don't just match the ordinary margin
    above the first line of text). Only meant to be used right before
    splitting an image into several tiles (see _pdf_image_pages): a
    fig.suptitle() is centered on the *whole original* figure, so slicing
    the image cuts that text in half, and each half then lands off-centre
    on its own page. Per-axes subplot titles sit below the suptitle and are
    left untouched, so they still render intact on every tile.
    """
    h = img.shape[0]
    row_all_white = (img[..., :3].min(axis=2) > 0.97).all(axis=1)
    limit = max(1, int(h * search_frac))

    runs, cur_start = [], None
    for i in range(limit):
        if row_all_white[i]:
            cur_start = i if cur_start is None else cur_start
        elif cur_start is not None:
            runs.append((cur_start, i))
            cur_start = None
    if cur_start is not None:
        runs.append((cur_start, limit))

    for start, end in runs:
        if start > 0 and (end - start) >= min_gap:
            return img[end:]
    return img


def _pdf_image_pages_from_array(pdf, img, name, max_aspect=2.6):
    """
    One page per tile of `img` (a pre-rendered chart, as a numpy array),
    landscape, split into consecutive vertical strips if its aspect ratio
    is too wide for one page - e.g. species_per_gene's one-panel-per-gene
    layout would go illegible squeezed onto a single page otherwise.
    """
    h, w = img.shape[0], img.shape[1]
    aspect = w / h

    n_tiles = max(1, math.ceil(aspect / max_aspect))
    if n_tiles > 1:
        img = _crop_suptitle_band(img)
        h = img.shape[0]
    tile_w = math.ceil(w / n_tiles)
    for i in range(n_tiles):
        x0, x1 = i * tile_w, min(w, (i + 1) * tile_w)
        tile = img[:, x0:x1]
        title = name if n_tiles == 1 else f"{name}  (part {i + 1}/{n_tiles})"
        fig = plt.figure(figsize=(11, 8.5))
        ax = fig.add_axes([0.02, 0.02, 0.96, 0.92])
        ax.imshow(tile)
        ax.axis("off")
        fig.text(0.02, 0.97, title, fontsize=8, color="#555555")
        pdf.savefig(fig, dpi=200)
        plt.close(fig)


def _savefig_to_pdf(pdf, fig, name, max_aspect=2.6):
    """
    Rasterize `fig` in memory (no PNG file written to disk) and add it to
    `pdf` via _pdf_image_pages_from_array - this is what save() uses while
    a PDF report is being built.
    """
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    img = mpimg.imread(buf)
    _pdf_image_pages_from_array(pdf, img, name, max_aspect=max_aspect)




def _pdf_event_distribution_combined(pdf, label_dfs, page_title, mixed_threshold_pct=5):
    """
    One page: a 2 x len(label_dfs) grid of event-type pies (row 1 =
    unanalyzable clusters, row 2 = analyzable clusters; one column per
    label) - a GridSpec keeps every pie the same size, so columns line up
    evenly and both rows stay vertically aligned across labels.

    "mixed" is broken into its individual combinations wherever a
    combination is itself >= mixed_threshold_pct of that panel's total (see
    _mixed_combo_breakdown) - a single "mixed" wedge covering 90%+ of the
    unanalyzable pie isn't useful, even though that's common in real data.
    Anything smaller stays folded into a residual "other mixed" wedge.

    Colors for these combo wedges are assigned by combo *identity*, computed
    once across all labels before drawing (pass 1), so the same combination
    gets the same color in every label's column - unlike event_distribution()'s
    own per-file pie, which shades by rank within that one pie (fine there
    since there's only one panel; here the same shade would otherwise mean a
    different combo in each column).
    """
    n = len(label_dfs)
    fig = plt.figure(figsize=(4.3 * n, 9))
    gs = fig.add_gridspec(4, n, height_ratios=[3, 0.6, 3, 0.6], hspace=0.4)
    fig.suptitle(page_title, fontsize=14, fontweight="bold")

    row_specs = [
        (0, 1, "unanalyzable", UNANALYZABLE_TYPES, "Unanalyzable clusters"),
        (2, 3, "analyzable",   ANALYZED_TYPES + ["domain swap"], "Analyzable clusters"),
    ]
    for pie_row, legend_row, status, group_types, row_title in row_specs:
        # ── pass 1: each label's combo breakdown, and a stable combo -> color map ──
        per_label = {}
        combo_order = []
        for label, df in label_dfs.items():
            cluster_labels = _cluster_event_labels(df)
            group_df = cluster_labels[cluster_labels["cluster_status"] == status]
            un_sets, an_sets = _mixed_combo_sets(df)
            combo_sets = un_sets if status == "unanalyzable" else an_sets
            items = _mixed_combo_breakdown(combo_sets, len(group_df), mixed_threshold_pct,
                                            max_items=len(MIXED_COMBO_SHADES) - 1)
            per_label[label] = (group_df, items)
            for combo_label, _ in items:
                if not combo_label.startswith("other mixed") and combo_label not in combo_order:
                    combo_order.append(combo_label)

        combo_color = {c: MIXED_COMBO_SHADES[i] for i, c in enumerate(combo_order)}

        # ── pass 2: draw, using the shared combo_color map ──
        legend_colors = {}
        for col_idx, label in enumerate(label_dfs):
            group_df, items = per_label[label]
            ax = fig.add_subplot(gs[pie_row, col_idx])

            labels, sizes, colors = _event_type_pie_data(group_df, group_types)
            for combo_label, count in items:
                labels.append(combo_label)
                sizes.append(count)
                colors.append(combo_color.get(combo_label, MIXED_COMBO_SHADES[-1]))

            _draw_event_type_pie(ax, labels, sizes, colors)
            ax.set_title(f"{label}\n{row_title} (n={len(group_df)})", fontsize=9)
            for lab, col in zip(labels, colors):
                # "other mixed (N combos)" - N is panel-specific, so it
                # doesn't belong as a distinct legend entry per label; one
                # shared "other mixed" row covers all of them.
                legend_key = "other mixed" if lab.startswith("other mixed") else lab
                legend_colors.setdefault(legend_key, col)

        legend_ax = fig.add_subplot(gs[legend_row, :])
        legend_ax.axis("off")
        ordered = (
            [dl for dl in (display_label(et) for et in group_types) if dl in legend_colors]
            + [c for c in combo_order if c in legend_colors]
            + (["other mixed"] if "other mixed" in legend_colors else [])
        )
        handles = [Patch(color=legend_colors[dl]) for dl in ordered]
        if handles:
            # Long "reason + reason" combo labels can make this legend
            # wider than the figure - ncol=1 for anything that long so it
            # wraps into rows tall enough to fit instead of clipping at the
            # page edge (see also event_distribution()'s own per-pie legend).
            ncol = 1 if any(len(o) > 25 for o in ordered) else min(len(ordered), 5)
            legend_ax.legend(handles, ordered, loc="center", ncol=ncol,
                              fontsize=8, frameon=False)

    pdf.savefig(fig, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _wrap_cell(value, width):
    text = "" if pd.isna(value) else str(value)
    return "\n".join(textwrap.wrap(text, width=width)) or text


def _pdf_csv_pages(pdf, csv_path, max_lines_per_page=45, table_char_width=150, max_rows=300):
    """Read `csv_path` and render it via _pdf_csv_pages_from_df - see that function for the actual rendering."""
    df_full = pd.read_csv(csv_path)
    _pdf_csv_pages_from_df(pdf, df_full, os.path.basename(csv_path),
                            max_lines_per_page=max_lines_per_page, table_char_width=table_char_width,
                            max_rows=max_rows)


def _pdf_csv_pages_from_df(pdf, df_full, name, max_lines_per_page=45, table_char_width=150, max_rows=300):
    """
    Render `df_full` as one or more table pages, paginated by cumulative
    wrapped-text line count (not a fixed row count) so a table with long
    text (e.g. InterPro descriptions) doesn't spill off the page while a
    compact one (e.g. mixed_combinations) still fits many rows per page.

    Each column gets its own wrap width, proportional to how wide its
    content naturally wants to be (capped, so one very long free-text
    column - e.g. "description" - can't starve the others down to
    illegibly narrow wrapping). A single fixed wrap width for every column
    was wrapping wide columns far narrower than their actual rendered
    table width, producing far more lines - and PDF pages - than needed,
    with each row mostly whitespace.

    Tables above max_rows (e.g. species_per_gene_comparison, which can have
    tens of thousands of rows - one per gene x label) are truncated to the
    first max_rows with a note - a matplotlib table page per handful of rows
    doesn't scale to that, and the PDF is a browsable summary, not a data
    dump; _save_table() writes the full, untruncated table to a CSV
    alongside the PDF whenever it's bigger than max_rows.
    """
    if df_full.empty:
        return
    truncated = len(df_full) > max_rows
    df = df_full.head(max_rows)

    # Median rather than max - a single unusually long value (e.g. one long
    # "name" among mostly-short ones) shouldn't force that column wide for
    # every row; it just wraps onto an extra line instead. (75th percentile
    # was tried first but often still exceeded the longest actual value in
    # a column, so it never actually forced a wrap - no narrower in practice.)
    natural_widths = []
    for col in df.columns:
        lengths = [len(str(v)) for v in df[col] if pd.notna(v)] or [0]
        natural_widths.append(min(max(len(col), round(np.median(lengths))), 60))
    total_natural = sum(natural_widths) or 1
    wrap_widths = [max(8, round(w / total_natural * table_char_width)) for w in natural_widths]

    wrapped = pd.DataFrame({
        col: df[col].apply(lambda v, ww=ww: _wrap_cell(v, ww))
        for col, ww in zip(df.columns, wrap_widths)
    })
    line_counts = wrapped.map(lambda v: v.count("\n") + 1).max(axis=1)

    col_widths = []
    for col in df.columns:
        max_len = max([len(col)] + [len(line) for v in wrapped[col] for line in v.split("\n")])
        col_widths.append(max_len)
    total_w = sum(col_widths) or 1
    col_widths = [w / total_w for w in col_widths]

    pages, current, current_lines = [], [], 0
    for idx in range(len(df)):
        row_lines = line_counts.iloc[idx]
        if current and current_lines + row_lines > max_lines_per_page:
            pages.append(current)
            current, current_lines = [], 0
        current.append(idx)
        current_lines += row_lines
    if current:
        pages.append(current)

    for page_num, idxs in enumerate(pages, start=1):
        chunk = wrapped.iloc[idxs]
        row_lines = line_counts.iloc[idxs]
        note = f"  (page {page_num}/{len(pages)})" if len(pages) > 1 else ""
        total_note = (f"showing first {max_rows} of {len(df_full)} row(s) - see the CSV for the rest"
                      if truncated else f"{len(df_full)} row(s) total")

        fig, ax = plt.subplots(figsize=(11, 8.5))
        ax.axis("off")
        ax.set_title(f"{name}{note} — {total_note}", fontsize=9, loc="left")

        table = ax.table(cellText=chunk.values, colLabels=df.columns, loc="upper center", cellLoc="left")
        table.auto_set_font_size(False)
        table.set_fontsize(6)

        total_lines = row_lines.sum() + 2  # header
        for c, w in enumerate(col_widths):
            table[(0, c)].set_height(2 / total_lines)
            table[(0, c)].set_width(w)
        for i, rl in enumerate(row_lines):
            for c, w in enumerate(col_widths):
                table[(i + 1, c)].set_height(rl / total_lines)
                table[(i + 1, c)].set_width(w)

        pdf.savefig(fig)
        plt.close(fig)


_COMBINE_PAGE_TITLES = {
    "event_distribution": "Event-Type Distribution",
    "length_change_shorter": "Length Change — Shorter",
    "length_change_longer": "Length Change — Longer",
}


def _run_pdf_report(runs, pdf_path, title, fetch_domain_descriptions=False, combine_sections=None):
    """
    Build one multi-page PDF covering one or more analyze_file() runs back
    to back, e.g. Hadas plus its human-only/mouse-only species split.
    `runs` is a list of (label, path, specie_filter) tuples.

    Unlike the old file-based build_pdf_report(), every chart/table is
    rendered *directly* into the open PdfPages as each analysis runs (see
    save()/_save_table(), which route to the currently-open PDF instead of
    OUT_DIR while it's set as _ACTIVE_PDF) - no intermediate PNG/CSV files
    for any of it, except tables too large for the PDF's row cap (see
    _save_table), which still get a full CSV alongside the PDF.

    combine_sections: names ("event_distribution", "length_change_shorter",
    "length_change_longer") to render once, up front, as a single page
    aligned across every label (see _pdf_event_distribution_combined /
    _pdf_length_change_combined) instead of once per label's own section -
    each label's own per-section analysis is then skipped for these names
    (see _run_analyses' skip_sections) so it isn't drawn twice.

    Returns {label: df}, matching what analyze_file() would return for each
    run - callers use this instead of a return value from analyze_file()
    itself, since analyze_file() isn't called directly here.
    """
    combine_sections = combine_sections or []
    prepared = [(*_load_and_prepare(path, label, specie), ) for label, path, specie in runs]
    label_dfs = {label: df for df, label in prepared}

    global _ACTIVE_PDF, _DEFERRED_TABLES
    with PdfPages(pdf_path) as pdf:
        _pdf_title_page(pdf, title, fontsize=24)
        _ACTIVE_PDF = pdf
        _DEFERRED_TABLES = []
        try:
            joined_labels = " / ".join(label_dfs)
            for prefix in combine_sections:
                page_title = f"{_COMBINE_PAGE_TITLES.get(prefix, prefix)} — {joined_labels}"
                if prefix == "event_distribution":
                    _pdf_event_distribution_combined(pdf, label_dfs, page_title)
                elif prefix in ("length_change_shorter", "length_change_longer"):
                    event_type = "shorter" if prefix == "length_change_shorter" else "longer"
                    _pdf_length_change_combined(pdf, label_dfs, event_type, page_title)

            for df, label in prepared:
                _pdf_title_page(pdf, label, fontsize=20)
                _run_analyses(df, label, fetch_domain_descriptions=fetch_domain_descriptions,
                               skip_sections=combine_sections)

            # Every table collected above (see _save_table/_DEFERRED_TABLES),
            # in the same order it was encountered - all together at the end
            # instead of interleaved into each section.
            if _DEFERRED_TABLES:
                _pdf_title_page(pdf, "Tables", fontsize=24)
                for table_df, table_title in _DEFERRED_TABLES:
                    _pdf_csv_pages_from_df(pdf, table_df, table_title)
        finally:
            _ACTIVE_PDF = None
            _DEFERRED_TABLES = None

    print(f"  → {pdf_path}")
    return label_dfs


# ══════════════════════════════════════════════════════════════════════════════
# Top-level entry point - importable (e.g. from domas.py) as well as CLI
# ══════════════════════════════════════════════════════════════════════════════

def generate_report(ioe_file=None, hadas_file=None, out_dir=None, fetch_domain_descriptions=False,
                     hadas_species_split=False, pdf_report=False):
    """
    Analyze whichever of ioe_file/hadas_file is given, plus compare_files()
    between them if both are given. At least one of ioe_file/hadas_file is
    required.

    fetch_domain_descriptions=True looks up each file's top domains on
    InterPro (see fetch_domain_descriptions) - off by default since it hits
    the network.

    hadas_species_split=True additionally runs the full analysis on Hadas
    restricted to just human, then just mouse (see analyze_file's
    specie_filter) - IOE has no 'specie' column, so this only affects Hadas.

    pdf_report=True (see _run_pdf_report) renders every chart/table straight
    into one PDF per file - Hadas_report.pdf covers "Hadas" plus, if
    hadas_species_split is also set, "Hadas Human"/"Hadas Mouse" as further
    sections of that same PDF; IOE_report.pdf covers just "IOE" - with no
    intermediate PNG/CSV files (aside from a handful of tables too large to
    fit the PDF's row cap, which still get a full CSV alongside it). With
    pdf_report=False, every chart/table is instead written as its own
    scattered PNG/CSV file in OUT_DIR, as analyze_file() always has.

    Returns a {label: df} dict of whatever was analyzed, e.g. {"IOE": df} or
    {"IOE": df, "Hadas": df}. Species-split subsets aren't included in this
    dict (they're not useful for compare_files(), which is IOE-vs-Hadas).
    """
    if not ioe_file and not hadas_file:
        raise ValueError("generate_report() needs at least one of ioe_file/hadas_file.")

    if out_dir:
        global OUT_DIR
        OUT_DIR = out_dir
        os.makedirs(OUT_DIR, exist_ok=True)

    named_dfs = {}
    if ioe_file:
        if pdf_report:
            ioe_dfs = _run_pdf_report([("IOE", ioe_file, None)], os.path.join(OUT_DIR, "IOE_report.pdf"),
                                       title="IOE — Full Report", fetch_domain_descriptions=fetch_domain_descriptions)
            named_dfs["IOE"] = ioe_dfs["IOE"]
        else:
            named_dfs["IOE"] = analyze_file(ioe_file, "IOE", fetch_domain_descriptions=fetch_domain_descriptions)

    if hadas_file:
        if pdf_report:
            # label is "Hadas" for every run here, not pre-suffixed with the
            # specie - _load_and_prepare() appends "Human"/"Mouse" itself
            # when specie_filter is given, so pre-suffixing here would
            # double up into "Hadas Human Human".
            runs = [("Hadas", hadas_file, None)]
            if hadas_species_split:
                runs += [("Hadas", hadas_file, "human"), ("Hadas", hadas_file, "mouse")]
            combine_sections = (
                ["event_distribution", "length_change_longer", "length_change_shorter"]
                if len(runs) > 1 else None
            )
            hadas_dfs = _run_pdf_report(runs, os.path.join(OUT_DIR, "Hadas_report.pdf"),
                                         title="Hadas — Full Report",
                                         fetch_domain_descriptions=fetch_domain_descriptions,
                                         combine_sections=combine_sections)
            named_dfs["Hadas"] = hadas_dfs["Hadas"]
        else:
            named_dfs["Hadas"] = analyze_file(hadas_file, "Hadas", fetch_domain_descriptions=fetch_domain_descriptions)
            if hadas_species_split:
                for specie in ("human", "mouse"):
                    analyze_file(hadas_file, "Hadas", fetch_domain_descriptions=fetch_domain_descriptions,
                                 specie_filter=specie)

    if len(named_dfs) > 1:
        compare_files(named_dfs)

    return named_dfs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ioe_file", default=IOE_FILE, help=f"Path to the IOE results.csv (default: {IOE_FILE})")
    parser.add_argument("--hadas_file", default=HADAS_FILE, help=f"Path to the Hadas results.csv (default: {HADAS_FILE})")
    parser.add_argument("--out_dir", default=OUT_DIR, help=f"Directory to write outputs to (default: {OUT_DIR})")
    parser.add_argument("--fetch_domain_descriptions", action="store_true",
                         help="Also look up each file's top domains on InterPro and write a description table (hits the network).")
    parser.add_argument("--hadas_species_split", action="store_true",
                         help="Also run the full analysis on Hadas restricted to just human, then just mouse.")
    parser.add_argument("--pdf_report", action="store_true",
                         help="Render every chart/table straight into one PDF per file, instead of scattered PNG/CSV files.")
    args = parser.parse_args()

    try:
        generate_report(ioe_file=args.ioe_file, hadas_file=args.hadas_file, out_dir=args.out_dir,
                         fetch_domain_descriptions=args.fetch_domain_descriptions,
                         hadas_species_split=args.hadas_species_split, pdf_report=args.pdf_report)
    except ValueError as e:
        print(f"\nERROR: {e}")
        raise SystemExit(1)

    print("\nDone.")
