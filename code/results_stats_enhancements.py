"""Enhanced multi-source reporting on top of results_stats.py.

Adds, relative to the stock report:
  1. Event-type distribution: pies WITH legends (the stock combined page dropped
     them) plus a grouped bar chart per status, for every source at once.
  2. AS-type x outcome heatmap: event counts per row (and per cell), not just %.
  3. Candidate-transcript profile: how many clusters have 1, 2, 3, ... possible
     transcripts, and for the multi-transcript ones which tie-break rule picked
     the representative (is_most_like_canonical vs is_longest_cds).
  4. Top domains across all sources: one bar per source per domain (grouped, not
     stacked), count + % printed at the end of each bar.
  5. Statistical enrichment/depletion between sources (log2FC, two-proportion
     z-test, BH q-values, chi-square + Cramer's V).
  6. A companion table for every analysis above, carrying the row/cluster counts
     and the statistical outcome.

'Hadas' is displayed as 'Reads' everywhere (the underlying labels are untouched).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
from scipy.stats import chi2_contingency, norm

import results_stats as rs

# ── display ──────────────────────────────────────────────────────────────────
SOURCE_COLORS = {
    "Genome Human": "#5B9BD5",
    "Genome Mouse": "#2E75B6",
    "Reads Human": "#ED7D31",
    "Reads Mouse": "#C55A11",
}
ANALYZED_LABELS = list(rs.ANALYZED_TYPES) + ["domain swap"]
EPS = 0.5  # pseudo-count in %-space so log2FC stays finite for zero-count labels


def disp(label):
    """Display name for a source label - Hadas is presented as Reads."""
    return label.replace("Hadas", "Reads")


def _short(text, width=44):
    text = str(text)
    return text if len(text) <= width else text[: width - 1] + "…"


def _fmt_p(p):
    p = float(p)
    if p == 0:
        return "<1e-300"
    if p >= 0.001:
        return f"{p:.3f}"
    return f"{p:.1e}"


# ── statistics ───────────────────────────────────────────────────────────────
def two_prop_z(c1, n1, c2, n2):
    """Two-proportion z-test. Returns (z, p)."""
    if n1 == 0 or n2 == 0:
        return 0.0, 1.0
    p1, p2 = c1 / n1, c2 / n2
    pool = (c1 + c2) / (n1 + n2)
    se = np.sqrt(pool * (1 - pool) * (1 / n1 + 1 / n2))
    if se == 0:
        return 0.0, 1.0
    z = (p1 - p2) / se
    return float(z), float(2 * norm.sf(abs(z)))


def bh(pvals):
    """Benjamini-Hochberg FDR correction."""
    p = np.asarray(pvals, float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def chi2_and_v(table):
    """Chi-square + Cramer's V for a contingency table. Returns (chi2, p, V)."""
    t = np.asarray(table, float)
    t = t[:, t.sum(axis=0) > 0]
    t = t[t.sum(axis=1) > 0, :]
    if t.shape[0] < 2 or t.shape[1] < 2:
        return np.nan, np.nan, np.nan
    chi2, p, _dof, _exp = chi2_contingency(t)
    n = t.sum()
    k = min(t.shape) - 1
    return float(chi2), float(p), float(np.sqrt(chi2 / (n * k))) if n * k else np.nan


# ── shared per-source preparation ────────────────────────────────────────────
def cluster_label_table(df):
    """One row per cluster: cluster_status + label (see rs._cluster_event_labels)."""
    return rs._cluster_event_labels(df)


def transcript_choice_table(df):
    """
    One row per cluster describing the candidate-transcript situation:
      n_transcripts - how many distinct transcripts were comparable (analyzed rows)
      rule          - which tie-break selected the representative:
                      'most_like_canonical', 'longest_cds',
                      'single candidate', or 'ambiguous (dropped)'
    Mirrors select_representative_transcript()'s priority exactly.
    """
    d = rs.normalize_event_types(df.copy())
    analyzed = d[d["event_type"].isin(rs.ANALYZED_TYPES)]
    if analyzed.empty:
        return pd.DataFrame(columns=["n_transcripts", "rule"])

    gcols = ["specie", "cluster"] if "specie" in analyzed.columns else ["cluster"]
    nuniq = analyzed.groupby(gcols, observed=True, dropna=False)["transcript_id"].nunique()

    # A CSV written without -write_all_comparable carries neither tag column: only
    # the chosen transcript was compared, so every cluster holds exactly one and
    # there is no rule breakdown to report.
    if not {"is_most_like_canonical", "is_longest_cds"} <= set(analyzed.columns):
        return pd.DataFrame({"n_transcripts": nuniq.to_numpy(),
                             "rule": "pre-selected at write time"}, index=nuniq.index)

    # Clusters carrying each tag; .size().index is cheaper than .groups.
    most_like_idx = (analyzed[analyzed["is_most_like_canonical"] == True]
                     .groupby(gcols, observed=True, dropna=False).size().index)
    longest_idx = (analyzed[analyzed["is_longest_cds"] == True]
                   .groupby(gcols, observed=True, dropna=False).size().index)

    idx = nuniq.index
    n_vals = nuniq.to_numpy()
    rule = np.where(
        idx.isin(most_like_idx), "most_like_canonical",
        np.where(idx.isin(longest_idx), "longest_cds",
                 np.where(n_vals == 1, "single candidate", "ambiguous (dropped)")))

    return pd.DataFrame({"n_transcripts": n_vals, "rule": rule})


# ══════════════════════════════════════════════════════════════════════════════
# 1. EVENT-TYPE DISTRIBUTION - pies (with legends) + grouped bars + table
# ══════════════════════════════════════════════════════════════════════════════
def _status_slice(clabels, status):
    return clabels[clabels["cluster_status"] == status]


def event_distribution_pies_page(pdf, label_dfs, clabels_by_label, status, page_title):
    """One page: a pie per source for `status` clusters, each WITH its legend."""
    n = len(label_dfs)
    fig, axes = plt.subplots(1, n, figsize=(5.5 * n, 10))
    if n == 1:
        axes = [axes]
    group_types = rs.UNANALYZABLE_TYPES if status == "unanalyzable" else ANALYZED_LABELS

    for ax, (label, df) in zip(axes, label_dfs.items()):
        group_df = _status_slice(clabels_by_label[label], status)
        counts = group_df["label"].value_counts()
        names = [et for et in group_types if et in counts.index]
        sizes = [int(counts[et]) for et in names]
        colors = [rs.event_color(et) for et in names]

        un_sets, an_sets = rs._mixed_combo_sets(df)
        combo_sets = un_sets if status == "unanalyzable" else an_sets
        m_labels, m_sizes, m_colors = rs._mixed_wedges(combo_sets, len(group_df))

        all_labels = [rs.display_label(et) for et in names] + list(m_labels)
        all_sizes = sizes + list(m_sizes)
        all_colors = colors + list(m_colors)

        wedges = rs._draw_event_type_pie(ax, all_labels, all_sizes, all_colors, autopct_threshold=4)
        total = sum(all_sizes)
        ax.set_title(f"{disp(label)}\n(n={total:,} clusters)", fontsize=11, fontweight="bold")
        if wedges:
            ax.legend(
                wedges, [_short(l) for l in all_labels],
                loc="upper center", bbox_to_anchor=(0.5, -0.02),
                fontsize=7, ncol=1, frameon=False,
            )

    fig.suptitle(page_title, fontsize=14, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    rs._savefig_to_pdf(pdf, fig, page_title)


def event_distribution_bars_page(pdf, label_dfs, clabels_by_label, status, page_title):
    """One page: grouped horizontal bars (one bar per source) per outcome type."""
    group_types = rs.UNANALYZABLE_TYPES if status == "unanalyzable" else ANALYZED_LABELS
    labels = list(label_dfs)

    counts = {}
    totals = {}
    for label in labels:
        gdf = _status_slice(clabels_by_label[label], status)
        vc = gdf["label"].value_counts()
        counts[label] = vc
        totals[label] = len(gdf)

    present = [et for et in list(group_types) + ["mixed"]
               if any(counts[l].get(et, 0) > 0 for l in labels)]
    if not present:
        return None

    n_src = len(labels)
    y = np.arange(len(present))
    height = 0.8 / n_src
    fig, ax = plt.subplots(figsize=(13, max(4.5, len(present) * 0.62 * n_src * 0.42)))

    for i, label in enumerate(labels):
        vals = [100 * counts[label].get(et, 0) / totals[label] if totals[label] else 0
                for et in present]
        raw = [int(counts[label].get(et, 0)) for et in present]
        offs = y + (i - (n_src - 1) / 2) * height
        ax.barh(offs, vals, height=height, color=SOURCE_COLORS.get(disp(label), "#999"),
                edgecolor="white", linewidth=0.4, label=disp(label))
        for yy, v, c in zip(offs, vals, raw):
            if v > 0:
                ax.text(v + max(vals) * 0.012, yy, f"{c:,} ({v:.1f}%)",
                        va="center", fontsize=6.5)

    ax.set_yticks(y)
    ax.set_yticklabels([rs.display_label(et) for et in present], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("% of that source's clusters in this status", fontsize=10)
    ax.set_xlim(0, max(1.0, ax.get_xlim()[1] * 1.18))
    ax.legend(fontsize=8.5, loc="lower right")
    ax.grid(axis="x", alpha=0.2, linestyle="--")
    ax.set_axisbelow(True)
    fig.suptitle(page_title, fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    rs._savefig_to_pdf(pdf, fig, page_title)

    rows = []
    for et in present:
        row = {"outcome": rs.display_label(et)}
        for label in labels:
            c = int(counts[label].get(et, 0))
            row[f"{disp(label)} n"] = c
            row[f"{disp(label)} %"] = round(100 * c / totals[label], 2) if totals[label] else 0.0
        rows.append(row)
    table = pd.DataFrame(rows)
    total_row = {"outcome": "TOTAL clusters"}
    for label in labels:
        total_row[f"{disp(label)} n"] = totals[label]
        total_row[f"{disp(label)} %"] = 100.0
    table = pd.concat([table, pd.DataFrame([total_row])], ignore_index=True)

    contingency = [[int(counts[l].get(et, 0)) for et in present] for l in labels]
    chi2, p, v = chi2_and_v(contingency)
    stat = pd.DataFrame([{
        "test": f"Chi-square, {status} outcome distribution x {n_src} sources",
        "chi2": round(chi2, 2) if np.isfinite(chi2) else "n/a",
        "p": _fmt_p(p) if np.isfinite(p) else "n/a",
        "Cramers_V": round(v, 4) if np.isfinite(v) else "n/a",
        "n_total_clusters": int(sum(totals.values())),
    }])
    rs._save_table(table, f"event_distribution_{status}_counts.csv",
                   title=f"Event distribution — {status} — counts & % per source")
    rs._save_table(stat, f"event_distribution_{status}_test.csv",
                   title=f"Event distribution — {status} — statistical outcome")
    return table


# ══════════════════════════════════════════════════════════════════════════════
# 2. CANDIDATE-TRANSCRIPT PROFILE
# ══════════════════════════════════════════════════════════════════════════════
RULE_COLORS = {
    "most_like_canonical": "#2E8B57",
    "longest_cds": "#7030A0",
    "single candidate": "#B0B7C0",
    "ambiguous (dropped)": "#C00000",
}
RULE_ORDER = ["single candidate", "most_like_canonical", "longest_cds", "ambiguous (dropped)"]
MAX_BUCKET = 8  # clusters with more candidates fold into "8+"


def _bucket(n):
    return str(int(n)) if n < MAX_BUCKET else f"{MAX_BUCKET}+"


def transcript_choice_page(pdf, label_dfs, page_title):
    """
    One page: per source, how many clusters have N possible transcripts, with the
    multi-transcript bars split by which tie-break rule chose the representative.
    """
    per_source = {label: transcript_choice_table(df) for label, df in label_dfs.items()}
    buckets = [str(i) for i in range(1, MAX_BUCKET)] + [f"{MAX_BUCKET}+"]

    n_src = len(label_dfs)
    ncols = 2
    nrows = int(np.ceil(n_src / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4.6 * nrows), squeeze=False)

    table_rows = []
    for idx, (label, tct) in enumerate(per_source.items()):
        ax = axes[idx // ncols][idx % ncols]
        if tct.empty:
            ax.set_visible(False)
            continue
        tct = tct.assign(bucket=tct["n_transcripts"].map(_bucket))
        pivot = (tct.groupby(["bucket", "rule"], observed=True).size()
                 .unstack(fill_value=0).reindex(index=buckets, fill_value=0))
        for r in RULE_ORDER:
            if r not in pivot.columns:
                pivot[r] = 0
        pivot = pivot[RULE_ORDER]

        bottom = np.zeros(len(buckets))
        x = np.arange(len(buckets))
        for rule in RULE_ORDER:
            vals = pivot[rule].to_numpy(dtype=float)
            if vals.sum() == 0:
                continue
            ax.bar(x, vals, bottom=bottom, color=RULE_COLORS[rule], label=rule,
                   edgecolor="white", linewidth=0.4)
            bottom += vals

        for xi, tot in zip(x, bottom):
            if tot > 0:
                ax.text(xi, tot, f"{int(tot):,}", ha="center", va="bottom", fontsize=6.5)

        ax.set_xticks(x)
        ax.set_xticklabels(buckets, fontsize=9)
        ax.set_xlabel("Possible (comparable) transcripts per cluster", fontsize=9)
        ax.set_ylabel("Clusters", fontsize=9)
        ax.set_title(f"{disp(label)}  (n={len(tct):,} clusters)", fontsize=10, fontweight="bold")
        ax.grid(axis="y", alpha=0.2, linestyle="--")
        ax.set_axisbelow(True)

        multi = tct[tct["n_transcripts"] > 1]
        rule_counts = multi["rule"].value_counts()
        n_multi = len(multi)
        table_rows.append({
            "source": disp(label),
            "clusters_total": len(tct),
            "clusters_1_transcript": int((tct["n_transcripts"] == 1).sum()),
            "clusters_multi_transcript": n_multi,
            "multi_%_of_total": round(100 * n_multi / len(tct), 2) if len(tct) else 0.0,
            "chosen_most_like_canonical": int(rule_counts.get("most_like_canonical", 0)),
            "chosen_longest_cds": int(rule_counts.get("longest_cds", 0)),
            "ambiguous_dropped": int(rule_counts.get("ambiguous (dropped)", 0)),
            "most_like_%_of_multi": round(100 * rule_counts.get("most_like_canonical", 0) / n_multi, 2) if n_multi else 0.0,
            "longest_cds_%_of_multi": round(100 * rule_counts.get("longest_cds", 0) / n_multi, 2) if n_multi else 0.0,
            "mean_transcripts": round(float(tct["n_transcripts"].mean()), 3),
            "max_transcripts": int(tct["n_transcripts"].max()),
        })

    for j in range(n_src, nrows * ncols):
        axes[j // ncols][j % ncols].set_visible(False)

    handles = [Patch(facecolor=RULE_COLORS[r], label=r) for r in RULE_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=len(RULE_ORDER), fontsize=9, frameon=False)
    fig.suptitle(page_title, fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    rs._savefig_to_pdf(pdf, fig, page_title)

    tbl = pd.DataFrame(table_rows)
    rs._save_table(tbl, "transcript_choice_summary.csv",
                   title="Candidate transcripts per cluster & tie-break rule used")

    # enrichment of the tie-break rule across sources (multi-transcript clusters only)
    labels = [r["source"] for r in table_rows]
    cont = [[r["chosen_most_like_canonical"], r["chosen_longest_cds"], r["ambiguous_dropped"]]
            for r in table_rows]
    chi2, p, v = chi2_and_v(cont)
    rs._save_table(pd.DataFrame([{
        "test": "Chi-square, tie-break rule (most_like/longest_cds/ambiguous) x sources",
        "sources": ", ".join(labels),
        "chi2": round(chi2, 2) if np.isfinite(chi2) else "n/a",
        "p": _fmt_p(p) if np.isfinite(p) else "n/a",
        "Cramers_V": round(v, 4) if np.isfinite(v) else "n/a",
        "n_multi_transcript_clusters": int(np.sum(cont)),
    }]), "transcript_choice_test.csv",
        title="Candidate transcripts — statistical outcome")
    return tbl


# ══════════════════════════════════════════════════════════════════════════════
# 3. AS-TYPE x OUTCOME HEATMAP WITH COUNTS
# ══════════════════════════════════════════════════════════════════════════════
def splice_type_vs_outcome_with_counts(df, label):
    """
    Replacement for rs.splice_type_vs_outcome: each row is annotated with the
    number of events on that line, and every cell shows % and the raw count.
    """
    if "as_type" not in df.columns:
        rs._warn(f"skipping AS splice type vs outcome for {label} — no as_type column.")
        return
    as_df = df[df["as_type"] != ""]
    if as_df.empty:
        rs._warn(f"skipping AS splice type vs outcome for {disp(label)} — no AS-type information in cluster ids.")
        return

    as_types = sorted(as_df["as_type"].unique())
    outcome_types = [et for et in rs.ANALYZED_TYPES if et in as_df["event_type"].values]

    counts = np.zeros((len(as_types), len(outcome_types)), dtype=int)
    row_totals = np.zeros(len(as_types), dtype=int)
    for i, at in enumerate(as_types):
        sub = as_df[as_df["as_type"] == at]
        row_totals[i] = len(sub)
        for j, et in enumerate(outcome_types):
            counts[i, j] = int((sub["event_type"] == et).sum())

    pct = np.where(row_totals[:, None] > 0, 100 * counts / np.maximum(row_totals[:, None], 1), 0.0)

    print(f"\n{'='*70}")
    print(f"AS SPLICE TYPE vs DOMAIN OUTCOME — {disp(label)}")
    print(f"{'='*70}")
    for i, at in enumerate(as_types):
        detail = ", ".join(f"{et}: {counts[i, j]}" for j, et in enumerate(outcome_types))
        print(f"  {at:<6s}  n={row_totals[i]:<8d} [{detail}]")

    # Keep the aspect under _savefig_to_pdf's 2.6 limit; past it the figure is
    # sliced across pages and the per-cell counts are cut in half.
    _w = 1.15 * len(outcome_types) + 4.0
    fig, ax = plt.subplots(figsize=(_w, max(_w / 2.3, len(as_types) * 0.95)))
    fig.suptitle(f"AS Splice Type vs Domain Outcome — {disp(label)}\n"
                 f"cell = % of that AS type's events (count below)",
                 fontsize=12, fontweight="bold")
    im = ax.imshow(pct, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(outcome_types)))
    ax.set_xticklabels([rs.display_label(et) for et in outcome_types], rotation=25, ha="right", fontsize=9)
    ax.set_yticks(range(len(as_types)))
    # Event count per line.
    ax.set_yticklabels([f"{at}  (n={row_totals[i]:,})" for i, at in enumerate(as_types)], fontsize=9)

    vmax = pct.max() if pct.size else 0
    for i in range(len(as_types)):
        for j in range(len(outcome_types)):
            if counts[i, j] == 0:
                continue
            color = "white" if pct[i, j] > 0.6 * vmax else "black"
            ax.text(j, i - 0.13, f"{pct[i, j]:.0f}%", ha="center", va="center",
                    fontsize=8, color=color, fontweight="bold")
            ax.text(j, i + 0.22, f"n={counts[i, j]:,}", ha="center", va="center",
                    fontsize=6.5, color=color)

    plt.colorbar(im, ax=ax, label="% of that AS type's events")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    rs.save(fig, f"splice_type_vs_outcome_{label.lower().replace(' ', '_')}.png")

    tbl = pd.DataFrame(counts, columns=[rs.display_label(et) for et in outcome_types])
    tbl.insert(0, "events_on_line", row_totals)
    tbl.insert(0, "AS_type", as_types)
    chi2, p, v = chi2_and_v(counts)
    tbl_stat = pd.DataFrame([{
        "test": "Chi-square, AS type x domain outcome",
        "source": disp(label),
        "chi2": round(chi2, 2) if np.isfinite(chi2) else "n/a",
        "p": _fmt_p(p) if np.isfinite(p) else "n/a",
        "Cramers_V": round(v, 4) if np.isfinite(v) else "n/a",
        "n_events": int(counts.sum()),
    }])
    rs._save_table(tbl, f"splice_type_counts_{label.lower().replace(' ', '_')}.csv",
                   title=f"AS type x outcome — event counts — {disp(label)}")
    rs._save_table(tbl_stat, f"splice_type_test_{label.lower().replace(' ', '_')}.csv",
                   title=f"AS type x outcome — statistical outcome — {disp(label)}")


# ══════════════════════════════════════════════════════════════════════════════
# 4. TOP DOMAINS ACROSS ALL SOURCES (grouped bars, count + % at bar end)
# ══════════════════════════════════════════════════════════════════════════════
def cross_source_domain_table(label_dfs, top_n=20):
    """Counts and % per domain per source, over representative-transcript rows."""
    counts, totals = {}, {}
    for label, df in label_dfs.items():
        rep = rs.select_representative_transcript(df, on_ambiguous="drop")
        rep = rep[rep["event_type"].isin(rs.ANALYZED_TYPES)]
        rep = rep[rep["domain_name"].notna()]
        counts[disp(label)] = rep["domain_name"].value_counts()
        totals[disp(label)] = int(len(rep))

    all_domains = set()
    for vc in counts.values():
        all_domains.update(vc.index)

    rows = []
    for dom in all_domains:
        row = {"domain": dom}
        tot = 0
        for src in counts:
            c = int(counts[src].get(dom, 0))
            row[f"{src} n"] = c
            row[f"{src} %"] = round(100 * c / totals[src], 4) if totals[src] else 0.0
            tot += c
        row["total_n"] = tot
        rows.append(row)

    df_out = pd.DataFrame(rows).sort_values("total_n", ascending=False).reset_index(drop=True)
    return df_out.head(top_n), totals


def cross_source_domain_pages(pdf, label_dfs, top_n=20, per_page=10):
    """
    Grouped horizontal bars: one bar per source for each top domain, with the
    count and percentage printed at the end of every bar. Paginated so the bars
    stay legible.
    """
    table, totals = cross_source_domain_table(label_dfs, top_n=top_n)
    sources = [disp(l) for l in label_dfs]
    n_src = len(sources)

    for start in range(0, len(table), per_page):
        chunk = table.iloc[start:start + per_page]
        y = np.arange(len(chunk))
        height = 0.8 / n_src
        # A minimum height keeps the aspect under _savefig_to_pdf's 2.6 limit.
        fig, ax = plt.subplots(figsize=(13, max(6.0, len(chunk) * n_src * 0.32 + 1.6)))

        xmax = max(
            (chunk[f"{s} %"].max() for s in sources if f"{s} %" in chunk), default=1.0
        ) or 1.0
        for i, src in enumerate(sources):
            pcts = chunk[f"{src} %"].to_numpy(dtype=float)
            raw = chunk[f"{src} n"].to_numpy(dtype=int)
            offs = y + (i - (n_src - 1) / 2) * height
            ax.barh(offs, pcts, height=height, color=SOURCE_COLORS.get(src, "#999"),
                    edgecolor="white", linewidth=0.4, label=src)
            for yy, p, c in zip(offs, pcts, raw):
                ax.text(p + xmax * 0.015, yy, f"{c:,} ({p:.2f}%)",
                        va="center", fontsize=6.8)

        ax.set_yticks(y)
        ax.set_yticklabels(chunk["domain"].tolist(), fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("% of that source's domain-impacting rows", fontsize=10)
        ax.set_xlim(0, xmax * 1.30)
        ax.legend(fontsize=8.5, loc="lower right")
        ax.grid(axis="x", alpha=0.2, linestyle="--")
        ax.set_axisbelow(True)
        page = start // per_page + 1
        n_pages = int(np.ceil(len(table) / per_page))
        title = f"Top domains across sources ({start + 1}-{min(start + per_page, len(table))} of {len(table)})"
        fig.suptitle(f"{title}   [page {page}/{n_pages}]", fontsize=13, fontweight="bold")
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        rs._savefig_to_pdf(pdf, fig, title)

    tot_row = {"domain": "TOTAL domain-impacting rows"}
    for s in sources:
        tot_row[f"{s} n"] = totals[s]
        tot_row[f"{s} %"] = 100.0
    tot_row["total_n"] = int(sum(totals.values()))
    out = pd.concat([table, pd.DataFrame([tot_row])], ignore_index=True)
    rs._save_table(out, "top_domains_cross_source.csv",
                   title="Top domains across sources — counts & % per source")

    cont = [[int(table[f"{s} n"].iloc[i]) for s in sources] for i in range(len(table))]
    chi2, p, v = chi2_and_v(np.array(cont).T)
    rs._save_table(pd.DataFrame([{
        "test": f"Chi-square, top-{len(table)} domain composition x sources",
        "sources": ", ".join(sources),
        "chi2": round(chi2, 2) if np.isfinite(chi2) else "n/a",
        "p": _fmt_p(p) if np.isfinite(p) else "n/a",
        "Cramers_V": round(v, 4) if np.isfinite(v) else "n/a",
        "n_rows_counted": int(np.sum(cont)),
    }]), "top_domains_test.csv", title="Top domains across sources — statistical outcome")
    return table, totals


# ══════════════════════════════════════════════════════════════════════════════
# 5. STATISTICAL ENRICHMENT / DEPLETION
# ══════════════════════════════════════════════════════════════════════════════
def enrichment_table(a_labels, b_labels, a_name, b_name, universe=ANALYZED_LABELS):
    """
    Enrichment of each outcome in `a` relative to `b`, over per-cluster labels.
    Positive log2FC = enriched in A.
    """
    ca = a_labels["label"].value_counts()
    cb = b_labels["label"].value_counts()
    names = [l for l in universe if (ca.get(l, 0) + cb.get(l, 0)) > 0]
    nA, nB = int(sum(ca.get(l, 0) for l in names)), int(sum(cb.get(l, 0) for l in names))

    rows = []
    for l in names:
        cA, cB = int(ca.get(l, 0)), int(cb.get(l, 0))
        pA = 100 * cA / nA if nA else 0.0
        pB = 100 * cB / nB if nB else 0.0
        z, p = two_prop_z(cA, nA, cB, nB)
        rows.append({
            "outcome": rs.display_label(l),
            f"{a_name} n": cA, f"{a_name} %": round(pA, 3),
            f"{b_name} n": cB, f"{b_name} %": round(pB, 3),
            "log2FC": round(float(np.log2((pA + EPS) / (pB + EPS))), 3),
            "direction": "enriched" if pA > pB else ("depleted" if pA < pB else "equal"),
            "z": round(z, 3), "p_raw": p,
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return out, (np.nan, np.nan, np.nan), (0, 0)
    out["q_BH"] = bh(out["p_raw"].to_numpy())
    out["significant_q<0.05"] = np.where(out["q_BH"] < 0.05, "yes", "no")
    out["p"] = out["p_raw"].map(_fmt_p)
    out["q"] = out["q_BH"].map(_fmt_p)
    out = out.drop(columns=["p_raw", "q_BH"])
    out = out.reindex(out["log2FC"].abs().sort_values(ascending=False).index).reset_index(drop=True)
    stats = chi2_and_v([[int(ca.get(l, 0)) for l in names], [int(cb.get(l, 0)) for l in names]])
    return out, stats, (nA, nB)


def enrichment_page(pdf, a_labels, b_labels, a_name, b_name, page_title):
    """Diverging log2FC bars + the full statistics table for one comparison."""
    tbl, (chi2, p, v), (nA, nB) = enrichment_table(a_labels, b_labels, a_name, b_name)
    if tbl.empty:
        return None

    plot_df = tbl.iloc[::-1]
    fig, ax = plt.subplots(figsize=(12, max(4.2, len(plot_df) * 0.46 + 2.0)))
    y = np.arange(len(plot_df))
    fcs = plot_df["log2FC"].to_numpy(dtype=float)
    sig = (plot_df["significant_q<0.05"] == "yes").to_numpy()
    colors = [SOURCE_COLORS.get(a_name, "#2f6fb0") if f > 0 else SOURCE_COLORS.get(b_name, "#c07219")
              for f in fcs]
    bars = ax.barh(y, fcs, color=colors, edgecolor="white", linewidth=0.5)
    for bar, is_sig in zip(bars, sig):
        bar.set_alpha(1.0 if is_sig else 0.38)
    ax.axvline(0, color="#444", linewidth=1)
    for yy, f, row_sig, a_n, b_n in zip(y, fcs, sig,
                                        plot_df[f"{a_name} n"], plot_df[f"{b_name} n"]):
        off = 0.03 * max(abs(fcs).max(), 0.5)
        ha = "left" if f >= 0 else "right"
        ax.text(f + (off if f >= 0 else -off), yy,
                f"{f:+.2f}  ({int(a_n):,} vs {int(b_n):,}){'' if row_sig else '  n.s.'}",
                va="center", ha=ha, fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["outcome"].tolist(), fontsize=9)
    lim = max(abs(fcs).max() * 1.75, 0.9)
    ax.set_xlim(-lim, lim)
    ax.set_xlabel(f"log2 fold-change   (>0 enriched in {a_name}   ·   <0 enriched in {b_name})",
                  fontsize=9.5)
    ax.grid(axis="x", alpha=0.2, linestyle="--")
    ax.set_axisbelow(True)
    subtitle = f"{a_name} n={nA:,} clusters · {b_name} n={nB:,} clusters"
    if np.isfinite(v):
        subtitle += f" · chi2 p={_fmt_p(p)}, Cramer's V={v:.3f}"
    fig.suptitle(f"{page_title}\n{subtitle}", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    rs._savefig_to_pdf(pdf, fig, page_title)

    slug = f"{a_name}_vs_{b_name}".lower().replace(" ", "_")
    rs._save_table(tbl, f"enrichment_{slug}.csv",
                   title=f"Enrichment — {a_name} vs {b_name} — counts & statistics")
    rs._save_table(pd.DataFrame([{
        "comparison": f"{a_name} vs {b_name}",
        "n_clusters_A": nA, "n_clusters_B": nB,
        "chi2": round(chi2, 2) if np.isfinite(chi2) else "n/a",
        "p": _fmt_p(p) if np.isfinite(p) else "n/a",
        "Cramers_V": round(v, 4) if np.isfinite(v) else "n/a",
        "n_outcomes_tested": len(tbl),
        "n_significant_q<0.05": int((tbl["significant_q<0.05"] == "yes").sum()),
    }]), f"enrichment_{slug}_test.csv",
        title=f"Enrichment — {a_name} vs {b_name} — statistical outcome")
    return tbl


# ══════════════════════════════════════════════════════════════════════════════
# Driver
# ══════════════════════════════════════════════════════════════════════════════
def build_enhanced_report(runs, pdf_path, title, enrichment_pairs=None, top_n_domains=20):
    """
    runs: list of (label, path, specie_filter) - same shape as rs._run_pdf_report.
    enrichment_pairs: list of (a_label, b_label) using DISPLAY names (e.g. "Reads Human").
    """
    prepared = [(*rs._load_and_prepare(path, label, specie),) for label, path, specie in runs]
    label_dfs = {disp(label): df for df, label in prepared}

    # Patch the heatmap so per-source analyses show counts per line.
    rs.splice_type_vs_outcome = splice_type_vs_outcome_with_counts

    print("\nCollapsing clusters for every source ...")
    clabels = {label: cluster_label_table(df) for label, df in label_dfs.items()}

    with PdfPages(pdf_path) as pdf:
        rs._pdf_title_page(pdf, title, fontsize=24)
        rs._ACTIVE_PDF = pdf
        rs._DEFERRED_TABLES = []
        try:
            # 1. event distribution: pies (with legends) then grouped bars
            event_distribution_pies_page(
                pdf, label_dfs, clabels, "unanalyzable",
                "Event-Type Distribution — Unanalyzable clusters")
            event_distribution_pies_page(
                pdf, label_dfs, clabels, "analyzable",
                "Event-Type Distribution — Analyzable clusters")
            event_distribution_bars_page(
                pdf, label_dfs, clabels, "unanalyzable",
                "Event-Type Distribution (bars) — Unanalyzable clusters")
            event_distribution_bars_page(
                pdf, label_dfs, clabels, "analyzable",
                "Event-Type Distribution (bars) — Analyzable clusters")

            # 2. candidate transcripts & tie-break rule
            transcript_choice_page(
                pdf, label_dfs,
                "Possible transcripts per cluster & how the representative was chosen")

            # 3. top domains across all sources
            cross_source_domain_pages(pdf, label_dfs, top_n=top_n_domains)

            # 4. enrichment / depletion
            if enrichment_pairs:
                rs._pdf_title_page(pdf, "Statistical enrichment / depletion", fontsize=20)
                for a_name, b_name in enrichment_pairs:
                    if a_name not in clabels or b_name not in clabels:
                        rs._warn(f"enrichment: unknown source pair {a_name} vs {b_name} - skipped.")
                        continue
                    an = clabels[a_name][clabels[a_name]["cluster_status"] == "analyzable"]
                    bn = clabels[b_name][clabels[b_name]["cluster_status"] == "analyzable"]
                    enrichment_page(pdf, an, bn, a_name, b_name,
                                    f"Enrichment — {a_name} vs {b_name} (analyzable clusters)")

            # 5. the standard per-source sections
            for df, label in prepared:
                d = disp(label)
                rs._pdf_title_page(pdf, d, fontsize=20)
                rs._run_analyses(df, d, skip_sections=["event_distribution"])

            if rs._DEFERRED_TABLES:
                rs._pdf_title_page(pdf, "Tables", fontsize=24)
                for tdf, ttitle in rs._DEFERRED_TABLES:
                    rs._pdf_csv_pages_from_df(pdf, tdf, ttitle)
        finally:
            rs._ACTIVE_PDF = None
            rs._DEFERRED_TABLES = None

    print(f"\n  → {pdf_path}")
    return label_dfs
