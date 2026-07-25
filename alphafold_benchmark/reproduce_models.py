#!/usr/bin/env python3
"""Reproduce and evaluate the two DOMAS calibrated scores end to end.

Retrieves every input needed to reproduce the benchmark EVALUATION (no DoChaP DB
required) and recomputes, for both models:
  * refit coefficients (gene-grouped) - should match the shipped utils models
  * overall AUC / accuracy@0.5 / R2
  * the same metrics WITH the suggested routing band applied (abstain on the middle)
  * a predicted-probability decile calibration table (+ in-bin accuracy)
  * calibration + per-decile accuracy graphs (PNG)

DATA SOURCES
  Downloaded fresh by this script into ./repro_data/ (small):
    - UniProt human reference proteome + varsplic  -> sequences (identity, protL, region)
    - UniProt humsavar                              -> pathogenic/benign labels (functional)
    - gnomAD v2.1.1 by-gene constraint              -> LOEUF
  Reused from the committed per-pair feature CSVs in this folder (each row = one
  canonical-vs-isoform splicing event; provenance in DATA_SOURCES.md):
    - full_pairs.csv          : tm (TM-score label), canon_lo/hi (changed region)
    - pae_features_full.csv   : pae_global (AFDB PAE), protL
    - bench_rich_results.csv  : region_am (AlphaMissense), max_cov_loss (Pfam)
    - full_features.csv       : buried_frac (AFDB RSA)

Run:  python3 reproduce_models.py
"""
import os, sys, gzip, re, io, urllib.request, collections, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, accuracy_score, r2_score

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "repro_data"); os.makedirs(DATA, exist_ok=True)

URLS = {
    "UP000005640_9606.fasta.gz": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz",
    "uniprot_sprot_varsplic.fasta.gz": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz",
    "humsavar.txt": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt",
    "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz": "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
}
def fetch(name):
    p = os.path.join(DATA, name)
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        print(f"  downloading {name} ...", flush=True)
        urllib.request.urlretrieve(URLS[name], p)
    return p

def load_fasta(p):
    out = {}; op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as f:
        cur = None; buf = []
        for ln in f:
            if ln.startswith(">"):
                if cur: out[cur] = "".join(buf)
                h = ln[1:]; cur = h.split("|")[1] if "|" in h else h.split()[0]; buf = []
            else: buf.append(ln.strip())
        if cur: out[cur] = "".join(buf)
    return out

def build_table():
    print("Retrieving inputs...", flush=True)
    canon = load_fasta(fetch("UP000005640_9606.fasta.gz"))
    iso = load_fasta(fetch("uniprot_sprot_varsplic.fasta.gz"))
    acc2gene = {}
    with gzip.open(fetch("UP000005640_9606.fasta.gz"), "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                m = re.search(r" GN=(\S+)", ln)
                if m: acc2gene[ln.split("|")[1]] = m.group(1)
    g2l = {}
    with gzip.open(fetch("gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"), "rt") as fh:
        h = fh.readline().rstrip().split("\t"); gi, li = h.index("gene"), h.index("oe_lof_upper")
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            try:
                if c[li] not in ("", "NA"): g2l[c[gi]] = float(c[li])
            except Exception: pass
    patho, benign = collections.defaultdict(set), collections.defaultdict(set)
    rx = re.compile(r"p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}"); started = False
    for ln in open(fetch("humsavar.txt")):
        if ln.startswith("_______"): started = True; continue
        if not started: continue
        f = ln.split()
        if len(f) < 5: continue
        m = rx.search(f[3])
        if not m: continue
        if f[4] == "LP/P": patho[f[1]].add(int(m.group(1)))
        elif f[4] == "LB/B": benign[f[1]].add(int(m.group(1)))

    fp = pd.read_csv(os.path.join(HERE, "full_pairs.csv")).rename(columns={"isoform": "iso"})
    pae = pd.read_csv(os.path.join(HERE, "pae_features_full.csv"))
    br = pd.read_csv(os.path.join(HERE, "bench_rich_results.csv"))[["iso", "region_am", "max_cov_loss"]]
    ff = pd.read_csv(os.path.join(HERE, "full_features.csv"))[["iso", "buried_frac"]]
    d = (fp[["iso", "acc", "tm", "canon_lo", "canon_hi", "kind"]]
         .merge(pae, on="iso").merge(br, on="iso").merge(ff, on="iso"))
    d = d[d["kind"] != "insertion"].dropna(subset=["canon_lo", "pae_global"]).reset_index(drop=True)
    d["canon_lo"] = d["canon_lo"].astype(int); d["canon_hi"] = d["canon_hi"].astype(int)

    def ident_rt(r):
        cs = canon.get(r["acc"]); a = iso.get(r["iso"])
        if not cs or not a: return np.nan
        shared = (r["canon_lo"] - 1) + (len(cs) - r["canon_hi"])
        return 100.0 * shared / max(len(cs), len(a))
    d["identity"] = d.apply(ident_rt, axis=1)
    d["loeuf"] = d["acc"].map(lambda a: g2l.get(acc2gene.get(a)))
    d["y_tm"] = (d["tm"] < 0.5).astype(int)
    def overlap(acc, lo, hi, tab): return int(acc in tab and any(lo <= p <= hi for p in tab[acc]))
    # impact_prob label (matches fit_calibrated.py): among regions overlapping SOME
    # humsavar variant, discriminate pathogenic (1) vs benign (0).
    d["pat"] = [overlap(a, lo, hi, patho) for a, lo, hi in zip(d["acc"], d["canon_lo"], d["canon_hi"])]
    d["ben"] = [overlap(a, lo, hi, benign) for a, lo, hi in zip(d["acc"], d["canon_lo"], d["canon_hi"])]
    d = d.dropna(subset=["identity"]).reset_index(drop=True)
    return d

# ---- gene-grouped logistic (numpy; matches utils fit) + Ridge for continuous R2 ----
def _fit_logistic(X, y, it=4000, lr=0.3, l2=1e-3):
    w = np.zeros(X.shape[1])
    for _ in range(it):
        p = 1 / (1 + np.exp(-np.clip(X @ w, -30, 30)))
        g = X.T @ (p - y) / len(y) + l2 * w; g[0] = (X[:, 0] @ (p - y)) / len(y); w -= lr * g
    return w

def evaluate(d, feats, bin_t, cont_t=None, sub=None):
    df = d if sub is None else d[sub].reset_index(drop=True)
    X = df[feats].copy(); med = X.median(); Xf = X.fillna(med)
    mu, sd = Xf.mean().values, Xf.std().values + 1e-9
    Xs = np.nan_to_num(((Xf - mu) / sd).values)   # residual NaN (all-NaN col in a subset) -> standardized mean
    yb = df[bin_t].values.astype(float); grp = df["acc"].values
    yc = df[cont_t].values.astype(float) if cont_t else None
    oof = np.full(len(df), np.nan); oofr = np.full(len(df), np.nan)
    for tr, te in GroupKFold(5).split(Xs, yb, grp):
        Xtr = np.column_stack([np.ones(len(tr)), Xs[tr]]); Xte = np.column_stack([np.ones(len(te)), Xs[te]])
        w = _fit_logistic(Xtr, yb[tr]); oof[te] = 1 / (1 + np.exp(-np.clip(Xte @ w, -30, 30)))
        if cont_t:
            b, *_ = np.linalg.lstsq(Xtr.T @ Xtr + 5 * np.eye(Xtr.shape[1]), Xtr.T @ yc[tr], rcond=None)
            oofr[te] = Xte @ b
    wf = _fit_logistic(np.column_stack([np.ones(len(df)), Xs]), yb)   # coefficients on all data
    coef = dict(zip(feats, wf[1:]))
    auc = roc_auc_score(yb, oof); acc = accuracy_score(yb, (oof >= 0.5))
    r2 = r2_score(yc, oofr) if cont_t else r2_score(yb, oof)
    return dict(df=df, oof=oof, oofr=oofr, yb=yb, yc=yc, coef=coef, intercept=wf[0], auc=auc, acc=acc, r2=r2)

def decile_table(oof, yb, yc):
    order = np.argsort(oof); bins = np.array_split(order, 10)
    rows = []
    for i, idx in enumerate(bins):
        p = oof[idx]; y = yb[idx]
        pred = (p >= 0.5).astype(int)
        rows.append(dict(decile=i + 1, lo=p.min(), hi=p.max(), mean_pred=p.mean(),
                         observed=y.mean(), n=len(idx), in_bin_acc=(pred == y).mean()))
    return pd.DataFrame(rows)

def banded(oof, oofr, yb, yc, lo, hi):
    """Metrics on the CALLED set after routing the [lo,hi] probability band out."""
    called = (oof <= lo) | (oof >= hi)
    y = yb[called]; p = oof[called]
    out = dict(pct_called=called.mean(),
               auc=roc_auc_score(y, p) if len(set(y)) > 1 else float("nan"),
               acc=accuracy_score(y, (p >= 0.5)),
               r2=(r2_score(yc[called], oofr[called]) if yc is not None else r2_score(y, p)))
    return out

def main():
    d = build_table()
    print(f"\nAssembled {len(d)} pairs, {d['acc'].nunique()} proteins.\n")
    report = {}

    # ---------- STRUCTURAL: fold_change_prob ----------
    F_STRUCT = ["pae_global", "identity", "max_cov_loss", "protL"]
    S = evaluate(d, F_STRUCT, "y_tm", cont_t="tm")
    print("=" * 70)
    print("fold_change_prob (STRUCTURAL, P(TM<0.5))")
    print(f"  n={len(S['df'])}  base rate(TM<0.5)={S['yb'].mean():.3f}  majority acc={max(S['yb'].mean(),1-S['yb'].mean()):.3f}")
    print(f"  coefficients (standardized): intercept={S['intercept']:+.3f}  " +
          "  ".join(f"{k}={v:+.3f}" for k, v in S["coef"].items()))
    print(f"  OVERALL (no band):  AUC={S['auc']:.3f}  acc={S['acc']:.3f}  R2={S['r2']:.3f}")
    b = banded(S["oof"], S["oofr"], S["yb"], S["yc"], 0.40, 0.60)
    print(f"  BAND 0.40-0.60 routed to folding -> CALLED set ({b['pct_called']*100:.0f}%): AUC={b['auc']:.3f}  acc={b['acc']:.3f}  R2={b['r2']:.3f}")
    dt = decile_table(S["oof"], S["yb"], S["yc"]); report["struct_decile"] = dt
    print(dt.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # ---------- FUNCTIONAL: impact_prob (pathogenic-vs-benign, matches fit_calibrated.py) ----------
    F_FUNC = ["region_am", "loeuf", "max_cov_loss", "buried_frac"]
    d["y_disc"] = d["pat"].astype(int)
    sub = ((d["pat"] == 1) | (d["ben"] == 1)).values     # regions overlapping SOME variant
    Fm = evaluate(d, F_FUNC, "y_disc", cont_t=None, sub=sub)
    print("\n" + "=" * 70)
    print("impact_prob (FUNCTIONAL, P(pathogenic | region overlaps a variant))")
    print(f"  n={len(Fm['df'])}  base rate={Fm['yb'].mean():.3f}  majority acc={max(Fm['yb'].mean(),1-Fm['yb'].mean()):.3f}")
    print(f"  coefficients (standardized): intercept={Fm['intercept']:+.3f}  " +
          "  ".join(f"{k}={v:+.3f}" for k, v in Fm["coef"].items()))
    print(f"  OVERALL (no band):  AUC={Fm['auc']:.3f}  acc={Fm['acc']:.3f}  R2={Fm['r2']:.3f}")
    bf = banded(Fm["oof"], Fm["oofr"], Fm["yb"], Fm["yc"], 0.40, 0.60)
    print(f"  BAND 0.40-0.60 abstain -> CALLED set ({bf['pct_called']*100:.0f}%): AUC={bf['auc']:.3f}  acc={bf['acc']:.3f}  R2={bf['r2']:.3f}")
    dtf = decile_table(Fm["oof"], Fm["yb"], None); report["func_decile"] = dtf
    print(dtf.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # ---------- graphs ----------
    for tag, res, dtab, title in [("foldchange", S, report["struct_decile"], "fold_change_prob (structural)"),
                                  ("impact", Fm, report["func_decile"], "impact_prob (functional)")]:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        ax[0].plot([0, 1], [0, 1], "--", color="#999", lw=1)
        ax[0].plot(dtab["mean_pred"], dtab["observed"], "o-", color="#2b6cb0")
        ax[0].set_xlabel("mean predicted probability (decile)"); ax[0].set_ylabel("observed rate")
        ax[0].set_title(f"{title}\ncalibration"); ax[0].set_xlim(0, 1); ax[0].set_ylim(0, 1)
        ax[1].bar(dtab["decile"], dtab["in_bin_acc"], color="#2b6cb0", alpha=0.85)
        ax[1].axhline(max(res["yb"].mean(), 1 - res["yb"].mean()), ls="--", color="#e53e3e",
                      label="majority baseline")
        ax[1].set_xlabel("predicted-probability decile (1=lowest)"); ax[1].set_ylabel("in-bin accuracy @0.5")
        ax[1].set_title("per-decile accuracy"); ax[1].set_ylim(0, 1); ax[1].legend(fontsize=8)
        fig.tight_layout(); out = os.path.join(HERE, f"model_card_{tag}.png")
        fig.savefig(out, dpi=130); plt.close(fig); print(f"\nsaved {out}")

if __name__ == "__main__":
    main()
