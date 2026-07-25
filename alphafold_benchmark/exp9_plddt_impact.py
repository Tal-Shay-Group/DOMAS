"""Does region_plddt earn a slot in the REAL impact_prob model?  Tested against the
actual functional feature set [region_am, loeuf, max_cov_loss, buried_frac] (NOT the
structural baseline exp8 used), with drop-column importance vs burial (the correlated
incumbent) and a benign specificity control (a real functional feature must enrich
pathogenic overlap but NOT benign)."""
import warnings; warnings.filterwarnings("ignore")
import sqlite3, numpy as np, pandas as pd
from exp_common import load_master, load_humsavar, add_overlap_label, cv_metrics, fmt, SP

ENR = "/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-db/data/enrichment.sqlite"
d = load_master()
d = d[d['kind'] != 'insertion'].dropna(subset=['canon_lo']).reset_index(drop=True)
d['canon_lo'] = d['canon_lo'].astype(int); d['canon_hi'] = d['canon_hi'].astype(int)
con = sqlite3.connect(ENR); cache = {}
def parr(acc):
    if acc in cache: return cache[acc]
    r = con.execute("SELECT plddt FROM afdb_plddt WHERE accession=?", (acc,)).fetchone()
    cache[acc] = [float(x) for x in r[0].split(',')] if r else None
    return cache[acc]
def rplddt(acc, lo, hi):
    a = parr(acc)
    if not a or len(a) < hi: return np.nan
    seg = a[lo-1:hi]; return round(sum(seg)/len(seg), 2) if seg else np.nan
d['region_plddt'] = [rplddt(a, lo, hi) for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]

IMPACT = ['region_am', 'loeuf', 'max_cov_loss', 'buried_frac']   # the REAL impact_prob features
patho, benign = load_humsavar()
print(f"corr(region_plddt, buried_frac) = {d[['region_plddt','buried_frac']].corr().iloc[0,1]:+.3f}  "
      f"(both = ordered/structured context)")

for name, tab in [('PATHOGENIC', patho), ('BENIGN (specificity control)', benign)]:
    sub = add_overlap_label(d, tab)
    base = cv_metrics(sub, IMPACT, 'y')
    plus = cv_metrics(sub, IMPACT + ['region_plddt'], 'y')
    # drop-column importance within the augmented model
    d_plddt = plus['auc'] - cv_metrics(sub, IMPACT, 'y')['auc']
    d_burial = plus['auc'] - cv_metrics(sub, [f for f in IMPACT if f != 'buried_frac'] + ['region_plddt'], 'y')['auc']
    print(f"\n=== {name}  n={len(sub)} pos={int(sub['y'].sum())} ===")
    print(f"  impact_prob features         {fmt(base)}")
    print(f"  + region_plddt               {fmt(plus)}")
    print(f"  drop-col dAUC region_plddt   {d_plddt:+.3f}")
    print(f"  drop-col dAUC buried_frac    {d_burial:+.3f}   (for contrast: burial's own marginal)")
    # region_plddt alone on this label
    print(f"  region_plddt alone           {fmt(cv_metrics(sub, ['region_plddt'], 'y'))}")
