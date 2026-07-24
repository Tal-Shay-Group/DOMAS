"""Task 1: does the POSITION of the changed region (N-term / C-term / internal)
add signal beyond the existing features?  Position is free from canon_lo/canon_hi/L.
Reports AUC, accuracy@0.5, R2 for both targets (structural R2 = continuous TM)."""
import numpy as np, pandas as pd
from exp_common import load_master, load_humsavar, add_overlap_label, cv_metrics, fmt, STRUCT

d = load_master()
d['rel_center'] = ((d['canon_lo'] + d['canon_hi']) / 2) / d['L']
d['dist_term'] = np.minimum(d['canon_lo'] - 1, d['L'] - d['canon_hi']) / d['L']
d['frac_before'] = (d['canon_lo'] - 1) / d['L']
d['is_terminal'] = (d['dist_term'] < 0.05).astype(float)
POS = ['rel_center', 'dist_term', 'frac_before', 'is_terminal']

print(f"master n={len(d)}  proteins={d['acc'].nunique()}")
print(f"terminal (<5pct from an end) = {100*d['is_terminal'].mean():.1f}pct")

print("\n=== STRUCTURAL (TM<0.5; R2 on continuous TM) — gene-grouped CV ===")
print(f"  baseline (6 feat)   {fmt(cv_metrics(d, STRUCT, 'y_tm', 'tm'))}")
print(f"  + dist_term         {fmt(cv_metrics(d, STRUCT+['dist_term'], 'y_tm', 'tm'))}")
print(f"  + all position      {fmt(cv_metrics(d, STRUCT+POS, 'y_tm', 'tm'))}")

patho, benign = load_humsavar()
dni = d[d['kind'] != 'insertion'].copy()
print("\n=== FUNCTIONAL (pathogenic overlap; R2 = label-variance) — gene-grouped CV ===")
for name, tab in [('PATHOGENIC', patho), ('BENIGN (control)', benign)]:
    sub = add_overlap_label(dni, tab); base = ['reglen', 'ident_rt']
    print(f"  {name:16} n={len(sub)} pos={int(sub['y'].sum())}")
    print(f"    base(len+id)      {fmt(cv_metrics(sub, base, 'y'))}")
    print(f"    + position        {fmt(cv_metrics(sub, base+POS, 'y'))}")
