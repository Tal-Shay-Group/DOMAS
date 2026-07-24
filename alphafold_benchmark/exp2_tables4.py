"""Task 2: the 21 Table-S4 columns dropped from table_s4_all.csv.
Two groups, very different status:
  (a) ISOFORM structural descriptors (helix/sheet/loop, Rg, charge, IDR%) — derived
      from the ISOFORM's own AF2 fold, so CIRCULAR for predicting TM (you'd have had
      to fold it). Reported as an upper bound, not a usable cheap feature.
  (b) PTM-change counts (added / missed / buried<->exposed) — canonical-vs-isoform
      deltas, NON-circular and functional. The genuinely useful test, vs pathogenic."""
import numpy as np, pandas as pd
from exp_common import load_master, load_humsavar, add_overlap_label, cv_metrics, fmt, STRUCT, SP

d = load_master()
s4 = pd.read_csv(f"{SP}/table_s4_full.csv").rename(columns={'isoform': 'iso'})
ren = {'helix':'iso_helix','sheet':'iso_sheet','loop':'iso_loop','surface charge':'iso_charge',
       'radius of gyration':'iso_rg','IDR percentage (pLDDT)':'iso_idr_plddt',
       'IDR percentage (SPOT)':'iso_idr_spot','non_IDR':'iso_nonidr','non_IDR_helix':'iso_nonidr_helix',
       'charge positive':'d_charge_pos','charge negative':'d_charge_neg',
       'radius positive':'d_rad_pos','radius negative':'d_rad_neg',
       'ptm added':'ptm_added','ptm buried to exposed':'ptm_b2e','ptm exposed to buried':'ptm_e2b',
       'ptm missed':'ptm_missed'}
s4 = s4.rename(columns=ren)
ISO_STRUCT = ['iso_helix','iso_sheet','iso_loop','iso_charge','iso_rg','iso_idr_plddt','iso_idr_spot','iso_nonidr']
PTM = ['ptm_added','ptm_b2e','ptm_e2b','ptm_missed']
for c in ISO_STRUCT+PTM+['d_charge_pos','d_charge_neg','d_rad_pos','d_rad_neg']:
    s4[c] = pd.to_numeric(s4[c], errors='coerce')
d = d.merge(s4[['iso']+ISO_STRUCT+PTM+['d_charge_pos','d_charge_neg','d_rad_pos','d_rad_neg']], on='iso', how='left')
print(f"master n={len(d)}  with Table-S4 structural cols: {d['iso_rg'].notna().sum()}  with PTM cols: {d['ptm_added'].notna().sum()}")

# ---- (a) STRUCTURAL — CIRCULAR upper bound ----
print("\n=== STRUCTURAL (TM<0.5; R2 on continuous TM) — CIRCULAR upper bound ===")
print(f"  baseline (6 cheap feat)            {fmt(cv_metrics(d, STRUCT, 'y_tm', 'tm'))}")
print(f"  + isoform SS/Rg/charge/IDR [CIRC]  {fmt(cv_metrics(d, STRUCT+ISO_STRUCT, 'y_tm', 'tm'))}")
print(f"  + charge/radius deltas [CIRC]      {fmt(cv_metrics(d, STRUCT+['d_charge_pos','d_charge_neg','d_rad_pos','d_rad_neg'], 'y_tm', 'tm'))}")
print("  (circular: these come from the isoform's own AF2 structure; a prospective")
print("   tool would have to fold the isoform to get them — i.e. it'd have TM already.)")

# ---- (b) FUNCTIONAL — PTM-change, NON-circular ----
patho, benign = load_humsavar()
dni = d[d['kind'] != 'insertion'].copy()
print("\n=== FUNCTIONAL (pathogenic overlap; R2 = label-variance) — PTM-change (non-circular) ===")
for name, tab in [('PATHOGENIC', patho), ('BENIGN (control)', benign)]:
    sub = add_overlap_label(dni, tab)
    base = ['reglen','ident_rt','region_am','loeuf']
    print(f"  {name:16} n={len(sub)} pos={int(sub['y'].sum())}")
    print(f"    base(len+id+am+loeuf)  {fmt(cv_metrics(sub, base, 'y'))}")
    print(f"    + PTM-change           {fmt(cv_metrics(sub, base+PTM, 'y'))}")
print("\n  PTM-change rates by pathogenic overlap (mean count in region):")
sub = add_overlap_label(dni, patho)
for c in PTM:
    print(f"    {c:12}  overlap=1: {sub[sub['y']==1][c].mean():.3f}   overlap=0: {sub[sub['y']==0][c].mean():.3f}")
