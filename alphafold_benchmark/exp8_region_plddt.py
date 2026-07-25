"""Does mean pLDDT of the AFFECTED (changed) region add to the PAE-era fold_change
model?  region_plddt = mean per-residue pLDDT over [canon_lo, canon_hi] of the
canonical AlphaFold model (from afdb_plddt). Region-specific AF confidence -- a
natural companion to pae_global (whole-structure). Tests the current shipped
4-feature model +/- region_plddt, both targets, AUC/acc/R2 + drop-column importance."""
import warnings; warnings.filterwarnings("ignore")
import sqlite3, numpy as np, pandas as pd
from exp_common import load_master, load_humsavar, add_overlap_label, cv_metrics, fmt, SP

ENR = "/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-db/data/enrichment.sqlite"
d = load_master()
d = d[d['kind'] != 'insertion'].dropna(subset=['canon_lo']).reset_index(drop=True)
d = d.drop(columns=['L']).merge(pd.read_csv(f"{SP}/pae_features_full.csv"), on='iso', how='inner').reset_index(drop=True)
d['canon_lo'] = d['canon_lo'].astype(int); d['canon_hi'] = d['canon_hi'].astype(int)

# region mean pLDDT from afdb_plddt (per-residue, comma-joined), cached per accession
con = sqlite3.connect(ENR); cache = {}
def plddt_arr(acc):
    if acc in cache: return cache[acc]
    row = con.execute("SELECT plddt FROM afdb_plddt WHERE accession=?", (acc,)).fetchone()
    cache[acc] = [float(x) for x in row[0].split(',')] if row else None
    return cache[acc]
def region_plddt(acc, lo, hi):
    a = plddt_arr(acc)
    if not a or len(a) < hi: return np.nan
    seg = a[lo-1:hi]
    return round(sum(seg)/len(seg), 2) if seg else np.nan
d['region_plddt'] = [region_plddt(a, lo, hi) for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]
d = d.rename(columns={'ident_rt': 'identity'})
print(f"n={len(d)}  region_plddt computed for {d['region_plddt'].notna().sum()}")
print(f"corr(region_plddt, pae_global)={d[['region_plddt','pae_global']].corr().iloc[0,1]:.3f}  "
      f"Spearman(region_plddt, TM)={d[['region_plddt','tm']].corr('spearman').iloc[0,1]:+.3f}")

SHIP = ['pae_global', 'identity', 'max_cov_loss', 'protL']   # current shipped model
print("\n=== STRUCTURAL (TM<0.5; R2 on continuous TM) — gene-grouped CV ===")
print(f"  shipped 4-feature            {fmt(cv_metrics(d, SHIP, 'y_tm', 'tm'))}")
print(f"  + region_plddt               {fmt(cv_metrics(d, SHIP+['region_plddt'], 'y_tm', 'tm'))}")
print(f"  region_plddt alone           {fmt(cv_metrics(d, ['region_plddt'], 'y_tm', 'tm'))}")
# drop-column importance of region_plddt within the augmented model
full = cv_metrics(d, SHIP+['region_plddt'], 'y_tm', 'tm')
drop = cv_metrics(d, SHIP, 'y_tm', 'tm')
print(f"  drop-column importance of region_plddt: dAUC={full['auc']-drop['auc']:+.3f}  dR2={full['r2']-drop['r2']:+.3f}")

patho, _ = load_humsavar()
sub = add_overlap_label(d, patho)
print("\n=== FUNCTIONAL (pathogenic overlap; R2 = label-variance) ===")
print(f"  n={len(sub)} pos={int(sub['y'].sum())}")
print(f"  shipped 4-feature            {fmt(cv_metrics(sub, SHIP, 'y'))}")
print(f"  + region_plddt               {fmt(cv_metrics(sub, SHIP+['region_plddt'], 'y'))}")
