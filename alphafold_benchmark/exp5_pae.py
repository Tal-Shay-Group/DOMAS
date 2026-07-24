"""Task 5: does PAE (predicted aligned error) add signal beyond burial?  The docs
ASSUMED PAE is redundant (extrapolated from contacts) but never measured it. PAE is
AF's confidence in the RELATIVE positioning of residue pairs -> whether the changed
region is rigidly docked to the core (low region<->rest PAE) or flexibly hung off
(high). That is inter-domain rigidity, which local burial can't see.
Structural target; compares baseline (incl. burial) vs +PAE. AUC/acc/R2."""
import warnings; warnings.filterwarnings("ignore")
import json, os, numpy as np, pandas as pd
from exp_common import load_master, load_humsavar, add_overlap_label, cv_metrics, fmt, STRUCT, SP

PAE = f"{SP}/pae"
samp = pd.read_csv(f"{SP}/pae_sample.csv")
samp['canon_lo'] = samp['canon_lo'].astype(int); samp['canon_hi'] = samp['canon_hi'].astype(int)

def load_pae(acc):
    p = f"{PAE}/{acc}.json"
    if not os.path.exists(p) or os.path.getsize(p) < 100: return None
    try:
        j = json.load(open(p))
        if isinstance(j, list): j = j[0]
        M = j.get('predicted_aligned_error') or j.get('pae')
        return np.asarray(M, dtype=np.float32) if M is not None else None
    except Exception: return None

rows = []
for _, r in samp.iterrows():
    M = load_pae(r['acc'])
    if M is None or M.ndim != 2: continue
    n = M.shape[0]; lo, hi = r['canon_lo'] - 1, r['canon_hi']
    if hi > n: continue
    reg = np.zeros(n, bool); reg[lo:hi] = True; rest = ~reg
    if reg.sum() == 0 or rest.sum() == 0: continue
    r2r = float((M[np.ix_(reg, rest)].mean() + M[np.ix_(rest, reg)].mean()) / 2)
    rows.append({'iso': r['iso'], 'pae_reg2rest': r2r,
                 'pae_within': float(M[np.ix_(reg, reg)].mean()),
                 'pae_global': float(M.mean())})
pf = pd.DataFrame(rows)
print(f"PAE features computed for {len(pf)} / {len(samp)} sampled pairs")

d = load_master().merge(pf, on='iso', how='inner')
print(f"merged n={len(d)} proteins={d['acc'].nunique()}  base(TM<0.5)={d['y_tm'].mean():.3f}")
print("\nSpearman(PAE feature, TM):")
for c in ['pae_reg2rest', 'pae_within', 'pae_global']:
    print(f"  {c:14} {d[[c,'tm']].corr('spearman').iloc[0,1]:+.3f}")
# partial vs burial: residualize both on buried_frac, correlate
from numpy.polynomial import polynomial as P
def resid(y, x):
    x = x.fillna(x.median()); b = np.polyfit(x, y, 1); return y - np.polyval(b, x)
dd = d.dropna(subset=['buried_frac', 'pae_reg2rest', 'tm'])
rp = resid(dd['pae_reg2rest'].reset_index(drop=True), dd['buried_frac'].reset_index(drop=True))
rt = resid(dd['tm'].reset_index(drop=True), dd['buried_frac'].reset_index(drop=True))
print(f"\npartial Spearman(pae_reg2rest, TM | buried_frac) = {pd.Series(rp).corr(pd.Series(rt),'spearman'):+.3f}")

PAEF = ['pae_reg2rest', 'pae_within', 'pae_global']
print("\n=== STRUCTURAL (TM<0.5; R2 on continuous TM) ===")
print(f"  baseline (6 feat, incl burial)   {fmt(cv_metrics(d, STRUCT, 'y_tm', 'tm'))}")
print(f"  + PAE (reg2rest/within/global)   {fmt(cv_metrics(d, STRUCT + PAEF, 'y_tm', 'tm'))}")
print(f"  PAE only (no baseline)           {fmt(cv_metrics(d, PAEF, 'y_tm', 'tm'))}")

patho, _ = load_humsavar()
sub = add_overlap_label(d[d['kind'] != 'insertion'], patho)
print("\n=== FUNCTIONAL (pathogenic overlap; R2 = label-variance) ===")
print(f"  n={len(sub)} pos={int(sub['y'].sum())}")
print(f"  baseline (6 feat)                {fmt(cv_metrics(sub, STRUCT, 'y'))}")
print(f"  + PAE                            {fmt(cv_metrics(sub, STRUCT + PAEF, 'y'))}")
