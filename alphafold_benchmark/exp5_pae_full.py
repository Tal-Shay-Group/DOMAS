"""Task 5 at FULL scale: PAE features for every downloaded canonical structure,
merged with the full benchmark. Confirms the sample result on all pairs.
Reports AUC / acc / R2 for both targets."""
import warnings; warnings.filterwarnings("ignore")
import json, os, numpy as np, pandas as pd
from exp_common import load_master, load_humsavar, add_overlap_label, cv_metrics, fmt, STRUCT, SP

PAE = f"{SP}/pae"
d = load_master()
d = d[d['kind'] != 'insertion'].dropna(subset=['canon_lo']).reset_index(drop=True)
d['canon_lo'] = d['canon_lo'].astype(int); d['canon_hi'] = d['canon_hi'].astype(int)

cache = {}
def get_M(acc):
    if acc in cache: return cache[acc]
    p = f"{PAE}/{acc}.json"; M = None
    if os.path.exists(p) and os.path.getsize(p) > 100:
        try:
            j = json.load(open(p)); j = j[0] if isinstance(j, list) else j
            A = j.get('predicted_aligned_error') or j.get('pae')
            if A is not None: M = np.asarray(A, np.float32)
        except Exception: M = None
    cache[acc] = M; return M

feats = []
for acc, grp in d.groupby('acc'):
    M = get_M(acc)
    if M is None or M.ndim != 2:
        cache[acc] = None; continue
    n = M.shape[0]; g = M.mean()
    for _, r in grp.iterrows():
        lo, hi = r['canon_lo'] - 1, r['canon_hi']
        if hi > n: continue
        reg = np.zeros(n, bool); reg[lo:hi] = True; rest = ~reg
        if reg.sum() == 0 or rest.sum() == 0: continue
        feats.append({'iso': r['iso'],
                      'pae_reg2rest': float((M[np.ix_(reg, rest)].mean() + M[np.ix_(rest, reg)].mean()) / 2),
                      'pae_global': float(g), 'protL': n})
    cache[acc] = None  # free
pf = pd.DataFrame(feats); pf.to_csv(f"{SP}/pae_features_full.csv", index=False)
d = d.merge(pf, on='iso', how='inner').reset_index(drop=True)
print(f"FULL PAE set: n={len(d)} pairs  proteins={d['acc'].nunique()}  base(TM<0.5)={d['y_tm'].mean():.3f}")

PAEF = ['pae_reg2rest', 'pae_global']
print("\n=== STRUCTURAL (TM<0.5; R2 on continuous TM) — gene-grouped CV, FULL ===")
print(f"  baseline (6 feat, incl burial)   {fmt(cv_metrics(d, STRUCT, 'y_tm', 'tm'))}")
print(f"  + pae_reg2rest                   {fmt(cv_metrics(d, STRUCT+['pae_reg2rest'], 'y_tm', 'tm'))}")
print(f"  + pae_global                     {fmt(cv_metrics(d, STRUCT+['pae_global'], 'y_tm', 'tm'))}")
print(f"  + both PAE                       {fmt(cv_metrics(d, STRUCT+PAEF, 'y_tm', 'tm'))}")
print(f"  + both PAE + protL               {fmt(cv_metrics(d, STRUCT+PAEF+['protL'], 'y_tm', 'tm'))}")

patho, _ = load_humsavar()
sub = add_overlap_label(d, patho)
print("\n=== FUNCTIONAL (pathogenic overlap; R2 = label-variance) — FULL ===")
print(f"  n={len(sub)} pos={int(sub['y'].sum())}")
print(f"  baseline (6 feat)                {fmt(cv_metrics(sub, STRUCT, 'y'))}")
print(f"  + both PAE                       {fmt(cv_metrics(sub, STRUCT+PAEF, 'y'))}")
open(f"{SP}/exp5full_done.txt", "w").write("done")
