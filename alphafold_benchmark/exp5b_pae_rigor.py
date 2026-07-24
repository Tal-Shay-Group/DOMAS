"""Task 5 rigor: the +0.131 AUC PAE jump is huge -> rule out confounds.
TM-score has a known length normalisation; pae_global is a whole-protein property that
could proxy protein size/architecture. Control for identity, protein length L, region
length; decompose which PAE feature drives it; check region-specific PAE survives."""
import warnings; warnings.filterwarnings("ignore")
import json, os, numpy as np, pandas as pd
from exp_common import load_master, cv_metrics, fmt, STRUCT, SP

PAE = f"{SP}/pae"; samp = pd.read_csv(f"{SP}/pae_sample.csv")
samp['canon_lo'] = samp['canon_lo'].astype(int); samp['canon_hi'] = samp['canon_hi'].astype(int)
def load_pae(acc):
    p = f"{PAE}/{acc}.json"
    if not os.path.exists(p) or os.path.getsize(p) < 100: return None
    try:
        j = json.load(open(p)); j = j[0] if isinstance(j, list) else j
        M = j.get('predicted_aligned_error') or j.get('pae')
        return np.asarray(M, np.float32) if M is not None else None
    except Exception: return None
rows = []
for _, r in samp.iterrows():
    M = load_pae(r['acc'])
    if M is None or M.ndim != 2: continue
    n = M.shape[0]; lo, hi = r['canon_lo'] - 1, r['canon_hi']
    if hi > n: continue
    reg = np.zeros(n, bool); reg[lo:hi] = True; rest = ~reg
    if reg.sum() == 0 or rest.sum() == 0: continue
    rows.append({'iso': r['iso'], 'protL': n,
                 'pae_reg2rest': float((M[np.ix_(reg, rest)].mean() + M[np.ix_(rest, reg)].mean()) / 2),
                 'pae_within': float(M[np.ix_(reg, reg)].mean()), 'pae_global': float(M.mean())})
pf = pd.DataFrame(rows)
d = load_master().merge(pf, on='iso', how='inner').reset_index(drop=True)
print(f"n={len(d)}")

def sp(a, b): return d[[a, b]].corr('spearman').iloc[0, 1]
print("\nSpearman with TM, and with potential confounds:")
for c in ['pae_reg2rest', 'pae_within', 'pae_global', 'protL', 'reglen', 'ident_rt', 'buried_frac']:
    print(f"  {c:13} vs TM {sp(c,'tm'):+.3f}   vs protL {sp(c,'protL'):+.3f}   vs reglen {sp(c,'reglen'):+.3f}")

# partial Spearman of pae_reg2rest vs TM, controlling a set of covariates (rank-linear residuals)
def partial(y, x, ctrls):
    Z = d[ctrls].fillna(d[ctrls].median()).rank().values
    Z = np.column_stack([np.ones(len(d)), Z])
    def res(v):
        v = pd.Series(v).rank().values; b, *_ = np.linalg.lstsq(Z, v, rcond=None); return v - Z @ b
    return pd.Series(res(d[y].values)).corr(pd.Series(res(d[x].values)))
for ctrl in [['buried_frac'], ['ident_rt'], ['ident_rt', 'reglen', 'buried_frac'],
             ['ident_rt', 'reglen', 'buried_frac', 'protL']]:
    print(f"  partial(pae_reg2rest, TM | {'+'.join(ctrl)}) = {partial('tm','pae_reg2rest',ctrl):+.3f}")
    print(f"  partial(pae_global,   TM | {'+'.join(ctrl)}) = {partial('tm','pae_global',ctrl):+.3f}")

print("\n=== STRUCTURAL incremental (AUC / acc / R2) ===")
combos = [("baseline (6 feat)", STRUCT),
          ("baseline + protL", STRUCT + ['protL']),
          ("baseline + pae_global only", STRUCT + ['pae_global']),
          ("baseline + pae_reg2rest only", STRUCT + ['pae_reg2rest']),
          ("baseline + pae_within only", STRUCT + ['pae_within']),
          ("baseline + protL + all PAE", STRUCT + ['protL', 'pae_reg2rest', 'pae_within', 'pae_global'])]
for name, F in combos:
    print(f"  {name:32} {fmt(cv_metrics(d, F, 'y_tm', 'tm'))}")
