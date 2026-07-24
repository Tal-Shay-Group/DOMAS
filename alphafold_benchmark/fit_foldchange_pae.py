"""Refit fold_change_prob on the PAE-era feature set: [pae_global, identity,
max_cov_loss, protL] (chosen by exp7 importance -- pae_global dominates, identity #2,
protL + max_cov_loss small real adds; burial/region_am/loeuf/mean_rsa/pae_reg2rest/
dist_term were ~0 marginal and dropped). Emits coefficients for utils._FOLD_CHANGE_MODEL
and gene-grouped CV AUC/acc/R2. Runtime `identity` = trimmed-region identity (ident_rt)."""
import warnings; warnings.filterwarnings("ignore")
import json, numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, accuracy_score, r2_score
from exp_common import load_master, SP

d = load_master()
d = d[d['kind'] != 'insertion'].dropna(subset=['canon_lo']).reset_index(drop=True)
pf = pd.read_csv(f"{SP}/pae_features_full.csv")
d = d.drop(columns=['L']).merge(pf, on='iso', how='inner').reset_index(drop=True)  # pf carries protL
d = d.rename(columns={'ident_rt': 'identity'})

FEATURES = ['pae_global', 'identity', 'max_cov_loss', 'protL']
d = d.dropna(subset=['pae_global', 'identity', 'protL']).reset_index(drop=True)
X = d[FEATURES].copy()
med = {c: float(X[c].median()) for c in FEATURES}
Xf = X.fillna(pd.Series(med))
y = d['y_tm'].values.astype(float); tm = d['tm'].values.astype(float); grp = d['acc'].values
mu = Xf.mean().values; sd = Xf.std().values + 1e-9
Xs = ((Xf - mu) / sd).values

def fit(Xd, yv, it=4000, lr=0.3, l2=1e-3):
    w = np.zeros(Xd.shape[1])
    for _ in range(it):
        p = 1 / (1 + np.exp(-np.clip(Xd @ w, -30, 30)))
        g = Xd.T @ (p - yv) / len(yv) + l2 * w; g[0] = (Xd[:, 0] @ (p - yv)) / len(yv); w -= lr * g
    return w

# gene-grouped CV (AUC/acc + Ridge R2 on continuous TM)
oof = np.full(len(y), np.nan); oofr = np.full(len(y), np.nan)
for tr, te in GroupKFold(5).split(Xs, y, grp):
    Xtr = np.column_stack([np.ones(len(tr)), Xs[tr]]); Xte = np.column_stack([np.ones(len(te)), Xs[te]])
    w = fit(Xtr, y[tr]); oof[te] = 1 / (1 + np.exp(-np.clip(Xte @ w, -30, 30)))
    # ridge for R2 on continuous TM
    b, *_ = np.linalg.lstsq(Xtr.T @ Xtr + 5.0 * np.eye(Xtr.shape[1]), Xtr.T @ tm[tr], rcond=None)
    oofr[te] = Xte @ b
print(f"n={len(y)} proteins={d['acc'].nunique()}  base(TM<0.5)={y.mean():.3f}  majority-acc={max(y.mean(),1-y.mean()):.3f}")
print(f"gene-grouped CV:  AUC={roc_auc_score(y,oof):.3f}  acc={accuracy_score(y,(oof>=0.5)):.3f}  R2={r2_score(tm,oofr):.3f}")

# final fit on all data -> shipped coefficients
wf = fit(np.column_stack([np.ones(len(y)), Xs]), y)
model = {'features': FEATURES,
         'median': [round(med[c], 4) for c in FEATURES],
         'mean': [round(float(x), 4) for x in mu],
         'std': [round(float(x), 4) for x in sd],
         'coef': [round(float(x), 4) for x in wf[1:]],
         'intercept': round(float(wf[0]), 4)}
open(f"{SP}/foldchange_model_pae.json", "w").write(json.dumps(model, indent=2))
open("/Users/arielmelchior/Documents/projects/DOMAS/alphafold_benchmark/foldchange_model.json", "w").write(json.dumps(model, indent=2))
print("\nMODEL =", json.dumps(model))
print("standardized coefficients:")
for c, w in zip(FEATURES, wf[1:]): print(f"  {c:14} {w:+.3f}")
