"""Combined model (our 6 features + ESM region & canonical-alt-difference embeddings,
PCA-50): linear (Ridge) regression on continuous TM, logistic classifier, decile
calibration table, and abstention-band table for routing to real folding.
Requires the embeddings produced by esm_embedding_full.py (Ec.npy, Ed.npy,
esm_full_meta.csv). Gene-grouped 5-fold CV throughout; PCA fit inside each fold."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, r2_score, mean_squared_error
SP = "<dir with Ec.npy / Ed.npy / esm_full_meta.csv>"
d = pd.read_csv(f"{SP}/esm_full_meta.csv"); Ec = np.load(f"{SP}/Ec.npy"); Ed = np.load(f"{SP}/Ed.npy")
F = ['ident', 'buried_frac', 'mean_rsa', 'region_am', 'loeuf', 'max_cov_loss']
Xb = d[F].values; tm = d['tm'].values.astype(float); y = (tm < 0.5).astype(int); grp = d['acc'].values
gkf = GroupKFold(5)
def combined(tr, te):
    pc = PCA(50).fit(Ec[tr]); pdd = PCA(50).fit(Ed[tr])
    return (np.hstack([Xb[tr], pc.transform(Ec[tr]), pdd.transform(Ed[tr])]),
            np.hstack([Xb[te], pc.transform(Ec[te]), pdd.transform(Ed[te])]))
def base(tr, te): return Xb[tr], Xb[te]
def cv(build, est):
    oof = np.full(len(d), np.nan)
    for tr, te in gkf.split(Xb, y, grp):
        Xtr, Xte = build(tr, te)
        m = make_pipeline(SimpleImputer(strategy='median'), StandardScaler(), est())
        yy = tm if isinstance(est(), Ridge) else y
        m.fit(Xtr, yy[tr]); oof[te] = m.predict(Xte) if isinstance(est(), Ridge) else m.predict_proba(Xte)[:, 1]
    return oof
ridge = lambda: Ridge(alpha=5.0); logit = lambda: LogisticRegression(max_iter=3000, C=0.2)
for name, b in [("our 6 features", base), ("+ ESM region+diff", combined)]:
    o = cv(b, ridge); print(f"Ridge->TM  {name:20} R2={r2_score(tm,o):.3f} RMSE={mean_squared_error(tm,o)**0.5:.3f}")
p = cv(combined, logit); print(f"\ncombined logistic AUC={roc_auc_score(y,p):.3f}")
dd = pd.DataFrame({'p': p, 'y': y}); dd['dec'] = pd.qcut(dd['p'], 10, labels=False, duplicates='drop')
print("decile | mean pred | observed change"); 
for b, g in dd.groupby('dec'): print(f"  {int(b)+1:2}  {g['p'].mean():.2f}  {g['y'].mean():.2f}")
print("\nfold-band | %routed | accInBand | %called | accCalled")
for lo, hi in [(45,55),(40,60),(35,65),(30,70),(25,75)]:
    inb = (p>=lo/100)&(p<=hi/100); kept=~inb
    ai = ((p[inb]>=.5).astype(int)==y[inb]).mean(); ac=((p[kept]>=.5).astype(int)==y[kept]).mean()
    print(f"  {lo}-{hi}  {inb.mean():.0%}  {ai:.0%}  {kept.mean():.0%}  {ac:.0%}")
