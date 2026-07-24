import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, accuracy_score
SP="/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"
meta=pd.read_csv(f"{SP}/esm_full_meta.csv")
rng=np.random.default_rng(1); idx=np.sort(rng.choice(len(meta),3000,replace=False))
d=meta.iloc[idx].reset_index(drop=True)
Ec150=np.load(f"{SP}/Ec.npy")[idx]; Ed150=np.load(f"{SP}/Ed.npy")[idx]
Ec650=np.load(f"{SP}/Ec650.npy"); Ed650=np.load(f"{SP}/Ed650.npy")
F=['ident','buried_frac','mean_rsa','region_am','loeuf','max_cov_loss']
Xb=d[F].values; tm=d['tm'].values.astype(float); y=(tm<0.5).astype(int); grp=d['acc'].values
folds=list(GroupKFold(5).split(Xb,y,grp))
LR=lambda:LogisticRegression(max_iter=3000,C=0.2)
def cv(EcM,EdM,nc):
    oof=np.full(len(d),np.nan)
    for tr,te in folds:
        ptr=[Xb[tr]]; pte=[Xb[te]]
        if EcM is not None:
            pc=PCA(nc).fit(EcM[tr]); pdd=PCA(nc).fit(EdM[tr])
            ptr+=[pc.transform(EcM[tr]),pdd.transform(EdM[tr])]; pte+=[pc.transform(EcM[te]),pdd.transform(EdM[te])]
        m=make_pipeline(SimpleImputer(strategy='median'),StandardScaler(),LR())
        m.fit(np.hstack(ptr),y[tr]); oof[te]=m.predict_proba(np.hstack(pte))[:,1]
    return oof
base=(1-y.mean()) if y.mean()<0.5 else y.mean()
print(f"n={len(d)}  base rate TM<0.5 = {y.mean():.3f}  majority-class accuracy = {base:.3f}\n")
print(f"{'model':<26}{'AUC':>7}{'acc@0.5':>10}")
for name,ec,ed,nc in [("baseline (6 feat)",None,None,0),
                      ("+150M region+diff",Ec150,Ed150,50),
                      ("+650M region+diff",Ec650,Ed650,50)]:
    p=cv(ec,ed,nc); pred=(p>=0.5).astype(int)
    print(f"{name:<26}{roc_auc_score(y,p):>7.3f}{accuracy_score(y,pred):>10.3f}")
