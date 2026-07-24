import warnings; warnings.filterwarnings("ignore")
import re, gzip, time, sys, numpy as np, pandas as pd, torch, esm
sys.path.insert(0,"/Users/arielmelchior/Documents/projects/DOMAS/code")
from enrichment import Enricher
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, r2_score
SP="/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"
DATA="/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-db/data"
def load_fa(p):
    o={}; op=gzip.open if p.endswith('.gz') else open
    with op(p,'rt') as f:
        cur=None;buf=[]
        for ln in f:
            if ln.startswith('>'):
                if cur:o[cur]=''.join(buf)
                cur=ln.split('|')[1] if '|' in ln else ln[1:].split()[0];buf=[]
            else:buf.append(ln.strip())
        if cur:o[cur]=''.join(buf)
    return o
canon=load_fa(f"{DATA}/uniprot/UP000005640_9606.fasta.gz"); iso=load_fa(f"{DATA}/uniprot/uniprot_sprot_varsplic.fasta.gz")
meta=pd.read_csv(f"{SP}/esm_full_meta.csv")
# fixed 3000-row subset (reuse indices so 150M embeddings align)
rng=np.random.default_rng(1); idx=np.sort(rng.choice(len(meta),3000,replace=False))
d=meta.iloc[idx].reset_index(drop=True)
Ec150=np.load(f"{SP}/Ec.npy")[idx]; Ed150=np.load(f"{SP}/Ed.npy")[idx]
print(f"subset n={len(d)} proteins={d['acc'].nunique()}",flush=True)
model,alphabet=esm.pretrained.esm2_t33_650M_UR50D(); model.eval(); bc=alphabet.get_batch_converter()
def emb(seq,lo,hi):
    a=max(0,lo-1-30); b=min(len(seq),hi+30); win=seq[a:b][:1022]
    if not win: return np.zeros(1280,np.float32)
    _,_,toks=bc([("x",win)])
    with torch.no_grad(): r=model(toks,repr_layers=[33])["representations"][33][0,1:len(win)+1]
    return r.mean(0).numpy()
t0=time.time(); Ec=np.zeros((len(d),1280),np.float32); Ed=np.zeros((len(d),1280),np.float32)
for i,row in d.iterrows():
    cs=canon.get(row['acc']); a=iso.get(row['iso']); reg=Enricher.changed_region(cs,a) if (cs and a) else None
    ec=emb(cs,*reg['canon']) if (reg and reg['canon']) else np.zeros(1280,np.float32)
    ea=emb(a,*reg['alt']) if (reg and reg['alt']) else ec
    Ec[i]=ec; Ed[i]=ec-ea
    if (i+1)%500==0: print(f"  650M embedded {i+1}/{len(d)} at {time.time()-t0:.0f}s",flush=True)
print(f"650M embedding done {time.time()-t0:.0f}s",flush=True)
np.save(f"{SP}/Ec650.npy",Ec); np.save(f"{SP}/Ed650.npy",Ed)
F=['ident','buried_frac','mean_rsa','region_am','loeuf','max_cov_loss']
Xb=d[F].values; tm=d['tm'].values.astype(float); y=(tm<0.5).astype(int); grp=d['acc'].values
folds=list(GroupKFold(5).split(Xb,y,grp))
def cv(EcM,EdM,nc,est,tgt):
    oof=np.full(len(d),np.nan)
    for tr,te in folds:
        parts_tr=[Xb[tr]]; parts_te=[Xb[te]]
        if EcM is not None:
            pc=PCA(nc).fit(EcM[tr]); pd_=PCA(nc).fit(EdM[tr])
            parts_tr+=[pc.transform(EcM[tr]),pd_.transform(EdM[tr])]; parts_te+=[pc.transform(EcM[te]),pd_.transform(EdM[te])]
        m=make_pipeline(SimpleImputer(strategy='median'),StandardScaler(),est())
        m.fit(np.hstack(parts_tr),tgt[tr]); pr=m.predict(np.hstack(parts_te)) if isinstance(est(),Ridge) else m.predict_proba(np.hstack(parts_te))[:,1]
        oof[te]=pr
    return oof
LR=lambda:LogisticRegression(max_iter=3000,C=0.2); RG=lambda:Ridge(alpha=5.0)
print("\n=== 150M vs 650M on identical 3000-row subset (gene-grouped CV) ===",flush=True)
print(f"AUC  baseline           {roc_auc_score(y,cv(None,None,0,LR,y)):.3f}",flush=True)
print(f"AUC  +150M region+diff  {roc_auc_score(y,cv(Ec150,Ed150,50,LR,y)):.3f}",flush=True)
print(f"AUC  +650M region+diff  {roc_auc_score(y,cv(Ec,Ed,50,LR,y)):.3f}",flush=True)
print(f"R2   baseline           {r2_score(tm,cv(None,None,0,RG,tm)):.3f}",flush=True)
print(f"R2   +150M region+diff  {r2_score(tm,cv(Ec150,Ed150,50,RG,tm)):.3f}",flush=True)
print(f"R2   +650M region+diff  {r2_score(tm,cv(Ec,Ed,50,RG,tm)):.3f}",flush=True)
open(f"{SP}/esm650_done.txt","w").write("done")
