import warnings; warnings.filterwarnings("ignore")
import re, gzip, time, numpy as np, pandas as pd, torch, esm
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
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
canon=load_fa(f"{DATA}/uniprot/UP000005640_9606.fasta.gz")
acc2gene={}
with gzip.open(f"{DATA}/uniprot/UP000005640_9606.fasta.gz","rt") as fh:
    for ln in fh:
        if ln.startswith('>'):
            m=re.search(r' GN=(\S+)',ln)
            if m: acc2gene[ln.split('|')[1]]=m.group(1)
g2l={}
with gzip.open(f"{SP}/gnomad_constraint.txt.bgz","rt") as fh:
    h=fh.readline().rstrip().split('\t'); gi,li=h.index('gene'),h.index('oe_lof_upper')
    for ln in fh:
        c=ln.rstrip('\n').split('\t')
        try:
            if c[li] not in('','NA'): g2l[c[gi]]=float(c[li])
        except: pass
bf=pd.read_csv(f"{SP}/full_features.csv"); br=pd.read_csv(f"{SP}/bench_rich_results.csv")[['iso','region_am','max_cov_loss']]
fp=pd.read_csv(f"{SP}/full_pairs.csv").rename(columns={'isoform':'iso'})
d=bf.merge(br,on='iso').merge(fp[['iso','canon_lo','canon_hi']],on='iso').dropna(subset=['mean_rsa']).copy()
d['canon_lo']=d['canon_lo'].astype(int); d['canon_hi']=d['canon_hi'].astype(int)
d['loeuf']=d['acc'].map(lambda a:g2l.get(acc2gene.get(a)))
d['y']=(d['tm']<0.5).astype(int)
# stratified sample ~3000
rng=np.random.default_rng(0)
d=d[d['acc'].isin(canon)].reset_index(drop=True)
chg=d[d['y']==1]; pre=d[d['y']==0]
samp=pd.concat([chg.sample(min(1500,len(chg)),random_state=0), pre.sample(min(1500,len(pre)),random_state=0)]).reset_index(drop=True)
print(f"sample n={len(samp)}  proteins={samp['acc'].nunique()}",flush=True)

model,alphabet=esm.pretrained.esm2_t30_150M_UR50D(); model.eval(); bc=alphabet.get_batch_converter()
def embed_region(acc,lo,hi):
    s=canon[acc]; a=max(0,lo-1-30); b=min(len(s),hi+30); win=s[a:b][:1022]
    _,_,toks=bc([("x",win)])
    with torch.no_grad(): r=model(toks,repr_layers=[30])["representations"][30][0,1:len(win)+1]
    return r.mean(0).numpy()
t0=time.time(); E=np.zeros((len(samp),640),dtype=np.float32)
for i,row in samp.iterrows():
    E[i]=embed_region(row['acc'],row['canon_lo'],row['canon_hi'])
    if (i+1)%500==0: print(f"  embedded {i+1}/{len(samp)} at {time.time()-t0:.0f}s",flush=True)
print(f"embedding done {time.time()-t0:.0f}s",flush=True)

F=['ident','buried_frac','mean_rsa','region_am','loeuf','max_cov_loss']
y=samp['y'].values.astype(float); grp=samp['acc'].values
Xb=samp[F].values
def cvauc(getX, esm_pca=None):
    oof=np.full(len(samp),np.nan)
    for tr,te in GroupKFold(5).split(Xb,y,grp):
        pipe=make_pipeline(SimpleImputer(strategy='median'),StandardScaler(),LogisticRegression(max_iter=3000,C=0.2))
        Xtr,Xte=getX(tr,te)
        pipe.fit(Xtr,y[tr]); oof[te]=pipe.predict_proba(Xte)[:,1]
    return roc_auc_score(y,oof)
def base(tr,te): return Xb[tr], Xb[te]
def esm_only(tr,te): return E[tr], E[te]
def combo(tr,te): return np.hstack([Xb[tr],E[tr]]), np.hstack([Xb[te],E[te]])
def combo_pca(tr,te):
    p=PCA(50).fit(E[tr]); return np.hstack([Xb[tr],p.transform(E[tr])]), np.hstack([Xb[te],p.transform(E[te])])
print("\ngene-grouped 5-fold CV AUC (structural, TM<0.5):",flush=True)
print(f"  our 6 features (baseline)   {cvauc(base):.3f}",flush=True)
print(f"  ESM region embedding only   {cvauc(esm_only):.3f}",flush=True)
print(f"  6 features + ESM (640)      {cvauc(combo):.3f}",flush=True)
print(f"  6 features + ESM-PCA50      {cvauc(combo_pca):.3f}",flush=True)
open(f"{SP}/esm_done.txt","w").write("done")
