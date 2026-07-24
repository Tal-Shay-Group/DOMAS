import warnings; warnings.filterwarnings("ignore")
import re, gzip, time, sys, numpy as np, pandas as pd, torch, esm
sys.path.insert(0,"/Users/arielmelchior/Documents/projects/DOMAS/code")
from enrichment import Enricher
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
canon=load_fa(f"{DATA}/uniprot/UP000005640_9606.fasta.gz"); iso=load_fa(f"{DATA}/uniprot/uniprot_sprot_varsplic.fasta.gz")
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
d=bf.merge(br,on='iso').dropna(subset=['mean_rsa']).copy()
d=d[d['acc'].isin(canon)].reset_index(drop=True)
d['loeuf']=d['acc'].map(lambda a:g2l.get(acc2gene.get(a))); d['y']=(d['tm']<0.5).astype(int)
print(f"full n={len(d)} proteins={d['acc'].nunique()}",flush=True)
model,alphabet=esm.pretrained.esm2_t30_150M_UR50D(); model.eval(); bc=alphabet.get_batch_converter()
def emb_window(seq,lo,hi):
    a=max(0,lo-1-30); b=min(len(seq),hi+30); win=seq[a:b][:1022]
    if not win: return np.zeros(640,dtype=np.float32)
    _,_,toks=bc([("x",win)])
    with torch.no_grad(): r=model(toks,repr_layers=[30])["representations"][30][0,1:len(win)+1]
    return r.mean(0).numpy()
t0=time.time(); Ec=np.zeros((len(d),640),np.float32); Ed=np.zeros((len(d),640),np.float32)
for i,row in d.iterrows():
    cs=canon[row['acc']]; a=iso.get(row['iso'])
    reg=Enricher.changed_region(cs,a) if a else None
    if reg and reg['canon']:
        lo,hi=reg['canon']; ec=emb_window(cs,lo,hi)
    else: ec=np.zeros(640,np.float32)
    if reg and reg['alt']:
        alo,ahi=reg['alt']; ea=emb_window(a,alo,ahi)
    else: ea=ec   # deletion: no alt region -> diff 0
    Ec[i]=ec; Ed[i]=ec-ea
    if (i+1)%2000==0: print(f"  embedded {i+1}/{len(d)} at {time.time()-t0:.0f}s",flush=True)
print(f"embedding done {time.time()-t0:.0f}s",flush=True)
np.save(f"{SP}/Ec.npy",Ec); np.save(f"{SP}/Ed.npy",Ed); d.to_csv(f"{SP}/esm_full_meta.csv",index=False)
F=['ident','buried_frac','mean_rsa','region_am','loeuf','max_cov_loss']
y=d['y'].values.astype(float); grp=d['acc'].values; Xb=d[F].values
def cvauc(build):
    oof=np.full(len(d),np.nan)
    for tr,te in GroupKFold(5).split(Xb,y,grp):
        Xtr,Xte=build(tr,te); pipe=make_pipeline(SimpleImputer(strategy='median'),StandardScaler(),LogisticRegression(max_iter=3000,C=0.2))
        pipe.fit(Xtr,y[tr]); oof[te]=pipe.predict_proba(Xte)[:,1]
    return roc_auc_score(y,oof)
def pca_feat(M):
    def b(tr,te):
        p=PCA(50).fit(M[tr]); return np.hstack([Xb[tr],p.transform(M[tr])]),np.hstack([Xb[te],p.transform(M[te])])
    return b
def both(tr,te):
    pc=PCA(50).fit(Ec[tr]); pd_=PCA(50).fit(Ed[tr])
    return (np.hstack([Xb[tr],pc.transform(Ec[tr]),pd_.transform(Ed[tr])]),
            np.hstack([Xb[te],pc.transform(Ec[te]),pd_.transform(Ed[te])]))
print("\nFULL 11k gene-grouped CV AUC (structural TM<0.5):",flush=True)
print(f"  baseline (6 features)            {cvauc(lambda tr,te:(Xb[tr],Xb[te])):.3f}",flush=True)
print(f"  + region-embed PCA50             {cvauc(pca_feat(Ec)):.3f}",flush=True)
print(f"  + canon-alt DIFF PCA50           {cvauc(pca_feat(Ed)):.3f}",flush=True)
print(f"  + region + diff PCA50            {cvauc(both):.3f}",flush=True)
open(f"{SP}/esm_full_done.txt","w").write("done")
