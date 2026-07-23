import warnings; warnings.filterwarnings("ignore")
import re, gzip, collections, sqlite3, numpy as np, pandas as pd, sys
sys.path.insert(0, "/Users/arielmelchior/Documents/projects/DOMAS/code")
from enrichment import Enricher
SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"
DATA = "/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-db/data"

# humsavar pathogenic / benign positions
patho = collections.defaultdict(set); benign = collections.defaultdict(set)
rx = re.compile(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}'); st = False
for ln in open(f"{SP}/humsavar.txt"):
    if ln.startswith('_______'): st = True; continue
    if not st: continue
    f = ln.split()
    if len(f) < 5 or not rx.search(f[3]): continue
    (patho if f[4]=='LP/P' else benign if f[4]=='LB/B' else collections.defaultdict(set))[f[1]].add(int(rx.search(f[3]).group(1)))

# acc -> gene symbol from reference-proteome FASTA headers (GN=)
acc2gene = {}
with gzip.open(f"{DATA}/uniprot/UP000005640_9606.fasta.gz", "rt") as fh:
    for ln in fh:
        if ln.startswith('>'):
            acc = ln.split('|')[1]; m = re.search(r' GN=(\S+)', ln)
            if m: acc2gene[acc] = m.group(1)

# gene -> LOEUF (oe_lof_upper) + pLI from gnomAD
gene2loeuf, gene2pli = {}, {}
with gzip.open(f"{SP}/gnomad_constraint.txt.bgz", "rt") as fh:
    hdr = fh.readline().rstrip('\n').split('\t'); gi = hdr.index('gene'); li = hdr.index('oe_lof_upper'); pi = hdr.index('pLI')
    for ln in fh:
        c = ln.rstrip('\n').split('\t')
        try:
            if c[li] not in ('', 'NA'): gene2loeuf[c[gi]] = float(c[li])
            if c[pi] not in ('', 'NA'): gene2pli[c[gi]] = float(c[pi])
        except (ValueError, IndexError): pass
print(f"acc->gene {len(acc2gene)}, LOEUF genes {len(gene2loeuf)}")

# per-pair features
d = pd.read_csv(f"{SP}/bench_rich_results.csv")[['iso','acc','rank','region_am']]
fp = pd.read_csv(f"{SP}/full_pairs.csv").rename(columns={'isoform':'iso'})
d = d.merge(fp[['iso','canon_lo','canon_hi','kind']], on='iso')
d = d[d['kind']!='insertion'].dropna(subset=['canon_lo']).copy()
d['canon_lo']=d['canon_lo'].astype(int); d['canon_hi']=d['canon_hi'].astype(int)
d['reglen']=d['canon_hi']-d['canon_lo']+1

e = Enricher(f"{DATA}/enrichment.sqlite", f"{DATA}/pfam/Pfam-A.hmm", sqlite3.connect("/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-web/DB_merged.sqlite"))
amc = {}
def peak_am(acc, lo, hi):
    if acc not in amc: amc[acc] = e.am_pathogenicity(acc)
    arr = amc[acc]
    if not arr or len(arr) < hi: return np.nan
    seg = [v for v in arr[lo-1:hi] if v is not None]
    return max(seg) if seg else np.nan
d['peak_am'] = [peak_am(a,lo,hi) for a,lo,hi in zip(d['acc'],d['canon_lo'],d['canon_hi'])]
d['mean_am'] = d['region_am']
d['loeuf'] = d['acc'].map(lambda a: gene2loeuf.get(acc2gene.get(a)))
d['pli']   = d['acc'].map(lambda a: gene2pli.get(acc2gene.get(a)))

ov = lambda a,lo,hi,t: bool(t.get(a)) and any(lo<=p<=hi for p in t[a])
d['pat'] = [ov(a,lo,hi,patho) for a,lo,hi in zip(d['acc'],d['canon_lo'],d['canon_hi'])]
d['ben'] = [ov(a,lo,hi,benign) for a,lo,hi in zip(d['acc'],d['canon_lo'],d['canon_hi'])]

# discrimination task: among variant-overlapping regions, pathogenic vs benign-only
disc = d[(d['pat']) | (d['ben'])].copy()
disc['y'] = disc['pat'].astype(int)
disc = disc[~((disc['ben']) & (disc['pat']))].copy() if False else disc  # keep both-overlap in pos
print(f"discrimination set: {len(disc)}  pos(patho)={disc['y'].sum()}  neg(benign-only)={(~disc['pat']&disc['ben']).sum()}")
print(f"  feature coverage: peak_am {disc['peak_am'].notna().mean():.0%}  loeuf {disc['loeuf'].notna().mean():.0%}")

def auc(y,s):
    r=pd.Series(s).rank().values; np_=y.sum(); nn=len(y)-np_; return (r[y==1].sum()-np_*(np_+1)/2)/(np_*nn)
def cv_auc(df, cols, k=5):
    df = df.dropna(subset=cols+['y']); y=df['y'].values.astype(float); n=len(df)
    if n<50 or y.sum()<5 or y.sum()>n-5: return float('nan'), n
    rng=np.random.default_rng(0); idx=rng.permutation(n); folds=np.array_split(idx,k); oof=np.zeros(n); X=df[cols].values.astype(float)
    for i in range(k):
        te=folds[i]; tr=np.concatenate([folds[j] for j in range(k) if j!=i])
        mu,sd=X[tr].mean(0),X[tr].std(0)+1e-9
        Xtr=np.column_stack([np.ones(len(tr)),(X[tr]-mu)/sd]); Xte=np.column_stack([np.ones(len(te)),(X[te]-mu)/sd])
        w=np.zeros(Xtr.shape[1])
        for _ in range(800):
            p=1/(1+np.exp(-np.clip(Xtr@w,-30,30))); g=Xtr.T@(p-y[tr])/len(tr)+1e-3*w; g[0]=(Xtr[:,0]@(p-y[tr]))/len(tr); w-=0.5*g
        oof[te]=1/(1+np.exp(-np.clip(Xte@w,-30,30)))
    return auc(y[:], oof), n

print("\n=== patho-vs-benign discrimination, 5-fold CV AUC (higher = better) ===")
for name, cols in [
    ("region length only",        ['reglen']),
    ("DOMAS impact rank (current)",['rank']),
    ("mean AM (current region_am)",['mean_am']),
    ("PEAK AM (new)",             ['peak_am']),
    ("gnomAD LOEUF (new)",        ['loeuf']),
    ("gnomAD pLI (new)",          ['pli']),
    ("rank + mean_am [baseline]", ['rank','mean_am']),
    ("rank + PEAK_am",            ['rank','peak_am']),
    ("rank + PEAK_am + LOEUF",    ['rank','peak_am','loeuf']),
    ("rank + PEAK_am + LOEUF + pLI",['rank','peak_am','loeuf','pli']),
]:
    a, n = cv_auc(disc, cols)
    print(f"  {name:32} AUC={a:.3f}  (n={n})")
