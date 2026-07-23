import warnings; warnings.filterwarnings("ignore")
import re, gzip, collections, json, numpy as np, pandas as pd
SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"
DATA = "/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-db/data"

def load_fasta(p):
    out = {}; op = gzip.open if p.endswith('.gz') else open
    with op(p, 'rt') as f:
        cur = None; buf = []
        for ln in f:
            if ln.startswith('>'):
                if cur: out[cur] = ''.join(buf)
                h = ln[1:]; cur = h.split('|')[1] if '|' in h else h.split()[0]; buf = []
            else: buf.append(ln.strip())
        if cur: out[cur] = ''.join(buf)
    return out
canon = load_fasta(f"{DATA}/uniprot/UP000005640_9606.fasta.gz")
iso = load_fasta(f"{DATA}/uniprot/uniprot_sprot_varsplic.fasta.gz")
acc2gene = {}
with gzip.open(f"{DATA}/uniprot/UP000005640_9606.fasta.gz", "rt") as fh:
    for ln in fh:
        if ln.startswith('>'):
            m = re.search(r' GN=(\S+)', ln)
            if m: acc2gene[ln.split('|')[1]] = m.group(1)
g2l = {}
with gzip.open(f"{SP}/gnomad_constraint.txt.bgz", "rt") as fh:
    h = fh.readline().rstrip().split('\t'); gi, li = h.index('gene'), h.index('oe_lof_upper')
    for ln in fh:
        c = ln.rstrip('\n').split('\t')
        try:
            if c[li] not in ('', 'NA'): g2l[c[gi]] = float(c[li])
        except: pass

bf = pd.read_csv(f"{SP}/full_features.csv")            # iso,acc,tm,ident,kind,mean_rsa,buried_frac,...
br = pd.read_csv(f"{SP}/bench_rich_results.csv")[['iso', 'region_am', 'max_cov_loss']]
fp = pd.read_csv(f"{SP}/full_pairs.csv").rename(columns={'isoform': 'iso'})
d = bf.merge(br, on='iso').merge(fp[['iso', 'canon_lo', 'canon_hi']], on='iso').dropna(subset=['mean_rsa', 'canon_lo']).copy()
d['canon_lo'] = d['canon_lo'].astype(int); d['canon_hi'] = d['canon_hi'].astype(int)

# runtime-computable sequence identity from the changed region (shared prefix+suffix / longer seq)
def ident_rt(r):
    cs = canon.get(r['acc']); a = iso.get(r['iso'])
    if not cs or not a: return np.nan
    shared = (r['canon_lo'] - 1) + (len(cs) - r['canon_hi'])
    return 100.0 * shared / max(len(cs), len(a))
d['ident_rt'] = d.apply(ident_rt, axis=1)
d['loeuf'] = d['acc'].map(lambda a: g2l.get(acc2gene.get(a)))
d['y'] = (d['tm'] < 0.5).astype(int)
print(f"identity: runtime vs paper Pearson = {d[['ident_rt','ident']].dropna().corr().iloc[0,1]:.3f}")

FEATURES = ['ident_rt', 'buried_frac', 'mean_rsa', 'region_am', 'loeuf', 'max_cov_loss']
d = d.dropna(subset=['ident_rt']).reset_index(drop=True)
med = {c: d[c].median() for c in FEATURES}
X = d[FEATURES].fillna(med).values.astype(float); y = d['y'].values.astype(float); grp = d['acc'].values
mu, sd = X.mean(0), X.std(0) + 1e-9; Xstd = (X - mu) / sd
def fit(Xs, yv, it=3000, lr=0.3, l2=1e-3):
    w = np.zeros(Xs.shape[1])
    for _ in range(it):
        p = 1 / (1 + np.exp(-np.clip(Xs @ w, -30, 30))); g = Xs.T @ (p - yv) / len(yv) + l2 * w; g[0] = (Xs[:, 0] @ (p - yv)) / len(yv); w -= lr * g
    return w
def auc(yv, s):
    r = pd.Series(s).rank().values; np_ = yv.sum(); nn = len(yv) - np_; return (r[yv == 1].sum() - np_ * (np_ + 1) / 2) / (np_ * nn)
rng = np.random.default_rng(0); uacc = pd.unique(grp); gf = {a: i for a, i in zip(uacc, rng.integers(0, 5, len(uacc)))}
gass = np.array([gf[a] for a in grp]); oof = np.zeros(len(y))
for i in range(5):
    te = np.where(gass == i)[0]; tr = np.where(gass != i)[0]
    w = fit(np.column_stack([np.ones(len(tr)), Xstd[tr]]), y[tr])
    oof[te] = 1 / (1 + np.exp(-np.clip(np.column_stack([np.ones(len(te)), Xstd[te]]) @ w, -30, 30)))
print(f"\nn={len(y)}  base(TM<0.5)={y.mean():.3f}  majority-acc={max(y.mean(),1-y.mean()):.3f}")
print(f"gene-grouped CV AUC = {auc(y, oof):.3f}")

# confusion matrix @ 0.5, and accuracy/sens/spec
for thr in [0.5]:
    pred = (oof >= thr).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum()); fp = int(((pred == 1) & (y == 0)).sum())
    tn = int(((pred == 0) & (y == 0)).sum()); fn = int(((pred == 0) & (y == 1)).sum())
    acc = (tp + tn) / len(y); sens = tp / (tp + fn); spec = tn / (tn + fp); prec = tp / (tp + fp) if tp+fp else 0
    print(f"\nconfusion matrix @ P>={thr} (predict structural change):")
    print(f"                    actual CHANGE   actual PRESERVED")
    print(f"  pred CHANGE          TP={tp:5}          FP={fp:5}")
    print(f"  pred PRESERVED       FN={fn:5}          TN={tn:5}")
    print(f"  accuracy={acc:.3f}  sensitivity={sens:.3f}  specificity={spec:.3f}  precision={prec:.3f}  balanced-acc={(sens+spec)/2:.3f}")

# calibration
print("\ncalibration (OOF pred vs observed):")
dd = d.assign(p=oof); dd['b'] = pd.qcut(dd['p'], 10, labels=False, duplicates='drop')
for b, g in dd.groupby('b'): print(f"  decile {int(b)+1:2}  pred={g['p'].mean():.2f}  observed={g['y'].mean():.2f}")

wf = fit(np.column_stack([np.ones(len(y)), Xstd]), y)
model = {'features': FEATURES, 'median': [round(med[c], 4) for c in FEATURES],
         'mean': [round(float(x), 4) for x in mu], 'std': [round(float(x), 4) for x in sd],
         'intercept': round(float(wf[0]), 4), 'coef': [round(float(x), 4) for x in wf[1:]]}
open(f"{SP}/foldchange_model.json", "w").write(json.dumps(model, indent=2))
print("\nMODEL =", json.dumps(model))
print("standardized coefficients:")
for c, wv in zip(FEATURES, wf[1:]): print(f"  {c:14} {wv:+.3f}")
