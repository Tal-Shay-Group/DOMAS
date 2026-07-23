import warnings; warnings.filterwarnings("ignore")
import re, gzip, collections, sqlite3, json, numpy as np, pandas as pd, sys
sys.path.insert(0, "/Users/arielmelchior/Documents/projects/DOMAS/code")
from enrichment import Enricher
SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"
DATA = "/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-db/data"

# labels + gene constraint
patho = collections.defaultdict(set); benign = collections.defaultdict(set)
rx = re.compile(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}'); st = False
for ln in open(f"{SP}/humsavar.txt"):
    if ln.startswith('_______'): st = True; continue
    if not st: continue
    f = ln.split()
    if len(f) < 5 or not rx.search(f[3]): continue
    if f[4] == 'LP/P': patho[f[1]].add(int(rx.search(f[3]).group(1)))
    elif f[4] == 'LB/B': benign[f[1]].add(int(rx.search(f[3]).group(1)))
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

# features
d = pd.read_csv(f"{SP}/bench_rich_results.csv")[['iso', 'acc', 'region_am', 'max_cov_loss']]
ff = pd.read_csv(f"{SP}/full_features.csv")[['iso', 'buried_frac']]
fp = pd.read_csv(f"{SP}/full_pairs.csv").rename(columns={'isoform': 'iso'})
d = d.merge(ff, on='iso').merge(fp[['iso', 'canon_lo', 'canon_hi', 'kind']], on='iso')
d = d[d['kind'] != 'insertion'].dropna(subset=['canon_lo']).copy()
d['canon_lo'] = d['canon_lo'].astype(int); d['canon_hi'] = d['canon_hi'].astype(int)
d['loeuf'] = d['acc'].map(lambda a: g2l.get(acc2gene.get(a)))
ov = lambda a, lo, hi, t: bool(t.get(a)) and any(lo <= p <= hi for p in t[a])
d['pat'] = [ov(a, lo, hi, patho) for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]
d['ben'] = [ov(a, lo, hi, benign) for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]

FEATURES = ['region_am', 'loeuf', 'max_cov_loss', 'buried_frac']
disc = d[(d['pat']) | (d['ben'])].copy(); disc['y'] = disc['pat'].astype(int)
# impute missing with training median (stored for runtime); standardize
med = {c: disc[c].median() for c in FEATURES}
X = disc[FEATURES].fillna(med).values.astype(float); y = disc['y'].values.astype(float)
mu, sd = X.mean(0), X.std(0) + 1e-9

def fit(Xs, yv, iters=3000, lr=0.3, l2=1e-3):
    w = np.zeros(Xs.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-np.clip(Xs @ w, -30, 30)))
        g = Xs.T @ (p - yv) / len(yv) + l2 * w; g[0] = (Xs[:, 0] @ (p - yv)) / len(yv)
        w -= lr * g
    return w
def auc(yv, s):
    r = pd.Series(s).rank().values; np_ = yv.sum(); nn = len(yv) - np_
    return (r[yv == 1].sum() - np_ * (np_ + 1) / 2) / (np_ * nn)

# 5-fold CV for honest AUC + out-of-fold probs (for calibration/deciles)
n = len(disc); rng = np.random.default_rng(0); idx = rng.permutation(n); folds = np.array_split(idx, 5); oof = np.zeros(n)
Xstd = (X - mu) / sd
for i in range(5):
    te = folds[i]; tr = np.concatenate([folds[j] for j in range(5) if j != i])
    w = fit(np.column_stack([np.ones(len(tr)), Xstd[tr]]), y[tr])
    oof[te] = 1 / (1 + np.exp(-np.clip(np.column_stack([np.ones(len(te)), Xstd[te]]) @ w, -30, 30)))
print(f"n={n}  positives(patho)={int(y.sum())}  base={y.mean():.3f}")
print(f"5-fold CV AUC (pair-level) = {auc(y, oof):.3f}")

# gene-grouped CV: each protein confined to one fold (no gene in train AND test),
# the leakage-free estimate since LOEUF is gene-level and overlap is gene-correlated
grp = disc['acc'].values; uacc = pd.unique(grp)
gfold = {a: i for a, i in zip(uacc, rng.integers(0, 5, len(uacc)))}
gassign = np.array([gfold[a] for a in grp]); goof = np.zeros(n)
for i in range(5):
    te = np.where(gassign == i)[0]; tr = np.where(gassign != i)[0]
    w = fit(np.column_stack([np.ones(len(tr)), Xstd[tr]]), y[tr])
    goof[te] = 1 / (1 + np.exp(-np.clip(np.column_stack([np.ones(len(te)), Xstd[te]]) @ w, -30, 30)))
print(f"5-fold CV AUC (gene-grouped, leakage-free) = {auc(y, goof):.3f}  ({len(uacc)} proteins)")

# final model on all data (the shipped coefficients)
wf = fit(np.column_stack([np.ones(n), Xstd]), y)
print("\ncalibration (out-of-fold predicted prob vs observed pathogenic rate):")
disc = disc.assign(p=oof); disc['pbin'] = pd.qcut(disc['p'], 10, labels=False, duplicates='drop')
for b, g in disc.groupby('pbin'):
    print(f"  decile {int(b)+1:2}  pred={g['p'].mean():.2f}  observed={g['y'].mean():.2f}  n={len(g)}")

model = {'features': FEATURES, 'median': [round(med[c], 4) for c in FEATURES],
         'mean': [round(float(x), 4) for x in mu], 'std': [round(float(x), 4) for x in sd],
         'intercept': round(float(wf[0]), 4), 'coef': [round(float(x), 4) for x in wf[1:]]}
print("\nMODEL =", json.dumps(model))
open(f"{SP}/calib_model.json", "w").write(json.dumps(model, indent=2))
# feature contributions (standardized coef)
print("\nstandardized coefficients (impact on log-odds):")
for c, wv in zip(FEATURES, wf[1:]): print(f"  {c:14} {wv:+.3f}")
