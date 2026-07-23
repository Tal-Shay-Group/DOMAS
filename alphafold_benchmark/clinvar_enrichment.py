"""Functional-axis test: do DOMAS impact / AlphaFold burial predict whether an
isoform's changed region overlaps a disease variant?  Label from UniProt
humsavar (pathogenic LP/P vs benign LB/B), annotated in canonical coordinates.
Benign is the specificity control - a real functional score should enrich for
pathogenic overlap but NOT benign. Region length is the dominant confound and is
controlled (within-band tables + as a covariate)."""
import warnings; warnings.filterwarnings("ignore")
import re, collections, numpy as np, pandas as pd
SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"

# --- humsavar: acc -> canonical positions, by clinical significance ---
patho = collections.defaultdict(set); benign = collections.defaultdict(set)
rx = re.compile(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}'); started = False
for ln in open(f"{SP}/humsavar.txt"):
    if ln.startswith('_______'): started = True; continue
    if not started: continue
    f = ln.split()
    if len(f) < 5: continue
    m = rx.search(f[3])
    if not m: continue
    if f[4] == 'LP/P': patho[f[1]].add(int(m.group(1)))
    elif f[4] == 'LB/B': benign[f[1]].add(int(m.group(1)))
print(f"humsavar: LP/P in {len(patho)} proteins, LB/B in {len(benign)} proteins")

# --- pair-level data: changed-region span + burial + impact rank ---
fp = pd.read_csv(f"{SP}/full_pairs.csv").rename(columns={'isoform': 'iso'})
ff = pd.read_csv(f"{SP}/full_features.csv")[['iso', 'buried_frac']]
br = pd.read_csv(f"{SP}/bench_rich_results.csv")[['iso', 'rank']]
d = fp.merge(ff, on='iso').merge(br, on='iso')
d = d[d['kind'] != 'insertion'].dropna(subset=['canon_lo']).copy()
d['canon_lo'] = d['canon_lo'].astype(int); d['canon_hi'] = d['canon_hi'].astype(int)
d['reglen'] = d['canon_hi'] - d['canon_lo'] + 1

def overlap(acc, lo, hi, table):
    s = table.get(acc); return 1 if s and any(lo <= p <= hi for p in s) else 0
def auc(y, s):
    r = pd.Series(s).rank().values; np_ = y.sum(); nn = len(y) - np_
    return (r[y == 1].sum() - np_ * (np_ + 1) / 2) / (np_ * nn)
def cv_auc(df, cols, k=5):
    df = df.dropna(subset=cols + ['y'])                      # drop rows missing any feature (e.g. no structure)
    y = df['y'].values.astype(float); n = len(df); rng = np.random.default_rng(0)
    idx = rng.permutation(n); folds = np.array_split(idx, k); oof = np.zeros(n)
    X = df[cols].values.astype(float)
    for i in range(k):
        te = folds[i]; tr = np.concatenate([folds[j] for j in range(k) if j != i])
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
        Xtr = np.column_stack([np.ones(len(tr)), (X[tr] - mu) / sd]); Xte = np.column_stack([np.ones(len(te)), (X[te] - mu) / sd])
        w = np.zeros(Xtr.shape[1])
        for _ in range(800):
            p = 1 / (1 + np.exp(-np.clip(Xtr @ w, -30, 30)))
            g = Xtr.T @ (p - y[tr]) / len(tr) + 1e-3 * w; g[0] = (Xtr[:, 0] @ (p - y[tr])) / len(tr); w -= 0.5 * g
        oof[te] = 1 / (1 + np.exp(-np.clip(Xte @ w, -30, 30)))
    return auc(y, oof)

for name, tab in [("PATHOGENIC (LP/P)", patho), ("BENIGN (LB/B) - specificity control", benign)]:
    sub = d[d['acc'].isin(tab)].copy()
    sub['y'] = [overlap(a, lo, hi, tab) for a, lo, hi in zip(sub['acc'], sub['canon_lo'], sub['canon_hi'])]
    print(f"\n===== {name} =====  pairs={len(sub)}  proteins={sub['acc'].nunique()}  positives={sub['y'].sum()} ({100*sub['y'].mean():.1f}%)")
    print("  overlap rate by impact:  " + "  ".join(
        f"{lv}={sub[sub['rank']==r]['y'].mean():.3f}(n={len(sub[sub['rank']==r])})" for r, lv in [(0,'none'),(1,'low'),(2,'mod'),(3,'high')]))
    print("  within 1-50aa band:      " + "  ".join(
        f"{lv}={sub[(sub['rank']==r)&(sub['reglen']<50)]['y'].mean():.2f}(n={len(sub[(sub['rank']==r)&(sub['reglen']<50)])})"
        for r, lv in [(0,'none'),(2,'mod'),(3,'high')]))
    b0 = cv_auc(sub, ['reglen', 'identity'])
    print(f"  CV AUC  length={cv_auc(sub,['reglen']):.3f}  base(len+id)={b0:.3f}"
          f"  +impact={cv_auc(sub,['reglen','identity','rank']):.3f} (Δ{cv_auc(sub,['reglen','identity','rank'])-b0:+.3f})"
          f"  +burial={cv_auc(sub,['reglen','identity','buried_frac']):.3f} (Δ{cv_auc(sub,['reglen','identity','buried_frac'])-b0:+.3f})")
