import numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore")
SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"
rng = np.random.default_rng(0)

bf = pd.read_csv(f"{SP}/full_features.csv")
imp = pd.read_csv(f"{SP}/bench_rich_results.csv")[['iso', 'rank']]
d = bf.merge(imp, on='iso', how='left').dropna(subset=['mean_rsa']).reset_index(drop=True)
d['rank'] = d['rank'].fillna(0)
y = (d['tm'].values < 0.5).astype(float)           # CHANGE = positive
tm = d['tm'].values.astype(float)
n = len(d); base = y.mean()
print(f"n={n} (non-insertion pairs with burial)   base rate CHANGE(TM<0.5)={base:.3f}   majority-class acc={max(base,1-base):.3f}\n")

def auc(yv, s):
    r = pd.Series(s).rank().values
    npos = yv.sum(); nneg = len(yv) - npos
    return (r[yv == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg)

def fit_lr(X, yv, iters=800, lr=0.5, l2=1e-3):
    w = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-np.clip(X @ w, -30, 30)))
        g = X.T @ (p - yv) / len(yv) + l2 * w
        g[0] = (X[:, 0] @ (p - yv)) / len(yv)      # no reg on intercept
        w -= lr * g
    return w

def cv_classify(cols, k=5):
    idx = rng.permutation(n); folds = np.array_split(idx, k)
    oof = np.zeros(n)
    for i in range(k):
        te = folds[i]; tr = np.concatenate([folds[j] for j in range(k) if j != i])
        Xtr = d.loc[tr, cols].values.astype(float); Xte = d.loc[te, cols].values.astype(float)
        mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
        Xtr = np.column_stack([np.ones(len(tr)), (Xtr - mu) / sd])
        Xte = np.column_stack([np.ones(len(te)), (Xte - mu) / sd])
        w = fit_lr(Xtr, y[tr])
        oof[te] = 1 / (1 + np.exp(-np.clip(Xte @ w, -30, 30)))
    a = auc(y, oof)
    # accuracy at 0.5, and at Youden-optimal threshold
    thr05 = (oof >= 0.5).astype(float); acc05 = (thr05 == y).mean()
    ts = np.unique(oof); best = (0, 0.5)
    for t in ts[::max(1, len(ts) // 200)]:
        pr = (oof >= t).astype(float)
        sens = pr[y == 1].mean(); spec = 1 - pr[y == 0].mean()
        if sens + spec > best[0]: best = (sens + spec, t)
    t = best[1]; pr = (oof >= t).astype(float)
    sens = pr[y == 1].mean(); spec = 1 - pr[y == 0].mean(); accY = (pr == y).mean()
    return a, acc05, accY, sens, spec

def cv_regress(cols, k=5):
    idx = rng.permutation(n); folds = np.array_split(idx, k); oof = np.zeros(n)
    for i in range(k):
        te = folds[i]; tr = np.concatenate([folds[j] for j in range(k) if j != i])
        Xtr = np.column_stack([np.ones(len(tr))] + [d.loc[tr, c].values.astype(float) for c in cols])
        Xte = np.column_stack([np.ones(len(te))] + [d.loc[te, c].values.astype(float) for c in cols])
        w, *_ = np.linalg.lstsq(Xtr, tm[tr], rcond=None)
        oof[te] = Xte @ w
    r2 = 1 - ((tm - oof) ** 2).sum() / ((tm - tm.mean()) ** 2).sum()
    rmse = np.sqrt(((tm - oof) ** 2).mean())
    return r2, rmse

print("=== CLASSIFICATION: predict CHANGE (TM<0.5), 5-fold CV ===")
print(f"{'model':34} {'AUC':>6} {'acc@.5':>7} {'accY':>6} {'sens':>6} {'spec':>6}")
models = [
    ('DOMAS impact rank (current)', ['rank']),
    ('sequence identity', ['ident']),
    ('burial (buried_frac)', ['buried_frac']),
    ('burial (mean_rsa)', ['mean_rsa']),
    ('identity + burial', ['ident', 'buried_frac']),
    ('identity + burial + SS + rsa', ['ident', 'buried_frac', 'structured_frac', 'mean_rsa']),
    ('identity + impact (control)', ['ident', 'rank']),
]
for name, cols in models:
    a, acc05, accY, sens, spec = cv_classify(cols)
    print(f"{name:34} {a:6.3f} {acc05:7.3f} {accY:6.3f} {sens:6.3f} {spec:6.3f}")

print("\n=== REGRESSION: predict the TM value, 5-fold CV ===")
print(f"{'model':34} {'R2':>7} {'RMSE':>7}")
for name, cols in [('burial (buried_frac) alone', ['buried_frac']),
                   ('mean_rsa alone', ['mean_rsa']),
                   ('sequence identity alone', ['ident']),
                   ('identity + burial', ['ident', 'buried_frac']),
                   ('identity + burial + SS + rsa', ['ident', 'buried_frac', 'structured_frac', 'mean_rsa'])]:
    r2, rmse = cv_regress(cols)
    print(f"{name:34} {r2:7.3f} {rmse:7.3f}")
