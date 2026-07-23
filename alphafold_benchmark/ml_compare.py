import warnings; warnings.filterwarnings("ignore")
import re, gzip, collections, numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import (GradientBoostingClassifier, RandomForestClassifier,
                              HistGradientBoostingClassifier)
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score
SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"
DATA = "/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-db/data"

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

bf = pd.read_csv(f"{SP}/full_features.csv")
br = pd.read_csv(f"{SP}/bench_rich_results.csv")[['iso', 'region_am', 'max_cov_loss']]
fp = pd.read_csv(f"{SP}/full_pairs.csv").rename(columns={'isoform': 'iso'})
d = bf.merge(br, on='iso').merge(fp[['iso', 'canon_lo', 'canon_hi']], on='iso').dropna(subset=['mean_rsa']).copy()
d['loeuf'] = d['acc'].map(lambda a: g2l.get(acc2gene.get(a)))
d['canon_lo'] = d['canon_lo'].astype(int); d['canon_hi'] = d['canon_hi'].astype(int)
d['pat'] = [bool(patho.get(a)) and any(lo <= p <= hi for p in patho[a]) for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]
d['ben'] = [bool(benign.get(a)) and any(lo <= p <= hi for p in benign[a]) for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]

def models():
    imp = SimpleImputer(strategy='median')
    return {
        'logistic (current)':      make_pipeline(SimpleImputer(strategy='median'), StandardScaler(), LogisticRegression(max_iter=2000, C=1.0)),
        'logistic + interactions': make_pipeline(SimpleImputer(strategy='median'), StandardScaler(), PolynomialFeatures(2, interaction_only=True, include_bias=False), LogisticRegression(max_iter=3000, C=1.0)),
        'random forest':           make_pipeline(SimpleImputer(strategy='median'), RandomForestClassifier(n_estimators=400, min_samples_leaf=10, n_jobs=-1, random_state=0)),
        'gradient boosting':       make_pipeline(SimpleImputer(strategy='median'), GradientBoostingClassifier(random_state=0)),
        'hist grad boosting':      HistGradientBoostingClassifier(random_state=0),   # handles NaN
        'neural net (MLP 32,16)':  make_pipeline(SimpleImputer(strategy='median'), StandardScaler(), MLPClassifier(hidden_layer_sizes=(32, 16), max_iter=800, alpha=1e-3, random_state=0)),
    }

def evaluate(df, feats, y, groups, name):
    print(f"\n=== {name}  (n={len(df)}, positives={int(y.sum())}, {pd.Series(groups).nunique()} proteins) ===")
    print("  gene-grouped 5-fold CV AUC:")
    gkf = GroupKFold(n_splits=5)
    for mname, model in models().items():
        oof = np.full(len(df), np.nan)
        for tr, te in gkf.split(df[feats], y, groups):
            model.fit(df[feats].iloc[tr], y[tr])
            oof[te] = model.predict_proba(df[feats].iloc[te])[:, 1]
        print(f"    {mname:26} {roc_auc_score(y, oof):.3f}")

# functional target: pathogenic vs benign-only, among variant-overlapping
disc = d[d['pat'] | d['ben']].reset_index(drop=True)
evaluate(disc, ['region_am', 'loeuf', 'max_cov_loss', 'buried_frac'],
         disc['pat'].astype(int).values, disc['acc'].values, "FUNCTIONAL: pathogenic vs benign")

# structural target: TM < 0.5
d2 = d.reset_index(drop=True); y2 = (d2['tm'] < 0.5).astype(int).values
evaluate(d2, ['ident', 'buried_frac', 'mean_rsa', 'region_am', 'loeuf', 'max_cov_loss'],
         y2, d2['acc'].values, "STRUCTURAL: TM<0.5")
