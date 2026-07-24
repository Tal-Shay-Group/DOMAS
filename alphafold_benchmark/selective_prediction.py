"""(b) Uncertainty routing. The fold_change_prob accuracy-coverage tradeoff and the
profile of the 'uncertain' subset -- the refolding-ambiguous SRP9 class that only
actual folding resolves. Route those to folding instead of forcing a call.
Gene-grouped out-of-fold predictions of the shipped fold-change model."""
import warnings; warnings.filterwarnings("ignore")
import re, gzip, numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"
DATA = "/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-db/data"

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
bf = pd.read_csv(f"{SP}/full_features.csv"); br = pd.read_csv(f"{SP}/bench_rich_results.csv")[['iso', 'region_am', 'max_cov_loss']]
d = bf.merge(br, on='iso').dropna(subset=['mean_rsa']).copy()
d['loeuf'] = d['acc'].map(lambda a: g2l.get(acc2gene.get(a)))
F = ['ident', 'buried_frac', 'mean_rsa', 'region_am', 'loeuf', 'max_cov_loss']
y = (d['tm'] < 0.5).astype(int).values; grp = d['acc'].values

pipe = make_pipeline(SimpleImputer(strategy='median'), StandardScaler(), LogisticRegression(max_iter=3000))
oof = np.full(len(d), np.nan)
for tr, te in GroupKFold(5).split(d[F], y, grp):
    pipe.fit(d[F].iloc[tr], y[tr]); oof[te] = pipe.predict_proba(d[F].iloc[te])[:, 1]
d['p'] = oof; d['conf'] = np.abs(oof - 0.5)
print(f"n={len(d)}  base(TM<0.5)={y.mean():.3f}  overall AUC={roc_auc_score(y, oof):.3f}\n")

print("SELECTIVE PREDICTION (call most-confident X%, route rest to folding):")
order = d.sort_values('conf', ascending=False).reset_index(drop=True)
for cov in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]:
    k = int(cov * len(order)); sub = order.iloc[:k]
    acc = ((sub['p'] >= 0.5).astype(int) == (sub['tm'] < 0.5).astype(int)).mean()
    print(f"  coverage {cov:>4.0%}  n={k:>6}  accuracy={acc:.1%}")

unc = d[(d['p'] >= 0.4) & (d['p'] <= 0.6)]; conf = d[(d['p'] < 0.2) | (d['p'] > 0.8)]
print(f"\nUNCERTAIN band (0.4-0.6): {len(unc)} ({len(unc)/len(d):.0%}), observed change rate={(unc['tm']<0.5).mean():.2f}")
print(f"CONFIDENT band (<0.2 or >0.8): {len(conf)} ({len(conf)/len(d):.0%})")
print("feature profile (uncertain vs confident mean):")
for c in F:
    print(f"  {c:14} uncertain={unc[c].mean():.2f}  confident={conf[c].mean():.2f}")
