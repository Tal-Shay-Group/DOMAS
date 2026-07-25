"""Final consolidated importance: all 7 existing model features (union of both shipped
models) + region_plddt = 8 features, evaluated on BOTH targets, gene-grouped CV.
For each: standardized logistic coef + drop-column ΔAUC/Δacc/ΔR2. Concludes which
features belong in which model."""
import warnings; warnings.filterwarnings("ignore")
import sqlite3, numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from exp_common import load_master, load_humsavar, add_overlap_label, cv_metrics, fmt, SP

ENR = "/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-db/data/enrichment.sqlite"
d = load_master()
d = d[d['kind'] != 'insertion'].dropna(subset=['canon_lo']).reset_index(drop=True)
d['canon_lo'] = d['canon_lo'].astype(int); d['canon_hi'] = d['canon_hi'].astype(int)
d = d.drop(columns=['L']).merge(pd.read_csv(f"{SP}/pae_features_full.csv"), on='iso', how='inner').reset_index(drop=True)
con = sqlite3.connect(ENR); cache = {}
def rplddt(acc, lo, hi):
    if acc not in cache:
        r = con.execute("SELECT plddt FROM afdb_plddt WHERE accession=?", (acc,)).fetchone()
        cache[acc] = [float(x) for x in r[0].split(',')] if r else None
    a = cache[acc]
    if not a or len(a) < hi: return np.nan
    seg = a[lo-1:hi]; return round(sum(seg)/len(seg), 2) if seg else np.nan
d['region_plddt'] = [rplddt(a, lo, hi) for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]
d = d.rename(columns={'ident_rt': 'identity'})

FEATS = ['pae_global', 'identity', 'max_cov_loss', 'protL', 'region_am', 'loeuf', 'buried_frac', 'region_plddt']

def std_coefs(df, target):
    X = df[FEATS].copy(); med = X.median(); Xf = X.fillna(med)
    mu, sd = Xf.mean(), Xf.std() + 1e-9
    lr = LogisticRegression(max_iter=5000, C=0.2).fit(((Xf - mu) / sd).values, df[target].values.astype(float))
    return dict(zip(FEATS, lr.coef_[0]))

def report(df, bin_t, cont_t, title):
    full = cv_metrics(df, FEATS, bin_t, cont_t)
    coef = std_coefs(df, bin_t)
    print(f"\n===== {title} =====  n={len(df)}  FULL: {fmt(full)}")
    print(f"{'feature':<14}{'std_coef':>10}{'dropDAUC':>10}{'dropDacc':>10}{'dropDR2':>9}{'aloneAUC':>10}")
    rows = []
    for f in FEATS:
        m = cv_metrics(df, [x for x in FEATS if x != f], bin_t, cont_t)
        a = cv_metrics(df, [f], bin_t, cont_t)
        rows.append((f, coef[f], full['auc']-m['auc'], full['acc']-m['acc'], full['r2']-m['r2'], a['auc']))
    for f, c, dA, dAc, dR, al in sorted(rows, key=lambda r: -r[2]):
        print(f"{f:<14}{c:>+10.2f}{dA:>+10.3f}{dAc:>+10.3f}{dR:>+9.3f}{al:>10.3f}")

# structural: all pairs
report(d, 'y_tm', 'tm', "STRUCTURAL (TM<0.5; R2=continuous TM)")
# functional: disease proteins, pathogenic overlap
patho, _ = load_humsavar()
report(add_overlap_label(d, patho), 'y', None, "FUNCTIONAL (pathogenic overlap; R2=label-variance)")
