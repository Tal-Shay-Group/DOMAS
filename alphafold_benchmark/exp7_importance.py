"""Feature importance for the combined structural model (baseline 6 + PAE + position).
Three views, gene-grouped CV throughout:
  (1) standardized logistic coefficients  (signed effect size, features z-scored)
  (2) drop-column importance               (ΔAUC / ΔR2 when each feature is removed)
  (3) single-feature AUC                    (each feature alone)"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from exp_common import load_master, cv_metrics, fmt, STRUCT, SP

d = load_master()
d = d[d['kind'] != 'insertion'].dropna(subset=['canon_lo']).reset_index(drop=True)
pf = pd.read_csv(f"{SP}/pae_features_full.csv")
d = d.merge(pf, on='iso', how='inner').reset_index(drop=True)
d['dist_term'] = np.minimum(d['canon_lo'] - 1, d['L'] - d['canon_hi']) / d['L']

F = STRUCT + ['pae_global', 'pae_reg2rest', 'protL', 'dist_term']
print(f"combined structural model: n={len(d)} pairs, {d['acc'].nunique()} proteins")
print("features:", F)
full = cv_metrics(d, F, 'y_tm', 'tm')
print(f"\nFULL combined model:  {fmt(full)}\n")

# (1) standardized logistic coefficients (fit on all data, standardized)
y = d['y_tm'].values.astype(float)
Xdf = d[F].copy()
med = Xdf.median(); Xf = Xdf.fillna(med)
mu, sd = Xf.mean(), Xf.std() + 1e-9
Xs = (Xf - mu) / sd
lr = LogisticRegression(max_iter=5000, C=0.2).fit(Xs.values, y)
coef = dict(zip(F, lr.coef_[0]))

# (2) drop-column importance (gene-grouped CV): ΔAUC, Δacc, ΔR2 vs full
print(f"{'feature':<14}{'std_coef':>10}{'dropΔAUC':>10}{'dropΔacc':>10}{'dropΔR2':>9}{'aloneAUC':>10}{'aloneAcc':>10}")
rows = []
for f in F:
    rest = [x for x in F if x != f]
    m_drop = cv_metrics(d, rest, 'y_tm', 'tm')
    m_alone = cv_metrics(d, [f], 'y_tm', 'tm')
    rows.append((f, coef[f], full['auc'] - m_drop['auc'], full['acc'] - m_drop['acc'],
                 full['r2'] - m_drop['r2'], m_alone['auc'], m_alone['acc']))
# sort by drop-AUC importance (largest loss when removed = most important)
for f, c, dauc, dacc, dr2, aA, aAcc in sorted(rows, key=lambda r: -r[2]):
    print(f"{f:<14}{c:>+10.3f}{dauc:>+10.3f}{dacc:>+10.3f}{dr2:>+9.3f}{aA:>10.3f}{aAcc:>10.3f}")

print("\nReading: std_coef = signed effect size (z-scored); dropΔAUC/dropΔacc/dropΔR2 = how")
print("much the metric FALLS when that feature is removed from the full model (bigger =")
print("more important / less redundant); aloneAUC/aloneAcc = that feature by itself.")
print(f"(majority-class accuracy baseline = {max(d['y_tm'].mean(), 1-d['y_tm'].mean()):.3f})")
