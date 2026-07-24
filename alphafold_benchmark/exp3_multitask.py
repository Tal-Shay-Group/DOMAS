"""Task 3: multi-task joint model. Shared trunk + 3 heads (TM-class, pathogenic-class,
TM-regression), trained jointly with masked losses (functional loss only on disease
proteins). Question: does the auxiliary task regularise the shared features and beat
the two INDEPENDENT logistic models?  Gene-grouped 5-fold CV, evaluated per target on
the same held-out folds. R2: structural = continuous TM (reg head); functional =
label-variance."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, accuracy_score, r2_score
from exp_common import load_master, load_humsavar, cv_metrics, fmt, STRUCT

torch.manual_seed(0); np.random.seed(0)
d = load_master()
patho, _ = load_humsavar()
# functional label on all rows; mask = disease protein & non-insertion
d['is_disease'] = d['acc'].isin(patho) & (d['kind'] != 'insertion')
d['y_patho'] = [1 if (patho.get(a) and any(lo <= p <= hi for p in patho[a])) else 0
                for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]
d = d.reset_index(drop=True)
F = STRUCT
X = d[F].values.astype(float)
ytm = d['y_tm'].values.astype(float); tm = d['tm'].values.astype(float)
yp = d['y_patho'].values.astype(float); mask = d['is_disease'].values.astype(float)
grp = d['acc'].values
print(f"n={len(d)}  disease-protein rows (functional eval)={int(mask.sum())}  patho-pos={int(yp[mask==1].sum())}")

class MT(nn.Module):
    def __init__(s, p):
        super().__init__()
        s.trunk = nn.Sequential(nn.Linear(p, 32), nn.ReLU(), nn.Linear(32, 16), nn.ReLU())
        s.h_tm = nn.Linear(16, 1); s.h_pa = nn.Linear(16, 1); s.h_reg = nn.Linear(16, 1)
    def forward(s, x):
        z = s.trunk(x); return s.h_tm(z).squeeze(-1), s.h_pa(z).squeeze(-1), s.h_reg(z).squeeze(-1)

bce = nn.BCEWithLogitsLoss(reduction='none'); mse = nn.MSELoss(reduction='mean')
oof_tm = np.full(len(d), np.nan); oof_pa = np.full(len(d), np.nan); oof_reg = np.full(len(d), np.nan)
for tr, te in GroupKFold(5).split(X, ytm, grp):
    mu = np.nanmedian(X[tr], 0); Xi = np.where(np.isnan(X), mu, X)
    m_, s_ = Xi[tr].mean(0), Xi[tr].std(0) + 1e-9; Xs = (Xi - m_) / s_
    Xtr = torch.tensor(Xs[tr], dtype=torch.float32); Xte = torch.tensor(Xs[te], dtype=torch.float32)
    net = MT(X.shape[1]); opt = torch.optim.Adam(net.parameters(), lr=0.01, weight_decay=1e-3)
    ytm_t = torch.tensor(ytm[tr]); yp_t = torch.tensor(yp[tr]); tm_t = torch.tensor(tm[tr], dtype=torch.float32)
    mk_t = torch.tensor(mask[tr])
    for ep in range(400):
        opt.zero_grad(); lt, lp, lr_ = net(Xtr)
        loss = bce(lt, ytm_t).mean() + (bce(lp, yp_t) * mk_t).sum() / (mk_t.sum() + 1e-9) + mse(lr_, tm_t)
        loss.backward(); opt.step()
    net.eval()
    with torch.no_grad():
        lt, lp, lr_ = net(Xte)
        oof_tm[te] = torch.sigmoid(lt).numpy(); oof_pa[te] = torch.sigmoid(lp).numpy(); oof_reg[te] = lr_.numpy()

print("\n=== STRUCTURAL (TM<0.5) ===")
print(f"  single-task logistic   {fmt(cv_metrics(d, F, 'y_tm', 'tm'))}")
print(f"  multi-task net         AUC={roc_auc_score(ytm,oof_tm):.3f}  acc={accuracy_score(ytm,(oof_tm>=0.5)):.3f}  R2={r2_score(tm,oof_reg):.3f}")

di = d[mask == 1]
print("\n=== FUNCTIONAL (pathogenic overlap, disease proteins) ===")
print(f"  single-task logistic   {fmt(cv_metrics(di, F, 'y_patho'))}")
print(f"  multi-task net         AUC={roc_auc_score(yp[mask==1],oof_pa[mask==1]):.3f}  acc={accuracy_score(yp[mask==1],(oof_pa[mask==1]>=0.5)):.3f}  R2={r2_score(yp[mask==1],oof_pa[mask==1]):.3f}")
