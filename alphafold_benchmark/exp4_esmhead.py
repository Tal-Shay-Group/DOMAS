"""Task 4: does a pooling that preserves POSITIONAL / shape structure beat mean-pooling
of the same ESM-2 150M per-residue embeddings?  Mean-pool (the shipped approach) vs
multi-pool = concat(mean, max, std, first-token, last-token). If richer pooling of the
same info doesn't help, a heavier attention/CNN head won't either (the E32 ceiling).
3000-row subset, gene-grouped CV, both targets, AUC/acc/R2."""
import warnings; warnings.filterwarnings("ignore")
import gzip, time, numpy as np, pandas as pd, torch, esm, sys
sys.path.insert(0, "/Users/arielmelchior/Documents/projects/DOMAS/code")
from enrichment import Enricher
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, accuracy_score, r2_score
from exp_common import load_fasta, load_humsavar, DATA, SP

canon = load_fasta(f"{DATA}/uniprot/UP000005640_9606.fasta.gz")
iso = load_fasta(f"{DATA}/uniprot/uniprot_sprot_varsplic.fasta.gz")
meta = pd.read_csv(f"{SP}/esm_full_meta.csv")
fp = pd.read_csv(f"{SP}/full_pairs.csv").rename(columns={'isoform': 'iso'})[['iso', 'canon_lo', 'canon_hi']]
meta = meta.merge(fp, on='iso', how='left')
meta['y_tm'] = meta['y']   # meta 'y' is the structural TM<0.5 label
rng = np.random.default_rng(1); idx = np.sort(rng.choice(len(meta), 3000, replace=False))
d = meta.iloc[idx].reset_index(drop=True)
patho, _ = load_humsavar()
d['y_patho'] = [1 if (patho.get(a) and pd.notna(lo) and any(int(lo) <= p <= int(hi) for p in patho[a])) else 0
                for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]
d['is_disease'] = d['acc'].isin(patho) & (d['kind'] != 'insertion')
print(f"subset n={len(d)} proteins={d['acc'].nunique()} disease-rows={int(d['is_disease'].sum())}", flush=True)

model, alphabet = esm.pretrained.esm2_t30_150M_UR50D(); model.eval(); bc = alphabet.get_batch_converter()
def toks_of(seq, lo, hi):
    a = max(0, int(lo) - 1 - 30); b = min(len(seq), int(hi) + 30); win = seq[a:b][:1022]
    if not win: return None
    _, _, t = bc([("x", win)])
    with torch.no_grad(): r = model(t, repr_layers=[30])["representations"][30][0, 1:len(win) + 1]
    return r  # (L,640)
def pools(r):
    if r is None: return np.zeros(640, np.float32), np.zeros(640 * 5, np.float32)
    a = r.numpy()
    mean = a.mean(0)
    multi = np.concatenate([mean, a.max(0), a.std(0), a[0], a[-1]]).astype(np.float32)
    return mean.astype(np.float32), multi

t0 = time.time()
Rc_m = np.zeros((len(d), 640), np.float32); Rc_M = np.zeros((len(d), 640 * 5), np.float32)
Rd_m = np.zeros((len(d), 640), np.float32); Rd_M = np.zeros((len(d), 640 * 5), np.float32)
for i, row in d.iterrows():
    cs = canon.get(row['acc']); a = iso.get(row['iso'])
    reg = Enricher.changed_region(cs, a) if (cs and a) else None
    tc = toks_of(cs, *reg['canon']) if (reg and reg['canon']) else None
    ta = toks_of(a, *reg['alt']) if (reg and reg['alt']) else None
    cm, cM = pools(tc); Rc_m[i] = cm; Rc_M[i] = cM
    if ta is None: Rd_m[i] = 0; Rd_M[i] = 0
    else:
        am, aM = pools(ta); Rd_m[i] = cm - am; Rd_M[i] = cM - aM
    if (i + 1) % 500 == 0: print(f"  embedded {i+1}/{len(d)} at {time.time()-t0:.0f}s", flush=True)
print(f"embedding done {time.time()-t0:.0f}s", flush=True)
np.save(f"{SP}/Rc_multi.npy", Rc_M); np.save(f"{SP}/Rd_multi.npy", Rd_M)

F = ['ident', 'buried_frac', 'mean_rsa', 'region_am', 'loeuf', 'max_cov_loss']
Xb = d[F].values.astype(float); grp = d['acc'].values
ytm = d['y_tm'].values.astype(float) if 'y_tm' in d else (d['tm'].values < 0.5).astype(float)
tm = d['tm'].values.astype(float)
folds = list(GroupKFold(5).split(Xb, ytm, grp))

def run(emb_c, emb_d, nc, yb, cont=None, sub=None):
    keep = np.ones(len(d), bool) if sub is None else sub
    pcl = np.full(len(d), np.nan); preg = np.full(len(d), np.nan)
    for tr, te in folds:
        tr = tr[keep[tr]]; te2 = te[keep[te]]
        def build(ix):
            parts = [Xb[ix]]
            if emb_c is not None:
                pc = PCA(nc).fit(emb_c[tr]); pd_ = PCA(nc).fit(emb_d[tr])
                parts += [pc.transform(emb_c[ix]), pd_.transform(emb_d[ix])]
            return np.hstack(parts)
        cl = make_pipeline(SimpleImputer(strategy='median'), StandardScaler(), LogisticRegression(max_iter=3000, C=0.2))
        cl.fit(build(tr), yb[tr]); pcl[te2] = cl.predict_proba(build(te2))[:, 1]
        if cont is not None:
            rg = make_pipeline(SimpleImputer(strategy='median'), StandardScaler(), Ridge(alpha=5.0))
            rg.fit(build(tr), cont[tr]); preg[te2] = rg.predict(build(te2))
    m = keep & ~np.isnan(pcl)
    r2 = r2_score(cont[m], preg[m]) if cont is not None else r2_score(yb[m], pcl[m])
    return roc_auc_score(yb[m], pcl[m]), accuracy_score(yb[m], (pcl[m] >= 0.5).astype(int)), r2

print("\n=== STRUCTURAL (TM<0.5; R2 on continuous TM) ===", flush=True)
for name, ec, ed in [("baseline (6 feat)", None, None),
                     ("+ mean-pool region+diff", Rc_m, Rd_m),
                     ("+ MULTI-pool region+diff", Rc_M, Rd_M)]:
    a, ac, r = run(ec, ed, 50, ytm, cont=tm); print(f"  {name:26} AUC={a:.3f}  acc={ac:.3f}  R2={r:.3f}", flush=True)

sub = d['is_disease'].values
yp = d['y_patho'].values.astype(float)
print("\n=== FUNCTIONAL (pathogenic overlap, disease proteins; R2 = label-variance) ===", flush=True)
for name, ec, ed in [("baseline (6 feat)", None, None),
                     ("+ mean-pool region+diff", Rc_m, Rd_m),
                     ("+ MULTI-pool region+diff", Rc_M, Rd_M)]:
    a, ac, r = run(ec, ed, 50, yp, cont=None, sub=sub); print(f"  {name:26} AUC={a:.3f}  acc={ac:.3f}  R2={r:.3f}", flush=True)
open(f"{SP}/exp4_done.txt", "w").write("done")
