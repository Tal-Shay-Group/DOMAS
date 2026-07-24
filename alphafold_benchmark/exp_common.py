"""Shared data assembly + gene-grouped CV for the follow-up experiments
(position, Table-S4 descriptors, multi-task, ESM head, PAE, etc.).
Builds one master dataframe with both targets and reusable CV helpers."""
import warnings; warnings.filterwarnings("ignore")
import re, gzip, collections, numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import roc_auc_score, accuracy_score, r2_score

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

def load_master():
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
    bf = pd.read_csv(f"{SP}/full_features.csv")
    br = pd.read_csv(f"{SP}/bench_rich_results.csv")[['iso', 'region_am', 'max_cov_loss', 'rank']]
    fp = pd.read_csv(f"{SP}/full_pairs.csv").rename(columns={'isoform': 'iso'})
    d = bf.merge(br, on='iso').merge(fp[['iso', 'canon_lo', 'canon_hi']], on='iso').copy()
    d = d.dropna(subset=['mean_rsa', 'canon_lo']).copy()
    d['canon_lo'] = d['canon_lo'].astype(int); d['canon_hi'] = d['canon_hi'].astype(int)
    d['reglen'] = d['canon_hi'] - d['canon_lo'] + 1
    d['L'] = d['acc'].map(lambda a: len(canon.get(a, '')) or np.nan)
    def ident_rt(r):
        cs = canon.get(r['acc']); a = iso.get(r['iso'])
        if not cs or not a: return np.nan
        shared = (r['canon_lo'] - 1) + (len(cs) - r['canon_hi'])
        return 100.0 * shared / max(len(cs), len(a))
    d['ident_rt'] = d.apply(ident_rt, axis=1)
    d['loeuf'] = d['acc'].map(lambda a: g2l.get(acc2gene.get(a)))
    d['y_tm'] = (d['tm'] < 0.5).astype(int)
    d = d.dropna(subset=['ident_rt', 'L']).reset_index(drop=True)
    d.attrs['canon'] = canon; d.attrs['iso'] = iso
    return d

# ---- humsavar pathogenic / benign positions (canonical coords) ----
def load_humsavar():
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
    return patho, benign

def add_overlap_label(d, table):
    """Return subset restricted to proteins in `table`, with y = region overlaps a variant."""
    sub = d[d['acc'].isin(table)].copy()
    sub['y'] = [1 if (table.get(a) and any(lo <= p <= hi for p in table[a])) else 0
                for a, lo, hi in zip(sub['acc'], sub['canon_lo'], sub['canon_hi'])]
    return sub

STRUCT = ['ident_rt', 'buried_frac', 'mean_rsa', 'region_am', 'loeuf', 'max_cov_loss']

def cv_auc(df, features, target, groups='acc', n=5, C=0.2, seed=0, extra=None):
    """Gene-grouped CV AUC. `extra` = dict name->ndarray aligned to df for embedding blocks."""
    df = df.reset_index(drop=True)
    y = df[target].values.astype(float); grp = df[groups].values
    Xb = df[features].values.astype(float) if features else np.empty((len(df), 0))
    oof = np.full(len(df), np.nan)
    for tr, te in GroupKFold(n).split(Xb, y, grp):
        Xtr, Xte = [Xb[tr]], [Xb[te]]
        if extra:
            for arr in extra.values():
                Xtr.append(arr[tr]); Xte.append(arr[te])
        m = make_pipeline(SimpleImputer(strategy='median'), StandardScaler(),
                          LogisticRegression(max_iter=3000, C=C))
        m.fit(np.hstack(Xtr), y[tr]); oof[te] = m.predict_proba(np.hstack(Xte))[:, 1]
    return roc_auc_score(y, oof)

def cv_metrics(df, features, bin_target, cont_target=None, groups='acc', n=5, C=0.2, alpha=5.0, extra=None):
    """Gene-grouped CV. Returns dict(auc, acc, r2).
    auc/acc from a logistic classifier on `bin_target` (acc @0.5).
    r2: if `cont_target` given, from a Ridge regressor on that continuous column
        (the structural TM axis); else label-variance R2 of the classifier prob vs
        the binary label (functional axis has no continuous target)."""
    df = df.reset_index(drop=True)
    yb = df[bin_target].values.astype(float); grp = df[groups].values
    yc = df[cont_target].values.astype(float) if cont_target else None
    Xb = df[features].values.astype(float) if features else np.empty((len(df), 0))
    def stack(idx):
        parts = [Xb[idx]]
        if extra:
            for arr in extra.values(): parts.append(arr[idx])
        return np.hstack(parts)
    pcl = np.full(len(df), np.nan); preg = np.full(len(df), np.nan)
    for tr, te in GroupKFold(n).split(Xb, yb, grp):
        cl = make_pipeline(SimpleImputer(strategy='median'), StandardScaler(),
                           LogisticRegression(max_iter=3000, C=C))
        cl.fit(stack(tr), yb[tr]); pcl[te] = cl.predict_proba(stack(te))[:, 1]
        if cont_target:
            rg = make_pipeline(SimpleImputer(strategy='median'), StandardScaler(), Ridge(alpha=alpha))
            rg.fit(stack(tr), yc[tr]); preg[te] = rg.predict(stack(te))
    auc = roc_auc_score(yb, pcl)
    acc = accuracy_score(yb, (pcl >= 0.5).astype(int))
    r2 = r2_score(yc, preg) if cont_target else r2_score(yb, pcl)
    return {'auc': auc, 'acc': acc, 'r2': r2}

def fmt(m):
    return f"AUC={m['auc']:.3f}  acc={m['acc']:.3f}  R2={m['r2']:.3f}"
