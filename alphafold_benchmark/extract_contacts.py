import csv, os, time
import numpy as np, pandas as pd

SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"

def cb_coords(acc):
    """Fast PDB parse: per-residue Cβ (Cα for Gly) coord, indexed by residue number.
    Returns dict {resnum: (x,y,z)} for chain A, first altloc only."""
    path = f"{SP}/afstruct/AF-{acc}.pdb"
    if not os.path.exists(path):
        return None
    ca, cb = {}, {}
    try:
        with open(path) as f:
            for ln in f:
                if not ln.startswith('ATOM'):
                    continue
                alt = ln[16]
                if alt not in (' ', 'A'):
                    continue
                atom = ln[12:16].strip()
                if atom not in ('CA', 'CB'):
                    continue
                resnum = int(ln[22:26])
                xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
                (cb if atom == 'CB' else ca)[resnum] = xyz
    except Exception:
        return None
    coord = {r: cb.get(r, ca.get(r)) for r in set(ca) | set(cb)}
    return coord

pairs = list(csv.DictReader(open(f"{SP}/full_pairs.csv")))
cache = {}
recs = []
t0 = time.time()
for i, p in enumerate(pairs):
    acc = p['acc']
    if acc not in cache:
        cache[acc] = cb_coords(acc)
    coord = cache[acc]
    rec = {'iso': p['isoform']}
    if coord and p['canon_lo']:
        lo, hi = int(p['canon_lo']), int(p['canon_hi'])
        region = [r for r in range(lo, hi + 1) if r in coord]
        rest = [r for r in coord if r < lo or r > hi]
        if region and rest:
            R = np.array([coord[r] for r in region])
            A = np.array([coord[r] for r in coord])           # all residues
            region_set = set(region)
            all_res = list(coord)
            # distances region x all
            d = np.sqrt(((R[:, None, :] - A[None, :, :]) ** 2).sum(-1))
            contact = d < 8.0
            is_region = np.array([r in region_set for r in all_res])
            ext = contact[:, ~is_region].sum()
            internal = contact[:, is_region].sum() - len(region)     # minus self-contacts
            rec['ext_contacts_per_res'] = ext / len(region)
            rec['contact_frac_external'] = ext / (ext + internal) if (ext + internal) > 0 else np.nan
    recs.append(rec)
    if (i + 1) % 2000 == 0:
        print(f"  {i+1}/{len(pairs)} at {time.time()-t0:.0f}s", flush=True)

cf = pd.DataFrame(recs)
cf.to_csv(f"{SP}/contact_features.csv", index=False)
print(f"\ncontact features for {cf['ext_contacts_per_res'].notna().sum()}/{len(cf)} in {time.time()-t0:.0f}s\n")

# merge with burial features + labels
bf = pd.read_csv(f"{SP}/full_features.csv")
d = bf.merge(cf, on='iso').dropna(subset=['mean_rsa', 'ext_contacts_per_res'])
print(f"merged n={len(d)}")

def spear(a, b):
    m = a.notna() & b.notna(); return float(np.corrcoef(pd.Series(a[m]).rank(), pd.Series(b[m]).rank())[0, 1])
def resid(y, cols):
    X = np.column_stack([np.ones(len(d))] + [d[c].values.astype(float) for c in cols])
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None); return y - X @ beta
def pcorr(a, cols):
    ra = resid(d[a].values.astype(float), cols); rt = resid(d['tm'].values.astype(float), cols)
    return float(np.corrcoef(ra, rt)[0, 1])
def r2(cols):
    X = np.column_stack([np.ones(len(d))] + [d[c].values.astype(float) for c in cols]); y = d['tm'].values.astype(float)
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None); pred = X @ beta
    return 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()

print("\n=== contact features vs TM ===")
for c in ['ext_contacts_per_res', 'contact_frac_external']:
    print(f"  {c:22} Spearman={spear(d[c], d['tm']):+.3f}  partial|ident={pcorr(c, ['ident']):+.3f}  "
          f"partial|ident+buried={pcorr(c, ['ident', 'buried_frac']):+.3f}")
print(f"\n  corr(ext_contacts, buried_frac) = {spear(d['ext_contacts_per_res'], d['buried_frac']):+.3f}  (expected high - both measure packing)")

print("\n=== incremental R^2 for TM: does contact add beyond identity+burial? ===")
for cols in [['ident'], ['ident', 'buried_frac'], ['ident', 'buried_frac', 'ext_contacts_per_res'],
             ['ident', 'buried_frac', 'contact_frac_external'],
             ['ident', 'buried_frac', 'ext_contacts_per_res', 'contact_frac_external', 'structured_frac']]:
    print(f"  R2 TM ~ {'+'.join(cols):55} = {r2(cols):.3f}")
open(f"{SP}/contact_done.txt", 'w').write("done")
