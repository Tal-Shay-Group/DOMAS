import warnings, csv, os, time
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import pydssp

SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"
MAXASA = dict(ALA=129, ARG=274, ASN=195, ASP=193, CYS=167, GLN=225, GLU=223, GLY=104, HIS=224,
              ILE=197, LEU=201, LYS=236, MET=224, PHE=240, PRO=159, SER=155, THR=172, TRP=285, TYR=263, VAL=174)
_parser = PDBParser(QUIET=True); _sr = ShrakeRupley()

def per_residue(acc):
    path = f"{SP}/afstruct/AF-{acc}.pdb"
    if not os.path.exists(path): return None
    try:
        st = _parser.get_structure(acc, path)
        res = [r for r in list(st[0].get_chains())[0] if r.id[0] == ' ']
        _sr.compute(st, level='R')
        rsa = [min(1.0, r.sasa / MAXASA.get(r.resname, 200)) for r in res]
        coords = []
        for r in res:
            try: coords.append([r['N'].coord, r['CA'].coord, r['C'].coord, r['O'].coord])
            except KeyError: coords = None; break
        ss = ''.join(pydssp.assign(np.array(coords, dtype=np.float32), out_type='c3')) if coords else None
        return rsa, ss
    except Exception:
        return None

pairs = list(csv.DictReader(open(f"{SP}/full_pairs.csv")))
cache = {}; recs = []; t0 = time.time()
for i, p in enumerate(pairs):
    acc = p['acc']
    if acc not in cache: cache[acc] = per_residue(acc)
    pr = cache[acc]
    rec = {'iso': p['isoform'], 'acc': acc, 'tm': float(p['tm']), 'ident': float(p['identity']),
           'kind': p['kind'], 'reglen': (int(p['canon_hi'])-int(p['canon_lo'])+1) if p['canon_lo'] else np.nan}
    if pr and p['canon_lo']:
        rsa, ss = pr; lo, hi = int(p['canon_lo']), int(p['canon_hi']); seg = rsa[lo-1:hi]
        if seg:
            rec['mean_rsa'] = float(np.mean(seg)); rec['buried_frac'] = float(np.mean([x < 0.25 for x in seg]))
            if ss:
                s = ss[lo-1:hi]; rec['helix_frac'] = s.count('H')/len(s); rec['strand_frac'] = s.count('E')/len(s)
                rec['structured_frac'] = (s.count('H')+s.count('E'))/len(s)
    recs.append(rec)
    if (i+1) % 2000 == 0: print(f"  {i+1}/{len(pairs)} at {time.time()-t0:.0f}s", flush=True)

df = pd.DataFrame(recs); df.to_csv(f"{SP}/full_features.csv", index=False)
have = df['mean_rsa'].notna().sum()
print(f"\nfeatures for {have}/{len(df)} pairs (rest insertion/no-structure), {len(cache)} structures, {time.time()-t0:.0f}s\n")

d = df.dropna(subset=['mean_rsa']).copy()
def pear(a,b): m=a.notna()&b.notna(); return float(np.corrcoef(a[m],b[m])[0,1])
def spear(a,b):
    m=a.notna()&b.notna(); return float(np.corrcoef(pd.Series(a[m]).rank(),pd.Series(b[m]).rank())[0,1])
def resid(y,cols):
    X=np.column_stack([np.ones(len(d))]+[d[c].values.astype(float) for c in cols])
    beta,_,_,_=np.linalg.lstsq(X,y,rcond=None); return y-X@beta
def pcorr(a,cols):
    ra=resid(d[a].values.astype(float),cols); rt=resid(d['tm'].values.astype(float),cols)
    return float(np.corrcoef(ra,rt)[0,1])
def r2(cols):
    X=np.column_stack([np.ones(len(d))]+[d[c].values.astype(float) for c in cols]); y=d['tm'].values.astype(float)
    beta,_,_,_=np.linalg.lstsq(X,y,rcond=None); pred=X@beta
    return 1-((y-pred)**2).sum()/((y-y.mean())**2).sum()

print(f"=== ALL 11k: burial/SS vs TM (n={len(d)}) ===")
for c in ['mean_rsa','buried_frac','helix_frac','strand_frac','structured_frac']:
    print(f"  {c:16} Spearman={spear(d[c],d['tm']):+.3f}  Pearson={pear(d[c],d['tm']):+.3f}  "
          f"partial|ident={pcorr(c,['ident']):+.3f}  partial|ident+len={pcorr(c,['ident','reglen']):+.3f}")
print(f"  reference identity vs TM Spearman={spear(d['ident'],d['tm']):+.3f}")

print("\n=== within identity bands ===")
for lo,hi in [(40,60),(60,70),(70,80),(80,85),(85,90),(90,95),(95,100.01)]:
    b=d[(d['ident']>=lo)&(d['ident']<hi)]
    if len(b)<20: continue
    print(f"  [{lo:>3},{hi:>3}) n={len(b):>5} Spear(rsa,TM)={spear(b['mean_rsa'],b['tm']):+.3f} Spear(buried,TM)={spear(b['buried_frac'],b['tm']):+.3f}")

print("\n=== incremental variance explained (linear R^2 for TM) ===")
for cols in [['ident'],['ident','reglen'],['ident','buried_frac'],['ident','buried_frac','structured_frac'],['ident','buried_frac','reglen']]:
    print(f"  R2 TM ~ {'+'.join(cols):40} = {r2(cols):.3f}")

chg=d[d['tm']<0.5]; pre=d[d['tm']>=0.85]
print(f"\nCHANGE(TM<0.5,n={len(chg)}): mean_rsa={chg['mean_rsa'].mean():.3f} buried={chg['buried_frac'].mean():.3f}")
print(f"PRESERVED(TM>=0.85,n={len(pre)}): mean_rsa={pre['mean_rsa'].mean():.3f} buried={pre['buried_frac'].mean():.3f}")
open(f"{SP}/full_features_done.txt",'w').write("done")
