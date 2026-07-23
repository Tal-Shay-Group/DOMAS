import warnings, csv, os, time
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import pydssp

SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"
MAXASA = dict(ALA=129, ARG=274, ASN=195, ASP=193, CYS=167, GLN=225, GLU=223, GLY=104, HIS=224,
              ILE=197, LEU=201, LYS=236, MET=224, PHE=240, PRO=159, SER=155, THR=172, TRP=285, TYR=263, VAL=174)
_parser = PDBParser(QUIET=True)
_sr = ShrakeRupley()

def per_residue(acc):
    """Return (rsa_list, ss_string) for the canonical AF structure, or None."""
    path = f"{SP}/afstruct/AF-{acc}.pdb"
    if not os.path.exists(path):
        return None
    try:
        st = _parser.get_structure(acc, path)
        chain = list(st[0].get_chains())[0]
        res = [r for r in chain if r.id[0] == ' ']
        _sr.compute(st, level='R')
        rsa = [min(1.0, r.sasa / MAXASA.get(r.resname, 200)) for r in res]
        coords = []
        for r in res:
            try:
                coords.append([r['N'].coord, r['CA'].coord, r['C'].coord, r['O'].coord])
            except KeyError:
                coords = None; break
        ss = ''.join(pydssp.assign(np.array(coords, dtype=np.float32), out_type='c3')) if coords else None
        return rsa, ss
    except Exception:
        return None

pairs = list(csv.DictReader(open(f"{SP}/sample_pairs.csv")))
cache = {}
t0 = time.time()
recs = []
for i, p in enumerate(pairs):
    acc = p['acc']
    if acc not in cache:
        cache[acc] = per_residue(acc)
    pr = cache[acc]
    rec = {'iso': p['isoform'], 'acc': acc, 'tm': float(p['tm']), 'ident': float(p['identity']), 'kind': p['kind']}
    if pr and p['canon_lo']:
        rsa, ss = pr
        lo, hi = int(p['canon_lo']), int(p['canon_hi'])
        seg = rsa[lo-1:hi]
        if seg:
            rec['mean_rsa'] = float(np.mean(seg))
            rec['buried_frac'] = float(np.mean([x < 0.25 for x in seg]))
            if ss:
                segss = ss[lo-1:hi]
                rec['helix_frac'] = segss.count('H')/len(segss)
                rec['strand_frac'] = segss.count('E')/len(segss)
                rec['structured_frac'] = (segss.count('H')+segss.count('E'))/len(segss)
    recs.append(rec)
    if (i+1) % 500 == 0:
        print(f"  {i+1}/{len(pairs)} at {time.time()-t0:.0f}s", flush=True)

df = pd.DataFrame(recs)
df.to_csv(f"{SP}/sample_features.csv", index=False)
have = df['mean_rsa'].notna().sum()
print(f"\nfeatures computed for {have}/{len(df)} pairs (rest = insertion/no-structure)  in {time.time()-t0:.0f}s\n")

d = df.dropna(subset=['mean_rsa']).copy()
def pear(a, b):
    m = a.notna() & b.notna(); return float(np.corrcoef(a[m], b[m])[0, 1])
def spear(a, b):
    m = a.notna() & b.notna()
    ra = pd.Series(a[m]).rank().values; rb = pd.Series(b[m]).rank().values
    return float(np.corrcoef(ra, rb)[0, 1])
def partial(col, ctrl='ident', y='tm'):
    x = d[ctrl].values.astype(float); a = d[col].values.astype(float); t = d[y].values.astype(float)
    m = ~np.isnan(a)
    x, a, t = x[m], a[m], t[m]
    ra = a - np.polyval(np.polyfit(x, a, 1), x)
    rt = t - np.polyval(np.polyfit(x, t, 1), x)
    return float(np.corrcoef(a, t)[0, 1]), float(np.corrcoef(ra, rt)[0, 1])

print("=== burial / SS features vs TM (n=%d) ===" % len(d))
print("expected: exposed(mean_rsa) POS with TM; buried_frac/structured NEG with TM\n")
for col in ['mean_rsa', 'buried_frac', 'helix_frac', 'strand_frac', 'structured_frac']:
    if col not in d: continue
    raw, part = partial(col)
    print(f"  {col:16} Spearman={spear(d[col], d['tm']):+.3f}  Pearson={pear(d[col], d['tm']):+.3f}   partial|identity: {raw:+.3f} -> {part:+.3f}")
print(f"\n  reference: identity vs TM  Spearman={spear(d['ident'], d['tm']):+.3f}")

print("\n=== within identity bands: Spearman(mean_rsa, TM) and (buried_frac, TM) ===")
for lo, hi in [(40,60),(60,70),(70,80),(80,85),(85,90),(90,95),(95,100.01)]:
    b = d[(d['ident']>=lo) & (d['ident']<hi)]
    if len(b) < 20: continue
    print(f"  [{lo:>3},{hi:>3})  n={len(b):>4}  Spear(rsa,TM)={spear(b['mean_rsa'],b['tm']):+.3f}  Spear(buried,TM)={spear(b['buried_frac'],b['tm']):+.3f}")

print("\n=== burial by TM verdict ===")
chg = d[d['tm'] < 0.5]; pre = d[d['tm'] >= 0.85]
print(f"  CHANGE  (TM<0.5,   n={len(chg)}): mean_rsa={chg['mean_rsa'].mean():.3f}  buried_frac={chg['buried_frac'].mean():.3f}")
print(f"  PRESERVED(TM>=0.85, n={len(pre)}): mean_rsa={pre['mean_rsa'].mean():.3f}  buried_frac={pre['buried_frac'].mean():.3f}")
open(f"{SP}/sample_features_done.txt", 'w').write("done")
