import sys, os, gzip, csv, sqlite3, collections
import numpy as np, pandas as pd
sys.path.insert(0, "/Users/arielmelchior/Documents/projects/DOMAS/code")
import build
from enrichment import Enricher
from utils import hmm_change_impact, insertion_impact

DATA = "/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-db/data"
SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"
PFAM = f"{DATA}/pfam/Pfam-A.hmm"
RANK = {'none': 0, 'gain': 1, 'low': 1, 'moderate': 2, 'high': 3}

def load_fasta(path):
    out = {}; op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        cur = None; buf = []
        for ln in f:
            if ln.startswith('>'):
                if cur: out[cur] = ''.join(buf)
                h = ln[1:]; cur = h.split('|')[1] if '|' in h else h.split()[0]; buf = []
            else: buf.append(ln.strip())
        if cur: out[cur] = ''.join(buf)
    return out

canon = load_fasta(f"{DATA}/uniprot/UP000005640_9606.fasta.gz")
iso = load_fasta(f"{DATA}/uniprot/uniprot_sprot_varsplic.fasta.gz")
rows = []
for r in csv.DictReader(open(f"{SP}/table_s4_all.csv")):
    if not r['tm']: continue
    acc = r['isoform'].split('-')[0]
    cs, as_ = canon.get(acc), iso.get(r['isoform'])
    if cs and as_ and r['identity']:
        rows.append({'iso': r['isoform'], 'acc': acc, 'tm': float(r['tm']),
                     'ident': float(r['identity']), 'cls': r['class'], 'cseq': cs, 'aseq': as_})

matches = collections.defaultdict(list)
for row in build._parse_hmmsearch_domtbl(f"{SP}/bench_full.domtbl"):
    matches[row[0]].append((row[1], row[10], row[12], row[3], row[4]))

db = sqlite3.connect("/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-web/DB_merged.sqlite")
e = Enricher(f"{DATA}/enrichment.sqlite", PFAM, db)
pl_cache, am_cache = {}, {}
def plddt(a):
    if a not in pl_cache: pl_cache[a] = e.plddt(a)
    return pl_cache[a]
def am(a):
    if a not in am_cache: am_cache[a] = e.am_pathogenicity(a)
    return am_cache[a]

def score(r):
    acc, iso_id = r['acc'], r['iso']
    cm = {m[0]: m for m in matches.get(acc, [])}
    am_m = {m[0]: m for m in matches.get(iso_id, [])}
    reg = Enricher.changed_region(r['cseq'], r['aseq'])
    fold = None; sites = []; ram = None; added = 'none'; net = 0
    if reg and reg['canon']:
        lo, hi = reg['canon']; fold = e.fold_state(acc, lo, hi); sites = e.functional_sites(acc, lo, hi)
        arr = am(acc)
        if arr:
            seg = [v for v in arr[lo-1:hi] if v is not None]; ram = round(sum(seg)/len(seg), 3) if seg else None
    if reg and reg['alt']:
        alo, ahi = reg['alt']; alen = ahi-alo+1
        clen = (reg['canon'][1]-reg['canon'][0]+1) if reg['canon'] else 0
        net = alen-clen
        if net > 0:
            indom = any(min(ahi, m[4])-max(alo, m[3])+1 >= 0.5*alen for m in matches.get(iso_id, []))
            if reg['canon']: flo, fhi = reg['canon']
            else:
                cl = len(r['cseq']); fl = alo-1; flo, fhi = max(1, fl), min(cl, fl+1)
            jf = e.fold_state(acc, flo, fhi) if fhi >= flo else None
            added = insertion_impact(net, indom, jf)
    hits = bool(sites)
    best = (2, 'moderate') if hits else (0, 'none')
    if RANK[added] > best[0]: best = (RANK[added], added)
    pl = plddt(acc); max_cov_loss = 0.0
    for a in set(cm) | set(am_m):
        ccov = cm[a][2] if a in cm else None
        acov = am_m[a][2] if a in am_m else None
        if ccov is not None:
            max_cov_loss = max(max_cov_loss, ccov - (acov or 0))
        cpl = None
        if a in cm and pl:
            s, en = cm[a][3], cm[a][4]; sg = pl[s-1:en] if en <= len(pl) else []
            cpl = round(sum(sg)/len(sg), 1) if sg else None
        lab = hmm_change_impact(ccov, acov, cpl, fold_state=fold, hits_functional_site=hits, region_am=ram)
        if RANK[lab] > best[0]: best = (RANK[lab], lab)
    return {'impact': best[1], 'rank': best[0], 'max_cov_loss': max_cov_loss,
            'region_am': ram if ram is not None else np.nan, 'func_site': int(hits),
            'net_added': max(0, net), 'kind': reg['kind'] if reg else 'identical'}

recs = []
for r in rows:
    s = score(r)
    recs.append({'iso': r['iso'], 'acc': r['acc'], 'tm': r['tm'], 'ident': r['ident'],
                 'cls': r['cls'], **s})
df = pd.DataFrame(recs)
df.to_csv(f"{SP}/bench_rich_results.csv", index=False)
print(f"scored {len(df)} pairs\n")

# constructed continuous DOMAS score (transparent combination of the drivers)
df['domas_score'] = (df['max_cov_loss'].clip(lower=0)
                     + 25*df['func_site']
                     + 50*(df['region_am'].fillna(0.5)-0.5).clip(lower=0)
                     + 0.3*df['net_added'])

def pear(a, b):
    m = a.notna() & b.notna(); a, b = a[m], b[m]
    return float(np.corrcoef(a, b)[0, 1])
def spear(a, b):
    m = a.notna() & b.notna()
    if m.sum() < 3: return float('nan')
    ra = pd.Series(a[m]).rank().values; rb = pd.Series(b[m]).rank().values
    return float(np.corrcoef(ra, rb)[0, 1])

print("=== (B) Correlation of DOMAS numeric scores with TM (expect NEGATIVE) ===")
for name, col in [('impact rank (0-3)', 'rank'), ('max Pfam coverage loss', 'max_cov_loss'),
                  ('region AlphaMissense', 'region_am'), ('constructed domas_score', 'domas_score')]:
    print(f"  {name:26} Spearman={spear(df[col], df['tm']):+.3f}  Pearson={pear(df[col], df['tm']):+.3f}")
print(f"  {'sequence identity':26} Spearman={spear(df['ident'], df['tm']):+.3f}  Pearson={pear(df['ident'], df['tm']):+.3f}  <- the confound")

print("\n=== (A) Control for identity: does score->TM survive within identity bands? ===")
bands = [(40,60),(60,70),(70,80),(80,85),(85,90),(90,95),(95,100.01)]
print(f"  {'identity band':13} {'n':>5} {'Spear(rank,TM)':>14} {'Spear(score,TM)':>15} {'meanTM none':>11} {'meanTM high':>11}")
for lo, hi in bands:
    b = df[(df['ident']>=lo) & (df['ident']<hi)]
    if len(b) < 20: continue
    sr = spear(b['rank'], b['tm']); ss = spear(b['domas_score'], b['tm'])
    tn = b[b['impact']=='none']['tm'].mean(); th = b[b['impact']=='high']['tm'].mean()
    print(f"  [{lo:>3},{hi:>3})     {len(b):>5} {sr:>+14.3f} {ss:>+15.3f} {tn:>11.3f} {th:>11.3f}")

# partial correlation: regress identity out of both (linear), correlate residuals
x = df['ident'].values.astype(float)
for col in ['rank', 'domas_score']:
    y = df[col].values.astype(float); t = df['tm'].values.astype(float)
    by = np.polyfit(x, y, 1); bt = np.polyfit(x, t, 1)
    ry = y - np.polyval(by, x); rt = t - np.polyval(bt, x)
    raw = float(np.corrcoef(y, t)[0,1]); part = float(np.corrcoef(ry, rt)[0,1])
    print(f"  partial Pearson({col},TM | identity): raw={raw:+.3f} -> controlled={part:+.3f}")

print("\n=== (C) Where do the errors fall relative to the TM=0.5 boundary? ===")
POS = {'moderate','high'}
df['pred_change'] = df['impact'].isin(POS)
df['truth_change'] = df['tm'] < 0.5
fp = df[(~df['truth_change']) & (df['pred_change'])]
fn = df[(df['truth_change']) & (~df['pred_change'])]
correct = df[df['pred_change']==df['truth_change']]
print(f"  FP={len(fp)}  FN={len(fn)}  correct={len(correct)}")
print(f"  median |TM-0.5|: correct={ (correct['tm']-0.5).abs().median():.3f}  FP={ (fp['tm']-0.5).abs().median():.3f}  FN={ (fn['tm']-0.5).abs().median():.3f}")
print(f"  FP TM distribution: 0.5-0.6={((fp['tm']>=0.5)&(fp['tm']<0.6)).sum()}  0.6-0.85={((fp['tm']>=0.6)&(fp['tm']<0.85)).sum()}  >=0.85 (egregious)={(fp['tm']>=0.85).sum()}")
print(f"  FN TM distribution: 0.4-0.5={((fn['tm']>=0.4)&(fn['tm']<0.5)).sum()}  0.3-0.4={((fn['tm']>=0.3)&(fn['tm']<0.4)).sum()}  <0.3 (egregious)={(fn['tm']<0.3).sum()}")
# adjacency: for the ordinal levels, how far is predicted rank from a TM-implied level?
def tm_level(t): return 3 if t<0.4 else 2 if t<0.6 else 1 if t<0.85 else 0
df['tm_lvl'] = df['tm'].apply(tm_level)
df['gap'] = (df['rank'] - df['tm_lvl']).abs()
print(f"\n  ordinal gap |DOMAS_rank - TM_level| : ==0 {(df['gap']==0).mean():.0%}  ==1 (adjacent) {(df['gap']==1).mean():.0%}  >=2 (far) {(df['gap']>=2).mean():.0%}")
print(f"  of all disagreements, share that are adjacent-only: {(df[df['gap']>0]['gap']==1).mean():.0%}")
open(f"{SP}/bench_rich_done.txt","w").write("done")
