"""Task 6: does cross-species conservation (phyloP) of the changed region add signal
BEYOND AlphaMissense (region_am)?  E26 assumed redundant (AM is trained on
conservation) but never measured it. Pipeline: EBI coordinates API (protein->genome)
-> map changed region to genomic intervals -> UCSC hg38 100-way phyloP bigWig (remote)
-> mean phyloP over region and over whole protein. Bounded sample; both targets."""
import warnings; warnings.filterwarnings("ignore")
import urllib.request, json, time, numpy as np, pandas as pd, pyBigWig
from exp_common import load_master, load_humsavar, add_overlap_label, cv_metrics, fmt, STRUCT, SP

BW = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw"
d = load_master()
d = d[d['kind'] != 'insertion'].dropna(subset=['canon_lo']).reset_index(drop=True)
d['canon_lo'] = d['canon_lo'].astype(int); d['canon_hi'] = d['canon_hi'].astype(int)
rng = np.random.default_rng(0)
samp = d.iloc[np.sort(rng.choice(len(d), 400, replace=False))].reset_index(drop=True)

def coords(acc):
    try:
        req = urllib.request.Request(f"https://www.ebi.ac.uk/proteins/api/coordinates/{acc}", headers={"Accept": "application/json"})
        r = json.load(urllib.request.urlopen(req, timeout=30))
        gn = r['gnCoordinate'][0]['genomicLocation']
        exons = []
        for e in gn['exon']:
            pb, pe = e['proteinLocation']['begin']['position'], e['proteinLocation']['end']['position']
            gb, ge = e['genomeLocation']['begin']['position'], e['genomeLocation']['end']['position']
            exons.append((pb, pe, min(gb, ge), max(gb, ge)))
        return {'chrom': 'chr' + str(gn['chromosome']), 'rev': bool(gn['reverseStrand']), 'exons': exons}
    except Exception:
        return None

def region_intervals(cm, lo, hi):
    """genomic [start,end] intervals covering protein positions [lo,hi]."""
    out = []
    for pb, pe, gb, ge in cm['exons']:
        a, b = max(lo, pb), min(hi, pe)
        if a > b: continue
        lo_off, hi_off = a - pb, b - pb   # aa offset from exon 5' (transcription) end
        if not cm['rev']:
            out.append((gb + lo_off * 3, gb + hi_off * 3 + 2))
        else:
            out.append((ge - hi_off * 3 - 2, ge - lo_off * 3))
    return out

cache = {}
bw = pyBigWig.open(BW)
def phylop_mean(chrom, ivs):
    vals = []; tot = 0
    for s, e in ivs:
        try:
            m = bw.stats(chrom, max(0, s - 1), e, type="mean")[0]
            if m is not None: vals.append(m * (e - s + 1)); tot += (e - s + 1)
        except Exception: pass
    return sum(vals) / tot if tot else np.nan

rows = []; t0 = time.time()
for i, r in samp.iterrows():
    acc = r['acc']
    if acc not in cache: cache[acc] = coords(acc)
    cm = cache[acc]
    if not cm: continue
    reg = phylop_mean(cm['chrom'], region_intervals(cm, r['canon_lo'], r['canon_hi']))
    whole = phylop_mean(cm['chrom'], [(gb, ge) for _, _, gb, ge in cm['exons']])
    rows.append({'iso': r['iso'], 'phylop_region': reg, 'phylop_whole': whole,
                 'phylop_rel': (reg - whole) if (reg == reg and whole == whole) else np.nan})
    if (i + 1) % 100 == 0: print(f"  {i+1}/{len(samp)} at {time.time()-t0:.0f}s", flush=True)
bw.close()
pf = pd.DataFrame(rows).dropna(subset=['phylop_region'])
pf.to_csv(f"{SP}/phylop_features.csv", index=False)
d2 = d.merge(pf, on='iso', how='inner').reset_index(drop=True)
print(f"\nphyloP computed for {len(d2)} pairs  proteins={d2['acc'].nunique()}", flush=True)
print(f"corr(phylop_region, region_am) = {d2[['phylop_region','region_am']].corr().iloc[0,1]:.3f}  "
      f"Spearman(phylop_region, TM) = {d2[['phylop_region','tm']].corr('spearman').iloc[0,1]:+.3f}", flush=True)

PH = ['phylop_region', 'phylop_whole', 'phylop_rel']
print("\n=== STRUCTURAL (TM<0.5; R2 on continuous TM) ===", flush=True)
print(f"  baseline (6 feat, incl region_am)  {fmt(cv_metrics(d2, STRUCT, 'y_tm', 'tm'))}", flush=True)
print(f"  + phyloP                           {fmt(cv_metrics(d2, STRUCT+PH, 'y_tm', 'tm'))}", flush=True)

patho, _ = load_humsavar()
sub = add_overlap_label(d2, patho)
print("\n=== FUNCTIONAL (pathogenic overlap; R2 = label-variance) ===", flush=True)
print(f"  n={len(sub)} pos={int(sub['y'].sum())}", flush=True)
print(f"  baseline (6 feat)                  {fmt(cv_metrics(sub, STRUCT, 'y'))}", flush=True)
print(f"  + phyloP                           {fmt(cv_metrics(sub, STRUCT+PH, 'y'))}", flush=True)
open(f"{SP}/exp6_done.txt", "w").write("done")
