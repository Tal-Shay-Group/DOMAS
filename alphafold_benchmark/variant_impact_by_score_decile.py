"""Pathogenic/benign variant overlap by DOMAS-score DECILE (10 equal-population
bins), the finer-grained companion to variant_impact_contingency.py. DOMAS emits
a categorical impact; this uses the continuous reconstruction (coverage-loss +
functional-site + AlphaMissense + insertion, same as the domas_vs_tm scatter)."""
import re, collections
import pandas as pd
SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"

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

d = pd.read_csv(f"{SP}/bench_rich_results.csv")
d['domas_score'] = (d['max_cov_loss'].clip(lower=0) + 25 * d['func_site']
                    + 50 * (d['region_am'].fillna(0.5) - 0.5).clip(lower=0) + 0.3 * d['net_added'])
fp = pd.read_csv(f"{SP}/full_pairs.csv").rename(columns={'isoform': 'iso'})
d = d.merge(fp[['iso', 'canon_lo', 'canon_hi']], on='iso')
d = d[d['kind'] != 'insertion'].dropna(subset=['canon_lo']).copy()
d['canon_lo'] = d['canon_lo'].astype(int); d['canon_hi'] = d['canon_hi'].astype(int)
ov = lambda a, lo, hi, t: 1 if t.get(a) and any(lo <= p <= hi for p in t[a]) else 0
d['pat'] = [ov(a, lo, hi, patho) for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]
d['ben'] = [ov(a, lo, hi, benign) for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]

rows = pd.concat([d[d['pat'] == 1][['domas_score']].assign(cls='Pathogenic'),
                  d[d['ben'] == 1][['domas_score']].assign(cls='Benign')]).reset_index(drop=True)
rows['bin'], edges = pd.qcut(rows['domas_score'], 10, retbins=True, duplicates='drop', labels=False)
tab = rows.groupby(['bin', 'cls']).size().unstack(fill_value=0)
totP, totB = tab['Pathogenic'].sum(), tab['Benign'].sum()
print(f"DOMAS-score deciles (equal count). patho={totP} benign={totB} base%patho={100*totP/(totP+totB):.0f}\n")
print(f"{'decile (score range)':24} {'Pathogenic(r%,c%)':>20} {'Benign(r%,c%)':>20}  binN  %patho")
for b in tab.index:
    p = int(tab.loc[b].get('Pathogenic', 0)); bn = int(tab.loc[b].get('Benign', 0)); cs = p + bn
    print(f"  {int(b)+1:2} [{edges[int(b)]:6.1f}-{edges[int(b)+1]:6.1f}) "
          f"{f'{p} ({100*p/totP:.0f}%,{100*p/cs:.0f}%)':>20} {f'{bn} ({100*bn/totB:.0f}%,{100*bn/cs:.0f}%)':>20} {cs:5}  {100*p/cs:.0f}%")
