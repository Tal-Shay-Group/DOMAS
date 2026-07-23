"""Per-pair contingency of variant class (pathogenic vs benign humsavar overlap)
against DOMAS impact level, with row% and col%. Unit = isoform-pair changed
region (one alternative-splicing event), the natural DOMAS unit."""
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

fp = pd.read_csv(f"{SP}/full_pairs.csv").rename(columns={'isoform': 'iso'})
br = pd.read_csv(f"{SP}/bench_rich_results.csv")[['iso', 'rank']]
d = fp.merge(br, on='iso'); d = d[d['kind'] != 'insertion'].dropna(subset=['canon_lo']).copy()
d['canon_lo'] = d['canon_lo'].astype(int); d['canon_hi'] = d['canon_hi'].astype(int)
ov = lambda a, lo, hi, t: 1 if t.get(a) and any(lo <= p <= hi for p in t[a]) else 0
d['patho'] = [ov(a, lo, hi, patho) for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]
d['benign'] = [ov(a, lo, hi, benign) for a, lo, hi in zip(d['acc'], d['canon_lo'], d['canon_hi'])]

lvls = [('none', 0), ('low', 1), ('moderate', 2), ('high', 3)]
counts = {cls: [len(d[(d[cls] == 1) & (d['rank'] == r)]) for _, r in lvls] for cls in ['patho', 'benign']}
colsum = [counts['patho'][i] + counts['benign'][i] for i in range(len(lvls))]
print(f"unit = isoform-pair changed regions overlapping a variant of that class")
print(f"proteins: Pathogenic={d[d['patho']==1]['acc'].nunique()}  Benign={d[d['benign']==1]['acc'].nunique()}\n")
print("class      | " + " | ".join(l.center(20) for l, _ in lvls) + " | total")
for cls, disp in [('patho', 'Pathogenic'), ('benign', 'Benign')]:
    c = counts[cls]; tot = sum(c)
    cells = [f"{c[i]} (r{100*c[i]/tot:.0f}% c{100*c[i]/colsum[i]:.0f}%)".center(20) for i in range(len(lvls))]
    print(f"{disp:10} | " + " | ".join(cells) + f" | {tot}")
print("col total  | " + " | ".join(str(x).center(20) for x in colsum) + f" | {sum(colsum)}")
