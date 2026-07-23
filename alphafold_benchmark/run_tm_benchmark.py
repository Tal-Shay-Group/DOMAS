import sys, os, gzip, csv, sqlite3, subprocess, time, collections
sys.path.insert(0, "/Users/arielmelchior/Documents/projects/DOMAS/code")
import build
from enrichment import Enricher
from utils import hmm_change_impact, insertion_impact

DATA = "/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-db/data"
SP = "/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad"
PFAM = f"{DATA}/pfam/Pfam-A.hmm"
RANK = {'none': 0, 'gain': 1, 'low': 1, 'moderate': 2, 'high': 3}
t0 = time.time()

def load_fasta(path, idfield=1):
    out = {}
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        cur = None; buf = []
        for ln in f:
            if ln.startswith('>'):
                if cur: out[cur] = ''.join(buf)
                h = ln[1:]
                cur = h.split('|')[idfield] if '|' in h else h.split()[0]
                buf = []
            else:
                buf.append(ln.strip())
        if cur: out[cur] = ''.join(buf)
    return out

print("loading sequences...", flush=True)
canon = load_fasta(f"{DATA}/uniprot/UP000005640_9606.fasta.gz")          # {acc: seq}
iso = load_fasta(f"{DATA}/uniprot/uniprot_sprot_varsplic.fasta.gz")      # {iso_id: seq}

rows = []
for r in csv.DictReader(open(f"{SP}/table_s4_all.csv")):
    if not r['tm']:
        continue
    acc = r['isoform'].split('-')[0]
    cs, as_ = canon.get(acc), iso.get(r['isoform'])
    if cs and as_:
        rows.append({'iso': r['isoform'], 'acc': acc, 'tm': float(r['tm']),
                     'ident': float(r['identity']) if r['identity'] else None,
                     'cls': r['class'], 'cseq': cs, 'aseq': as_})
print(f"runnable pairs (both seqs present): {len(rows)} / from table", flush=True)

# one hmmsearch over all unique sequences (canonical ids = acc, isoform ids = iso id)
uniq = {}
for r in rows:
    uniq[r['acc']] = r['cseq']
    uniq[r['iso']] = r['aseq']
domtbl = f"{SP}/bench_full.domtbl"
if not os.path.exists(domtbl):
    fa = f"{SP}/bench_full.fa"
    with open(fa, 'w') as f:
        for i, s in uniq.items():
            f.write(f">{i}\n{s}\n")
    print(f"hmmsearch over {len(uniq)} sequences...", flush=True)
    subprocess.run(['hmmsearch', '--cut_ga', '--cpu', '8', '--domtblout', domtbl, PFAM, fa],
                   check=True, capture_output=True)
    print(f"hmmsearch done at {time.time()-t0:.0f}s", flush=True)

matches = collections.defaultdict(list)
for row in build._parse_hmmsearch_domtbl(domtbl):   # (prot,acc,name,ali_s,ali_e,...,cov=12)
    matches[row[0]].append((row[1], row[10], row[12], row[3], row[4]))  # acc,bits,cov,ali_s,ali_e

db = sqlite3.connect("/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-web/DB_merged.sqlite")
e = Enricher(f"{DATA}/enrichment.sqlite", PFAM, db)
pl_cache, am_cache = {}, {}
def plddt(acc):
    if acc not in pl_cache: pl_cache[acc] = e.plddt(acc)
    return pl_cache[acc]
def am(acc):
    if acc not in am_cache: am_cache[acc] = e.am_pathogenicity(acc)
    return am_cache[acc]

def impact(r):
    canon_id, iso_id, acc = r['acc'], r['iso'], r['acc']
    cm = {m[0]: m for m in matches.get(canon_id, [])}
    am_m = {m[0]: m for m in matches.get(iso_id, [])}
    reg = Enricher.changed_region(r['cseq'], r['aseq'])
    fold = None; sites = []; ram = None; added = 'none'
    if reg and reg['canon']:
        lo, hi = reg['canon']
        fold = e.fold_state(acc, lo, hi)
        sites = e.functional_sites(acc, lo, hi)
        arr = am(acc)
        if arr:
            seg = [v for v in arr[lo-1:hi] if v is not None]
            ram = round(sum(seg)/len(seg), 3) if seg else None
    if reg and reg['alt']:
        alo, ahi = reg['alt']; alen = ahi-alo+1
        clen = (reg['canon'][1]-reg['canon'][0]+1) if reg['canon'] else 0
        net = alen-clen
        if net > 0:
            indom = any(min(ahi, m[4])-max(alo, m[3])+1 >= 0.5*alen for m in matches.get(iso_id, []))
            if reg['canon']:
                flo, fhi = reg['canon']
            else:
                cl = len(r['cseq']); fl = alo-1; flo, fhi = max(1, fl), min(cl, fl+1)
            jf = e.fold_state(acc, flo, fhi) if fhi >= flo else None
            added = insertion_impact(net, indom, jf)
    hits = bool(sites)
    best = (2, 'moderate') if hits else (0, 'none')
    if RANK[added] > best[0]: best = (RANK[added], added)
    pl = plddt(acc)
    for a in set(cm) | set(am_m):
        ccov = cm[a][2] if a in cm else None
        acov = am_m[a][2] if a in am_m else None
        cpl = None
        if a in cm and pl:
            s, en = cm[a][3], cm[a][4]; sg = pl[s-1:en] if en <= len(pl) else []
            cpl = round(sum(sg)/len(sg), 1) if sg else None
        lab = hmm_change_impact(ccov, acov, cpl, fold_state=fold, hits_functional_site=hits, region_am=ram)
        if RANK[lab] > best[0]: best = (RANK[lab], lab)
    kind = reg['kind'] if reg else 'identical'
    return best[1], kind

out = open(f"{SP}/bench_full_results.csv", 'w', newline='')
w = csv.writer(out); w.writerow(['isoform', 'acc', 'tm', 'identity', 'class', 'kind', 'impact'])
n = 0
for r in rows:
    imp, kind = impact(r)
    w.writerow([r['iso'], r['acc'], r['tm'], r['ident'], r['cls'], kind, imp])
    n += 1
    if n % 1000 == 0:
        print(f"  scored {n}/{len(rows)} at {time.time()-t0:.0f}s", flush=True)
out.close()
print(f"scored {n} pairs at {time.time()-t0:.0f}s", flush=True)

# confusion vs TM verdict (CHANGE if tm<0.5)
res = list(csv.DictReader(open(f"{SP}/bench_full_results.csv")))
def verdict(tm): return 'CHANGE' if float(tm) < 0.5 else 'PRESERVED'
for thr_name, positive in [('impact>=moderate', {'moderate', 'high'}),
                           ('impact>=low', {'low', 'moderate', 'high', 'gain'}),
                           ('impact==high', {'high'})]:
    tp = fp = tn = fn = 0
    for x in res:
        v = verdict(x['tm']); pos = x['impact'] in positive
        if v == 'CHANGE' and pos: tp += 1
        elif v == 'CHANGE' and not pos: fn += 1
        elif v == 'PRESERVED' and pos: fp += 1
        else: tn += 1
    sens = tp/(tp+fn) if tp+fn else 0
    spec = tn/(tn+fp) if tn+fp else 0
    print(f"\n[{thr_name}]  TP={tp} FN={fn} FP={fp} TN={tn}"
          f"  sensitivity={sens:.2f} specificity={spec:.2f}", flush=True)
# impact distribution by TM band
print("\nimpact distribution by TM band:", flush=True)
bands = [('<0.5 CHANGE', lambda t: t < 0.5), ('0.5-0.85', lambda t: 0.5 <= t < 0.85), ('>=0.85 PRESERVED', lambda t: t >= 0.85)]
for bn, bf in bands:
    c = collections.Counter(x['impact'] for x in res if bf(float(x['tm'])))
    print(f"  {bn:18} n={sum(c.values()):5}  {dict(c)}", flush=True)
open(f"{SP}/bench_full_done.txt", 'w').write(f"{time.time()-t0:.0f}s\n")
print(f"DONE {time.time()-t0:.0f}s", flush=True)
