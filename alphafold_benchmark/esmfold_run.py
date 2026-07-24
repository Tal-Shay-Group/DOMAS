"""Task 9: ESMFold the isoform AND canonical (method-matched), compute TM between
them with tmtools, and test whether ESMFold-TM reproduces the paper's AF2-TM label
-- especially on the hard high-identity/low-TM outliers (the SRP9 class) that cheap
features confidently mis-call. If folding recovers the label, it resolves the
uncertain band that the cheap-feature ceiling cannot."""
import warnings; warnings.filterwarnings("ignore")
import time, io, numpy as np, pandas as pd, torch
from transformers import AutoTokenizer, EsmForProteinFolding
from tmtools import tm_align
from exp_common import load_fasta, DATA, SP

canon = load_fasta(f"{DATA}/uniprot/UP000005640_9606.fasta.gz")
iso = load_fasta(f"{DATA}/uniprot/uniprot_sprot_varsplic.fasta.gz")
sel = pd.read_csv(f"{SP}/esmfold_pairs.csv")

t0 = time.time()
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True).eval()
model.esm = model.esm.float()
model.trunk.set_chunk_size(64)
print(f"model loaded {time.time()-t0:.0f}s", flush=True)

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
def fold_ca(seq):
    with torch.no_grad():
        ids = tok([seq], return_tensors="pt", add_special_tokens=False)['input_ids']
        out = model(ids)
    pdb = model.output_to_pdb(out)[0]
    ca, s = [], []
    for ln in pdb.splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
            ca.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            s.append(AA3.get(ln[17:20].strip(), 'X'))
    return np.array(ca, np.float64), ''.join(s)

cache = {}
def get(seq):
    if seq not in cache: cache[seq] = fold_ca(seq)
    return cache[seq]

rows = []; tf = time.time()
for i, r in sel.iterrows():
    cs = canon.get(r['acc']); a = iso.get(r['iso'])
    if not cs or not a: continue
    cca, cseq = get(cs); aca, aseq = get(a)
    res = tm_align(cca, aca, cseq, aseq)
    etm = max(res.tm_norm_chain1, res.tm_norm_chain2)
    rows.append({'iso': r['iso'], 'paper_tm': r['tm'], 'esmfold_tm': etm,
                 'ident': r['ident_rt'], 'outlier': int(r['ident_rt'] >= 80 and r['tm'] < 0.5)})
    print(f"  [{i+1}/{len(sel)}] {r['iso']} paperTM={r['tm']:.2f} esmTM={etm:.2f} ({time.time()-tf:.0f}s)", flush=True)
res = pd.DataFrame(rows); res.to_csv(f"{SP}/esmfold_results.csv", index=False)

print(f"\nfolded {len(cache)} unique seqs, {len(res)} pairs in {time.time()-tf:.0f}s", flush=True)
print(f"Pearson(esmfold_tm, paper_tm) = {res[['esmfold_tm','paper_tm']].corr().iloc[0,1]:.3f}", flush=True)
print(f"Spearman = {res[['esmfold_tm','paper_tm']].corr('spearman').iloc[0,1]:.3f}", flush=True)
# agreement on the CHANGE call (TM<0.5)
res['paper_change'] = res['paper_tm'] < 0.5; res['esm_change'] = res['esmfold_tm'] < 0.5
acc = (res['paper_change'] == res['esm_change']).mean()
print(f"CHANGE-call agreement (TM<0.5): {acc:.2f}  (n={len(res)})", flush=True)
o = res[res['outlier'] == 1]
print(f"\nhigh-id/low-TM outliers (n={len(o)}): ESMFold correctly calls CHANGE on "
      f"{int((o['esmfold_tm']<0.5).sum())}/{len(o)}", flush=True)
print(o[['iso','ident','paper_tm','esmfold_tm']].to_string(index=False), flush=True)
open(f"{SP}/esmfold_run_done.txt", "w").write("done")
