"""Feasibility probe: load HF ESMFold and time one short CPU fold.
Decides whether a real folding run is tractable here."""
import warnings; warnings.filterwarnings("ignore")
import time, torch
from transformers import AutoTokenizer, EsmForProteinFolding
t0 = time.time()
tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
model = model.eval()
model.esm = model.esm.float()  # CPU wants float32
print(f"model loaded in {time.time()-t0:.0f}s", flush=True)
seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"  # ~77 aa
t1 = time.time()
with torch.no_grad():
    inp = tok([seq], return_tensors="pt", add_special_tokens=False)['input_ids']
    out = model(inp)
print(f"folded {len(seq)} aa in {time.time()-t1:.0f}s  plddt_mean={out['plddt'].mean().item():.1f}", flush=True)
open("/private/tmp/claude-501/-Users-arielmelchior-Documents-projects/25923521-8868-4403-b287-93f9380532d6/scratchpad/esmfold_probe_done.txt","w").write("done")
