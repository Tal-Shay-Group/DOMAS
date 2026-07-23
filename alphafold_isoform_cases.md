# AlphaFold-assessed isoform structural changes — benchmark cases for DOMAS

Published cases where splice-isoform protein changes were assessed **visually /
structurally with AlphaFold**, each with an explicit verdict (real structural
change vs. structure preserved). The intent: run DOMAS on the same isoforms and
check whether our impact call agrees.

**Scope for the DOMAS run:** we have no junctions here, so treat the change as
the **whole protein** — compare the canonical (reference) isoform against the
alternative isoform end-to-end (prefix/suffix-trim = the changed region), and
read the DOMAS `impact` / `spade` / HMM-coverage output.

All genes below are human and present in DoChaP (verified), so DOMAS can run.
The one remaining setup task is resolving each paper's *specific* alternative
isoform (e.g. "MDH2 isoform 2", "BAX-δ") to the matching **Ensembl transcript**
in DoChaP.

---

## Cases

Ground truth = the exact UniProt isoform pair the papers assessed with AlphaFold,
with the published pairwise **TM-score** (>0.5 = structurally similar / fold
preserved; <0.5 = structurally different). Both the canonical and the alternate
UniProt isoform sequences were matched **by sequence** against DoChaP's Ensembl
peptides — **all 7 resolved to a DoChaP transcript at 100% identity (EXACT)**, so
every case is reproducible as a real canonical-vs-alt DOMAS run.

Source of TM-scores + identities: Miller et al. *Genome Biol* 2025 (13059-025-03744-x),
Tables S2 (`MOESM3`, "high-identity/low-TM outliers") and S4 (`MOESM5`).

| # | Gene | canonical (acc → DoChaP tx) | alternate isoform (acc → DoChaP tx) | seq id | TM-score | structural verdict | expected DOMAS |
|---|------|-----------------------------|--------------------------------------|--------|----------|--------------------|----------------|
| 1 | **CFHR3** | Q02985 → ENST00000367425 | **Q02985-2** → ENST00000391985 | 81.5% | **0.485** | **REAL CHANGE** (fold differs) | domain dropped/shorter → HIGH |
| 2 | **MDH2** | P40926 → ENST00000315758 | **P40926-2** → ENST00000432020 | 87.6% | **0.480** | **REAL CHANGE** (fold differs) | shorter + functional-site loss → HIGH |
| 3 | **BAX** | Q07812 → ENST00000345358 | **Q07812-4** → ENST00000354470 | 74.5% | **0.478** | **REAL CHANGE** (fold differs) | shorter/dropped → HIGH |
| 4 | **CDK1** | P06493 → ENST00000856537 | **P06493-2** → ENST00000316629 | 80.8% | **0.815** | preserved fold, **functional** (phospho-Thr15 / PTM) | LOW by structure — functional-site flag needed |
| 5 | **CXCR3** | P49682 → ENST00000373693 | **P49682-2** → ENST00000373691 | 87.7% | **0.870** | preserved fold, **functional** (ligand-binding pocket) | LOW by structure — functional-site flag needed |
| 6 | **SRP9** | P49458 → ENST00000304786 | **P49458-2** → ENST00000366839 | 49.5% | **0.758** | **NO CHANGE** (low id, fold preserved) | preserved → none/LOW *(true negative)* |
| 7 | **EEF1AKMT3** | Q96AZ1 → ENST00000300209 | **Q96AZ1-2** → ENST00000333012 | 43.4% | **0.636** | **PARTIAL** (retains 4/7 SAM sites, rearranged) | partial → LOW–MODERATE |

**Corrections vs. the earlier draft** (from pinning the exact isoform):
- **BAX** ground truth is **Q07812-4** with **TM 0.478 → a real structural change**
  (not merely a buried→exposed PTM as first written).
- **CXCR3** ground truth is **P49682-2** with **TM 0.870 → fold preserved**; the
  change is functional (binding pocket), *not* a structural "added region". This
  flips the earlier "REAL CHANGE / longer" expectation.
- **CDK1** ground truth **P06493-2**, **TM 0.815 → fold preserved**; the biology
  is a functional-residue (Thr15) effect, not a fold change.

### Cases without an exact isoform ID yet (separate source)

| # | Gene | ENSG | canonical tx | note |
|---|------|------|--------------|------|
| 8 | **SEPTIN9** | ENSG00000184640 | ENST00000427177 | Qualitative/evolutionary case: AS regions correlate with loops/IDRs across species homologs. No isoform-pair structural comparison. |
| 9 | **MAPT** (Tau) | ENSG00000186868 | ENST00000262410 | Actually a **disease-missense** case (L284R, S285R in the R3/R4 repeats near splice regions) — *not* an isoform-pair (3R/4R) structural comparison. DoChaP canonical (833aa) ≠ UniProt P10636 (758aa). |

**Cases 8–9 are NOT usable as quantitative benchmark cases.** Their source paper
(Yang et al., bioRxiv 2024.01.30.578053) provides **no structural ground truth**:
zero TM-scores, zero UniProt accessions, and no exact isoform pairs anywhere in
the text or its Tables S2–S4 (which cover cross-species homologs and aggregate
secondary-structure counts, not isoform-vs-isoform structure). Both case studies
are intrinsically-disordered-region events, where AlphaFold confidence is low by
construction — which is why that paper stays qualitative. There is nothing to
pin and nothing to compare DOMAS against, so they are dropped (not merely
"pending"). The 7 exact cases above are the entire quantitative benchmark.

---

## Results — whole-protein DOMAS on the exact pairs

Ran each **exact** canonical-tx vs alternate-tx pair (no junctions, whole-protein
scope) through the DOMAS enrichment impact (`hmm_change_impact`: Pfam coverage
loss × pLDDT × UniProt fold/functional-site × AlphaMissense constraint).

| Gene | canon→alt tx | Δlen | changed region | DOMAS **impact** | key signals | TM | structural truth | call |
|------|--------------|:----:|:--------------:|:----------------:|-------------|:----:|------------------|------|
| CFHR3 | 367425→391985 | −61 | 144–204 | **MODERATE** | Sushi PF00084, disulfides, AM 0.43 | 0.485 | **CHANGE** | ✓ caught (under-graded) |
| MDH2 | 315758→432020 | −42 | 144–185 | **HIGH** | folded, BINDING+2×MOD_RES, AM 0.74, Pfam→84–91% | 0.480 | **CHANGE** | ✓ caught |
| BAX | 345358→354470 | −49 | 30–78 | **HIGH** | folded, BH3 MOTIF, Pfam PF00452→88% | 0.478 | **CHANGE** | ✓ caught |
| CDK1 | 856537→316629 | −57 | 107–163 | **HIGH** | folded, ACT_SITE+phospho-Thr, AM 0.87, kinase→40% | 0.815 | preserved (functional) | ✗ over-call by structure / defensible by function |
| CXCR3 | 373693→373691 | **+47** | 2–4 | **NONE** | 7TM PF00001 100→100, no sites in trimmed region | 0.870 | preserved (functional) | ✓ on fold; **misses** the functional N-term insertion |
| SRP9 | 304786→366839 | −4 | 48–86 | **HIGH** | folded, acetyl-K, AM 0.67, Pfam PF05486→71% | 0.758 | **preserved (true negative)** | ✗ **false positive** |
| EEF1AKMT3 | 300209→333012 | −77 | 98–226 | **HIGH** | folded, 3× SAM BINDING, Pfam PF10294→36% | 0.636 | partial | ~ over-call (partial change is real) |

**Sensitivity (real structural change, TM < 0.5): 3/3 flagged** — CFHR3, MDH2, BAX
all called MODERATE/HIGH. DOMAS never missed a genuine fold change.

**Specificity is the weak axis, exactly as the professor predicted:**
- **SRP9 is the clean failure** — TM 0.758 (fold preserved), the paper's own
  true-negative, yet DOMAS says HIGH. The C-terminal AS region (48–86) diverges
  in sequence (≈50% id), damages Pfam coverage (71%), is folded and hits an
  acetyl-lysine — every cheap signal fires — but the divergent sequence **re-folds
  into the same shape**, which only the 3D structure reveals. No sequence/HMM/
  feature signal can see this.
- **CDK1 & EEF1AKMT3** are over-calls against *global fold* (TM > 0.5) but both
  have documented functional changes (CDK1 loses the Thr15 phospho-region / part
  of the kinase domain; EEF1AKMT3 loses 3 SAM-binding residues). So HIGH is wrong
  on the TM axis but not biologically empty — impact measures functional-site /
  domain-integrity disruption, a *different axis* than global fold similarity.
- **CXCR3 exposes the insertion blind spot**: the change is a **+47-aa N-terminal
  insertion**, but `changed_region` is anchored on canonical coordinates, where an
  insertion leaves almost nothing (trimmed to 2–4). Impact = NONE. DOMAS agrees
  with the preserved fold but is blind to the functional (ligand-binding) effect.

**Bottom line:** with the isoform pairing now verified (all 7 exact, 100% id),
the earlier qualitative read holds on solid ground — DOMAS has excellent
sensitivity to real structural change and a real specificity gap on
sequence-divergent-but-refolding isoforms (SRP9) and on insertions (CXCR3).
Impact is best read as a *functional-disruption* score, not a fold-change
predictor; the two agree on the clear cases and diverge exactly where the
professor said only visual AlphaFold inspection is definitive.

---

## Sources

- **Cases 1–7:** *Predicting the structural impact of human alternative
  splicing*, **Genome Biology 2025**, doi:10.1186/s13059-025-03744-x —
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC12442299/> (local copy:
  `papers/Predicting the structural impact of human alternative splicing.pdf`).
  AlphaFold2 structures of SwissProt + CHESS isoforms; cases 1–7 are its worked
  examples (CFHR3, MDH2, BAX, CDK1, CXCR3, SRP9, EEF1AKMT3).
- **Cases 8–9:** *Systematic characterization of protein structural features of
  alternative splicing isoforms using AlphaFold 2*, bioRxiv 2024.01.30.578053v2
  (local copy: `papers/2024.01.30.578053v2.full.pdf`). This paper is mostly a
  statistical enrichment analysis; its named case studies are **Septin-9** (AS
  in disordered loops) and **Tau/MAPT** (missense mutations in the 2N4R R1/R2
  repeats; the 3R/4R split is the classic exon-10 AS event). NOTE: it does *not*
  contain BAX — that example is from the Genome Biology 2025 paper.
- **Tool / method context:** *IsoformMapper* (protein-network community analysis
  of splice variants), bioRxiv 2025.03.05.641708 —
  <https://www.biorxiv.org/content/10.1101/2025.03.05.641708v1.full>.

## Overall finding from 578053 (context, not a single case)
Alternative-splicing regions are **significantly enriched in loops and highly
exposed / disordered residues** — i.e. most AS happens in disordered regions, so
it rarely disrupts the folded core but often exposes/buries functional sites.
This directly supports our impact design (pLDDT/fold downweight for disordered
regions, *rescued* by functional-site + AlphaMissense signals) and predicts that
a good chunk of real AS effects are the "disordered but functional" kind that a
pure fold-change metric would miss.

## Dataset-level context (for a possible bulk comparison later)
The Genome Biology 2025 study flags **328 isoforms with high sequence identity
(>70%) but low TM-score (<0.5)** — i.e. big structural change from a small
sequence change — of which **53 show clear structural domain alterations**.
These 328/53 are a ready-made "hard positive" set to test whether DOMAS's
impact score flags structure changes that sequence identity alone misses.

## Notes / caveats
- "No-change" examples are under-represented in the literature (papers emphasize
  changes); only SRP9 (clear) + EEF1AKMT3 (partial) here. Worth adding more
  preserved-structure cases before treating this as a balanced benchmark.
- Verdicts are the *papers'* structural calls, not experimental ground truth.
- Next step: map each alternative isoform to a DoChaP Ensembl transcript, run
  DOMAS whole-protein (canonical vs alt), and tabulate DOMAS impact vs the
  AlphaFold verdict here.
