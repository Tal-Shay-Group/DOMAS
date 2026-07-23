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

## Scaled benchmark — 11,068 isoform pairs vs TM-score

The 7 hand-picked cases are too few to measure specificity, so the whole
Genome Biology 2025 dataset was turned into a benchmark. Its supplementary
**Table S4** lists **11,161 isoforms (5,923 genes)** with an exact UniProt
isoform ID, TM-score, sequence identity, pLDDT and structural features;
**Table S2** is 328 curated high-identity / low-TM "real change" outliers. Data
in `alphafold_benchmark/` (`table_s4_all.csv`, `table_s2_changes.csv`), runner
`alphafold_benchmark/run_tm_benchmark.py`, per-pair output `bench_full_results.csv`.

No DoChaP transcript is needed to score DOMAS impact vs TM: impact needs only the
two protein sequences + the canonical UniProt accession (the enrichment key).
Isoform sequences came from UniProt `varsplic` (99% coverage), canonicals from
the human reference proteome (100%). One hmmsearch over 16,948 unique sequences,
then the updated (insertion-aware) impact scorer on **11,068 runnable pairs**.

**Verdict = TM<0.5 (CHANGE) vs TM≥0.5 (PRESERVED). Base rate of CHANGE = 32.8%.**

| operating point | sensitivity | specificity | precision | lift over base |
|-----------------|:-----------:|:-----------:|:---------:|:--------------:|
| impact ≥ moderate | 0.79 | 0.29 | 35% | 1.08× |
| impact = high     | 0.57 | 0.58 | 40% | 1.21× |

Sharper cuts:
- **Curated real-changes (Table S2, 327):** only **68%** flagged ≥moderate (28%
  scored "none"). The pilot's 3/3 sensitivity was small-sample luck.
- **Clean negatives (TM≥0.95, 392):** only **51%** correctly "none/low" — nearly
  half of near-identical-fold isoforms get flagged moderate/high.
- **Monotone but shallow:** mean impact-rank falls 2.26 (TM<0.4) → 1.22 (TM≥0.95).
  Signal exists; it is weak.

**The dominant confound is sequence identity.** Mean identity by impact label:
none 86% · moderate 81% · **high 62.6%**. Impact largely tracks how much sequence
diverged — and TM tracks identity too (the paper's central finding), so most of
the impact↔TM association is both variables following identity. Controlling for
identity, DOMAS adds little independent *structural* signal.

**Honest conclusion.** At scale, DOMAS impact is a **weak predictor of global
fold change** (high-impact enriches for real change only ~1.2× over base rate).
It is best read as a *functional-disruption / domain-integrity* flag whose main
correlate is sequence divergence — not a TM-score surrogate. Two caveats keep
this from being purely a negative result: (1) TM-score is global-fold
superposition and under-weights peripheral domain loss (a dropped peripheral
domain can still give TM>0.85), so some "false positives" are real domain losses
TM misses; (2) neither TM nor impact measures biological significance — the
structure-biology point that only visual AlphaFold inspection is definitive. The
SRP9 specificity gap from the pilot generalises across 11k pairs.

Independent cross-check not yet integrated: **ASpdb** (>7,200 alt-isoform AF2
structures with comparative structural-alteration calls) — a second, non-TM label
source to test whether the specificity gap is TM-metric-specific.

### AlphaFold burial (SASA/RSA) + DSSP — the first identity-orthogonal signal

Since impact (a UniProt/Pfam-derived score) collapses to zero once identity is
controlled, we tested whether the **canonical AlphaFold 3D structure** carries
signal identity can't — at **full scale**. For every runnable pair we downloaded
the canonical AFDB model (5,880 structures) and computed, over the changed region:
mean **RSA** (exposure, Biopython Shrake–Rupley), **buried fraction** (RSA<0.25),
and DSSP-style secondary structure (pydssp). Features land for **10,232 / 11,068
pairs** (the ~7.5% pure insertions have no canonical region). Data + extractor in
`alphafold_benchmark/` (`full_features.csv`, `extract_burial_dssp_full.py`).

**Burial predicts TM and survives identity control** (unlike impact) — full 11k:

| feature | Spearman vs TM | partial \| identity | partial \| identity + region-length |
|---------|:--------------:|:-------------------:|:-----------------------------------:|
| buried_frac | +0.13 | **+0.30** | **+0.30** |
| mean_rsa    | −0.14 | **−0.28** | −0.28 |
| structured_frac (H+E) | +0.08 | +0.17 | +0.17 |
| *impact rank (for contrast)* | −0.19 | *+0.05* | — |

Non-zero in **every** identity band (Spearman(buried,TM) +0.22…+0.37). Incremental
variance: R²(TM ~ identity) = **0.24** → R²(TM ~ identity + buried_frac) = **0.31**
(+0.069, ~+29%); region length adds nothing on top (0.310). So burial is a genuine
identity-independent structural signal — the one thing we found that improves
fold-change prediction, robust from the 2k sample to the full 11k.

**The direction is counter-intuitive and matters.** Real changes (TM<0.5) have
*more exposed* changed regions (mean RSA 0.50, buried 0.21); preserved isoforms
have *more buried* ones (0.44, buried 0.30). So **exposed/peripheral changes drive
low global TM**, not buried-core losses. Mechanistically: TM is dominated by the
compact core, so a small buried substitution is scaffolded and superposes (high
TM), whereas large alternative exons sit in exposed loops/termini and form a big
non-superposing fraction (low TM). This is the **opposite sign** to the current
impact scorer's "structured/folded → raise level" weighting — but note that
weighting targets *functional* importance, not TM, so it should not be blindly
flipped; the two targets differ.

**Net:** AlphaFold is not exhausted by pLDDT. Burial of the changed region adds
real, identity-orthogonal predictive power for fold change (impact did not). The
remaining ceiling (does the *remainder* refold — the SRP9 class) still needs the
isoform actually folded, which is what the TM ground truth already encodes.

**Predictive performance (5-fold CV, `predict_eval.py`).** Predicting the change
verdict (TM<0.5; base rate 33.8%, so majority-class accuracy 66.2%) on the 10,232
non-insertion pairs:

| model | AUC | accuracy | R²(TM) |
|-------|:---:|:--------:|:------:|
| DOMAS impact rank (current) | 0.568 | 0.662 | — |
| identity + impact *(control)* | 0.696 | 0.722 | — |
| sequence identity | 0.697 | 0.722 | 0.239 |
| burial alone (buried_frac) | 0.586 | 0.662 | 0.027 |
| **identity + burial** | **0.755** | **0.738** | **0.307** |
| identity + burial + SS + RSA | 0.757 | 0.738 | 0.306 |

Reads: (1) DOMAS impact alone is a near-useless structural classifier — AUC 0.568,
accuracy = the majority baseline; "identity + impact" (0.696) equals identity
alone (0.697), so impact adds nothing. (2) Burial alone is weak (AUC 0.586, R²
0.03) because identity dominates absolute variance, but as an orthogonal add-on it
lifts identity 0.697 → **0.755** AUC / 0.239 → **0.307** R². (3) The cheap ceiling
is AUC ≈ 0.755 / ~74%: whether the remainder *refolds* needs the isoform folded,
which is the TM label itself.

**Contact modularity is redundant with burial** (`extract_contacts.py`,
`contact_features.csv`). Contacts from the changed region to the rest of the fold
(Cβ<8Å) correlate with TM on their own (Spearman +0.46–0.48) and survive identity
control (partial +0.18–0.25), but once burial is also controlled the partial drops
to ~+0.10, and incremental R² over identity+burial is only +0.006 (0.308 → 0.314).
Contacts and burial measure the same packing (mutual Spearman +0.42), so **burial
is the compact summary of the canonical-structure geometry** — contacts (and, by
extension, likely PAE, the one untested AF signal) add little on top. Burial is the
signal worth carrying forward.

### Numeric score vs TM, identity control, and where the errors fall

(`alphafold_benchmark/analyze_identity_control.py`, `bench_rich_results.csv`.)

**(B) Do the DOMAS numeric scores correlate with TM?** Weakly, in the right
direction (higher impact → lower TM):

| DOMAS numeric | Spearman vs TM | Pearson vs TM |
|---------------|:--------------:|:-------------:|
| impact rank (0–3)         | −0.19 | −0.17 |
| max Pfam coverage loss    | −0.27 | −0.31 |
| region AlphaMissense      | −0.09 | −0.10 |
| constructed domas_score   | −0.28 | −0.31 |
| **sequence identity**     | **+0.46** | **+0.49** |

Sequence identity alone tracks TM ~2× more strongly than any DOMAS score.

**(A) Control for identity — the decisive test.** Within fixed identity bands the
score→TM relationship **vanishes**: Spearman(rank, TM) is ≈0 in every band
(−0.08…+0.09), and mean TM for "none" vs "high" calls is nearly equal within a
band (e.g. identity 85–90%: none 0.696 vs high 0.697). Partial correlation
collapses to zero:

- partial Pearson(rank, TM | identity): **−0.17 → +0.05**
- partial Pearson(domas_score, TM | identity): **−0.31 → −0.03**

So once sequence divergence is held constant, DOMAS's fold/domain/functional-site
signals carry **essentially no additional information about global fold change.**
The whole (weak) impact↔TM association was identity acting through both. For the
task "predict TM," DOMAS impact ≈ a sequence-identity threshold, no better.

**(C) Are the errors near the TM=0.5 boundary?** No — many are gross, not
borderline. At impact≥moderate: FP=5250, FN=770. Median |TM−0.5| for false
positives is **0.238** (a typical FP has TM≈0.74, well inside "preserved"), and
**1,279 FPs have TM≥0.85** (clearly-preserved folds flagged moderate/high). 198
FNs have TM<0.3 (clear changes missed entirely). Treating the levels as ordered
against TM-implied levels, only 27% match exactly, 36% are adjacent, and **37%
are off by ≥2 levels** — half of all disagreements are *not* adjacent.

**Bottom line.** DOMAS impact is a noisy proxy for sequence identity with no
demonstrable independent power to predict global fold change, and its errors are
frequently severe rather than borderline. The honest caveat still stands: TM is
global-fold superposition and is the wrong ground truth for *functional* /
domain-integrity disruption (a dropped peripheral domain can leave TM high). This
benchmark therefore refutes "impact predicts structure (TM)"; it cannot confirm
or refute "impact predicts functional significance," which needs functional
ground truth — the visual/experimental standard the structure biologist named.

---

## Functional axis — pathogenic-variant overlap (the first positive result)

TM measures *fold*, not *biological significance*; on TM, impact added nothing
beyond identity. So we tested the axis DOMAS is actually for: does the changed
region harbour disease variants? Label = the changed region overlaps a
**pathogenic/likely-pathogenic (LP/P)** variant from UniProt **humsavar** (33,233
LP/P variants, annotated in canonical coordinates — no genomic mapping needed).
Restricted to proteins that have ≥1 LP/P variant (so a positive is possible).
`alphafold_benchmark/clinvar_enrichment.py`.

Coverage: of 10,233 non-insertion pairs, **2,023 pairs (998 disease proteins)**
are testable; **907 (44.8%)** overlap a pathogenic variant.

**Impact tracks pathogenic overlap strongly** — raw rate rises none 8.5% → low
17.8% → moderate 34.1% → **high 63.5%** (7.5×). Region length is a big confound
(longer regions catch more variants; AUC 0.813 alone), but the enrichment
**survives length control** — within every length band:

| region length | none | moderate | high |
|---------------|:----:|:--------:|:----:|
| 1–50 aa   | 0.05 | 0.17 | 0.31 |
| 50–120 aa | 0.15 | 0.28 | 0.54 |
| 120–300 aa| 0.33 | 0.69 | 0.73 |
| 300+ aa   | —    | 0.72 | 0.91 |

Multivariate 5-fold CV AUC for predicting overlap (disease proteins):

| model | AUC |
|-------|:---:|
| region length alone | 0.813 |
| length + identity (baseline) | 0.817 |
| baseline + DOMAS impact | 0.833 |
| baseline + burial | 0.839 |
| baseline + impact + burial | 0.841 |

**The finding:** unlike on TM (where impact added ~0 after identity control),
impact adds real signal here (+0.016 AUC over baseline; strong within-band
separation), and burial adds even more (+0.022). This is consistent with impact
being a **functional-disruption** score, not a fold predictor — validated for the
first time on an independent, non-structural label.

**Benign contrast — the specificity check.** Repeating the test with *benign*
(LB/B) variants (5,924 pairs, 3,103 proteins) rules out "impact just finds
variant-dense regions." At fixed length (1–50 aa band), impact separates
**pathogenic** overlap 6× (none 0.05 → high 0.31) but is **flat for benign**
(none 0.17 → high 0.17). Multivariate, **both features are pathogenicity-specific**:
impact adds **+0.016 AUC for pathogenic** overlap but **−0.001 for benign**;
burial adds **+0.022 for pathogenic** but **+0.000 for benign**. Benign variants
are spread everywhere (none-impact regions overlap them 28% vs 8.5% for
pathogenic), and neither feature carries information about their location. So the
enrichment is **pathogenicity-specific**, not variant density — DOMAS is finding
functionally important regions, not just polymorphic ones.

**Caveats.** (1) Length dominates and must be controlled (done). (2) *Semi-circular*:
impact uses UniProt functional residues (`func_site`) as a feature, and pathogenic
variants cluster at functional sites — so impact's contribution is partly built
in; **burial is structure-only and non-circular**, making its +0.022 the cleaner
result. (3) Positive-only: "no known variant" is not proof of neutrality, so this
measures recall/enrichment on functional positives, not specificity.

**Per-pair contingency** (`variant_impact_contingency.py`; unit = one splicing
event; row% = share within the class, col% = pathogenic/benign split within that
impact level; proteins: 552 pathogenic, 1,674 benign):

| variant class \ DOMAS impact | none | low | moderate | high | total |
|------------------------------|:----:|:---:|:--------:|:----:|:-----:|
| **Pathogenic** | 28 (r3% c7%) | 13 (r1% c14%) | 188 (r21% c25%) | **678 (r75% c28%)** | 907 |
| **Benign** | 355 (r13% c93%) | 82 (r3% c86%) | 561 (r20% c75%) | 1748 (r64% c72%) | 2746 |
| column total | 383 | 95 | 749 | 2426 | 3653 |

75% of pathogenic-overlapping events are high-impact (vs 3% none); the pathogenic
share within a column climbs monotonically none 7% → high 28%, while benign stays
diffuse. (Rows are different-sized universes, so the column *trend* is the signal,
not the absolute split; a region can overlap both classes.)

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
