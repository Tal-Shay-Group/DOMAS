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

**Per-pair contingency** (`variant_impact_contingency.py`; impact scored *with
burial*, Option B; unit = one splicing event; row% = share within the class,
col% = pathogenic/benign split within that impact level; proteins: 552 pathogenic,
1,674 benign):

| variant class \ DOMAS impact | none | low | moderate | high | total |
|------------------------------|:----:|:---:|:--------:|:----:|:-----:|
| **Pathogenic** | 29 (r3% c7%) | 13 (r1% c14%) | 184 (r20% c23%) | **681 (r75% c29%)** | 907 |
| **Benign** | 369 (r13% c93%) | 77 (r3% c86%) | 599 (r22% c77%) | 1701 (r62% c71%) | 2746 |
| column total | 398 | 90 | 783 | 2382 | 3653 |

75% of pathogenic-overlapping events are high-impact (vs 3% none); the pathogenic
share within a column climbs monotonically none 7% → high 29%, while benign stays
diffuse. (Rows are different-sized universes, so the column *trend* is the signal,
not the absolute split; a region can overlap both classes.)

### An enrichment signal, NOT a classifier

The enrichment above is real but must not be over-read. Everything so far tests
*variant-overlap presence* (does the region overlap a pathogenic variant),
controlling for region length. A harder, clinically-useful question is: **given a
region overlaps some variant, is it the pathogenic or the benign one?** There the
score is near chance, because **both classes pile into moderate/high**:
`P(mod-or-high | pathogenic) = 95%` but `P(mod-or-high | benign) = 84%` — it says
"important" to almost everything, so specificity is only ~16%.

Under an assumed **50:50 pathogenic:benign prior** (using the class-conditional
row rates, so base-rate-independent):

| rule → call Pathogenic | precision (right) | wrong | sensitivity | specificity | balanced acc |
|------------------------|:-----------------:|:-----:|:-----------:|:-----------:|:------------:|
| impact ≥ moderate | 53.2% | 46.8% | 95.4% | 16.2% | 55.8% |
| impact = high     | 54.8% | 45.2% | 75.1% | 38.1% | 56.6% |

So as a pathogenic-vs-benign *classifier* DOMAS impact is barely above a coin flip
(~53% precision, ~47% wrong). The length-controlled enrichment (E16–E17) was
largely carried by region size; conditioning on "a variant is present" and forcing
class balance removes that crutch. **Correct framing: impact is a coarse
enrichment / prioritisation signal ("these high-impact events are enriched for
functional consequence — look here first"), never a per-event classification.**

**Finer view — by continuous DOMAS-score decile** (`variant_impact_by_score_decile.py`;
score = the reconstruction coverage-loss + functional-site + AM + insertion, binned
into 10 equal-count deciles; c% = pathogenic share within the bin; overall base 25%):

| score decile (range) | Pathogenic (r%, c%) | Benign (r%, c%) | bin N | % pathogenic |
|----------------------|:-------------------:|:---------------:|:-----:|:------------:|
| 1 · [0.0–1.2)   | 20 (2%, 5%)   | 346 (13%, 95%) | 366 | **5%** |
| 2 · [1.2–25)    | 122 (13%, 20%)| 475 (17%, 80%) | 597 | 20% |
| 3 · [25–29)     | 40 (4%, 30%)  | 93 (3%, 70%)   | 133 | 30% |
| 4 · [29–45)     | 117 (13%, 32%)| 253 (9%, 68%)  | 370 | 32% |
| 5 · [45–64)     | 106 (12%, 29%)| 256 (9%, 71%)  | 362 | 29% |
| 6 · [64–86)     | 117 (13%, 32%)| 247 (9%, 68%)  | 364 | 32% |
| 7 · [86–103)    | 78 (9%, 21%)  | 288 (10%, 79%) | 366 | 21% |
| 8 · [103–122)   | 95 (10%, 24%) | 296 (11%, 76%) | 391 | 24% |
| 9 · [122–126)   | 73 (8%, 21%)  | 267 (10%, 79%) | 340 | 21% |
| 10 · [126–190)  | 139 (15%, 38%)| 225 (8%, 62%)  | 364 | **38%** |

The score's real resolution is at the **bottom** — the lowest decile (score ≈ 0) is
95% benign, so a near-zero score reliably means benign. Above that it **plateaus and
is non-monotonic** (~20–38%, dips at deciles 7 and 9), and **even the top decile is
only 38% pathogenic** — a maximal DOMAS score is still more likely benign. Same
conclusion as the categorical view, sharpened: good at flagging the benign low end,
weak resolution across the moderate-to-high range.

### Calibrated continuous score (`impact_prob`)

The decile analysis showed the categorical impact throws away signal. So a
**calibrated continuous companion** was built: a logistic model (`utils.impact_probability`,
new `impact_prob` output column) over `region_am` + gnomAD **LOEUF** + `max_cov_loss`
+ `buried_frac`, trained on the humsavar patho-vs-benign label. Provisioning:
gnomAD gene constraint → `gene_constraint` table (`build_gnomad`, `Enricher.loeuf`).
Fit/eval in `alphafold_benchmark/fit_calibrated.py`, coefficients in `calib_model.json`.

What was tried and what survived (`improve_proto.py`):
- **Peak / mean-of-top-K AM: rejected.** For domain-scale regions (median 169 aa)
  the *mean* AM beats every peak variant monotonically (mean 0.754 vs top-1..top-50
  0.681–0.703) — pathogenicity here is about *overall* regional constraint, not a
  hotspot.
- **Not quantising: the big lever.** Categorical impact AUC 0.585 → continuous mean
  AM 0.754 for patho-vs-benign — same signal, un-bucketed.
- **gnomAD LOEUF: the useful new feature** (non-circular): 0.754 → 0.770.

The shipped model: **5-fold CV AUC 0.769** (pair-level), and **0.765 under
gene-grouped CV** — each protein confined to one fold, so no gene appears in both
train and test. The 0.003 gap confirms negligible leakage (1,979 proteins for
3,282 pairs, ~1.7 isoforms each; even LOEUF alone survives gene-grouped CV at
0.654). The model is well-calibrated (predicted prob tracks observed rate), and
the deciles are now **monotonic** with the top decile at **58% pathogenic** (vs
the categorical score's 38%). Full per-pair decile table
(r% = share within class, c% = pathogenic share within bin; patho 907, benign 2,746):

| impact_prob decile (range) | Pathogenic (r%, c%) | Benign (r%, c%) | bin N | % pathogenic |
|----------------------------|:-------------------:|:---------------:|:-----:|:------------:|
| 1 · [0.02–0.07) | 14 (2%, 4%) | 352 (13%, 96%) | 366 | 4% |
| 2 · [0.07–0.11) | 29 (3%, 8%) | 336 (12%, 92%) | 365 | 8% |
| 3 · [0.11–0.15) | 31 (3%, 8%) | 334 (12%, 92%) | 365 | 8% |
| 4 · [0.15–0.19) | 60 (7%, 16%) | 305 (11%, 84%) | 365 | 16% |
| 5 · [0.19–0.24) | 68 (7%, 19%) | 299 (11%, 81%) | 367 | 19% |
| 6 · [0.24–0.29) | 95 (10%, 26%) | 269 (10%, 74%) | 364 | 26% |
| 7 · [0.29–0.35) | 110 (12%, 30%) | 262 (10%, 70%) | 372 | 30% |
| 8 · [0.35–0.44) | 135 (15%, 38%) | 223 (8%, 62%) | 358 | 38% |
| 9 · [0.44–0.58) | 151 (17%, 41%) | 214 (8%, 59%) | 365 | 41% |
| 10 · [0.58–0.90) | 214 (24%, 58%) | 152 (6%, 42%) | 366 | **58%** |

Top-decile % pathogenic across the three scores: categorical high bucket 29% →
constructed `domas_score` 38% → calibrated `impact_prob` **58%**.

`impact_prob` is a companion to (not a replacement for) the categorical `impact`;
it stays honestly bounded (semi-circular via AM, positive-only labels), but it is
calibrated and monotonic where the categorical score was neither.

### Second calibrated score: `fold_change_prob` (structural axis)

The same calibration approach, applied to the **TM (fold-change)** target, gives a
second companion score `fold_change_prob` = P(TM < 0.5). `utils.fold_change_probability`
over identity + burial (buried_frac, mean_rsa) + region_am + LOEUF + max_cov_loss,
trained on the 10,232 TM-labelled changed regions. Fit/eval in
`alphafold_benchmark/fit_foldchange.py`, coefficients in `foldchange_model.json`.
Runtime `identity` is computed from the trimmed changed region (correlates 0.81
with the paper's alignment identity).

**Gene-grouped 5-fold CV AUC = 0.777**, well-calibrated. Notably AM helps here only
*in combination* — useless alone (0.54) or with identity (0.70), but identity +
burial + AM = 0.78 (a suppressor effect). Confusion matrix (out-of-fold, threshold
P ≥ 0.5; base rate TM<0.5 = 34%, majority-class accuracy 66%):

| | actual CHANGE | actual PRESERVED |
|--------------------|:---:|:---:|
| **predicted CHANGE** | TP 1,580 | FP 745 |
| **predicted PRESERVED** | FN 1,874 | TN 6,033 |

**Accuracy 74%**, sensitivity 0.46, specificity 0.89, precision 0.68, balanced 0.67.
(Threshold 0.5 is conservative given the 34% base rate — high specificity, half the
changes missed; the AUC 0.777 is the threshold-independent summary.)

**The two models are different questions — and burial flips sign** (standardised
logistic coefficients):

| feature | `impact_prob` (functional / pathogenic) | `fold_change_prob` (structural / TM) |
|--------------|:---:|:---:|
| identity     | —      | **−0.67** |
| region_am    | **+0.80** | +0.46 |
| buried_frac  | **+0.24** | **−0.93** |
| mean_rsa     | —      | +0.18 |
| loeuf        | **−0.42** | +0.11 |
| max_cov_loss | −0.07  | **+0.54** |
| intercept    | −1.19  | −0.85 |

The **buried_frac sign flip** is the biological duality of the whole
investigation: a buried changed region *harbours pathogenic variants* (functional,
+0.24) yet *preserves the fold* (structural, −0.93). DOMAS now emits both honest,
calibrated probabilities — `impact_prob` for functional relevance and
`fold_change_prob` for structural change — alongside the categorical `impact`.

**Why logistic regression** (`ml_compare.py`). Both scores use L2-regularised
logistic regression, chosen empirically over random forest, gradient boosting,
hist-gradient-boosting and a small neural net (all gene-grouped CV):

| model | Functional AUC | Structural AUC |
|-------|:---:|:---:|
| **logistic (chosen)** | **0.761** | 0.787 |
| logistic + interactions | 0.760 | 0.788 |
| random forest | 0.741 | 0.798 |
| gradient boosting | 0.748 | 0.796 |
| neural net (MLP 32,16) | 0.732 | 0.794 |

Logistic wins outright on the functional target and loses only ~0.01 on the
structural one, while keeping signed, interpretable coefficients (the sign-flip
story) and native calibration. With just 4–6 monotone tabular features, tree
ensembles and NNs have nothing extra to exploit — the ceiling is set by the
features/labels, not the algorithm. (A Bayesian logistic would give the same
point AUC and add only per-prediction uncertainty, not accuracy.)

### Uncertainty routing (`fold_change_call`)

Since the cheap-feature ceiling (~AUC 0.79) is set by the "does the remainder
refold?" residual that only actual folding resolves, the honest output is a
**triage**, not a forced call. `utils.fold_change_call` maps `fold_change_prob` →
`change` (≥0.6), `preserved` (≤0.4), or `uncertain` (route to folding/inspection);
wired into add_scores as a `fold_change_call` column. `selective_prediction.py`.

**Selective prediction** (gene-grouped OOF, structural model): abstaining on the
least-confident cases lifts accuracy on the *called* set from **75% (call all) →
80% at 80% coverage → 86% at 50% coverage**.

The **uncertain band (prob 0.4–0.6) is 17% of pairs** and has an observed change
rate of **0.47** — a genuine coin-flip, i.e. the model correctly expresses "don't
know," not error. Its feature signature is exactly the SRP9 class: more **divergent**
(identity 64 vs 80), more **exposed** (buried_frac 0.20 vs 0.34), more **domain
loss** (cov-loss 43 vs 26) — a big, exposed change whose refolding cheap features
can't judge. Those are the cases to fold (ColabFold) or inspect; the confident 41%
(prob <0.2 or >0.8) can be reported directly.

**Accuracy/AUC by abstention band** (route-band = probability range sent to
folding; the rest is "called"). Structural (`fold_change_prob`, base 34%):

| route-band | % routed | acc in band | acc called | AUC called |
|------------|:---:|:---:|:---:|:---:|
| none (call all) | 0% | — | 75.1% | 0.787 |
| 40–60 | 17% | 54% | 79.4% | 0.801 |
| 35–65 | 26% | 58% | 81.1% | 0.805 |
| 30–70 | 36% | 61% | 83.1% | 0.814 |
| 20–80 | 59% | 66% | 88.3% | 0.822 |

Functional (`impact_prob`, base 28%): 40–60 routes 16% → acc called 78.8%; 20–80
routes 55% → 89.0% — but AUC-on-called *falls* (0.769→0.709), since removing the
middle leaves confident-but-internally-less-rankable extremes. The `acc in band` ≈
coin flip everywhere confirms the routed cases are the genuinely hard ones.

**Widening the trigger for the confidently-wrong outliers** (`selective_prediction`
follow-up). The high-identity / low-TM outliers (identity ≥80% AND TM<0.5; 1,035
pairs = 10%) are *confidently* mis-called `preserved`, so the confidence gate
catches only **17%** of them. Domain-loss widening barely helps (~30%), but
**"high identity AND exposed changed region" (mean_rsa ≥0.5) catches ~73%** — the
outliers are exposed changes that flip the fold (the burial signal). Cost: fold
budget doubles (17%→39%), precision ~19%. A real fix for the blind spot, but one
for **high-stakes/clinical** use, not casual screening.

### Breaking the ceiling — protein language model embeddings

Every hand-crafted feature, interaction, model family, cascade and specialist stalled
at ~AUC 0.79 (the information ceiling). The one orthogonal data type — **pLM
embeddings** — broke it. Mean-pooled **ESM-2 (150M)** over the changed region, PCA-50,
gene-grouped CV, full 11k (`esm_embedding_full.py`):

| model | AUC |
|-------|:---:|
| baseline (our 6 features) | 0.787 |
| + ESM region embedding (PCA-50) | 0.816 |
| + canonical−alt **difference** embedding | 0.794 |
| **+ region + difference** | **0.822** (+0.035) |

Nuances: ESM *alone* (0.730) is weaker than our features (they're genuinely good);
the raw 640-dim vector *overfits* (0.782) and must be PCA-compressed; the
fold-change-targeted canonical−alt difference adds on top. This is the first thing
to exceed the ceiling — confirming that the limit was **information**, and pLMs carry
more of it than any scalar. It also answers "run our model on the VEPs' data": DOMAS
already consumes the #1 VEP (AlphaMissense = `region_am`); adding ESM **embeddings**
helps (+0.035), whereas scalar VEP scores (EVE etc.) are redundant with AlphaMissense.
The end-game is **ESMFold** — folding the isoform directly from the same pLM — for the
routed uncertain cases.

**Combined model — regression, calibration, and the fold band**
(`esm_combined_analysis.py`, gene-grouped CV, PCA fit per fold). The joint fit
(our 6 features + ESM region + difference PCA-50, one model) also lifts the
*continuous* TM prediction — Ridge R² **0.360 → 0.461** (+0.10, RMSE 0.178→0.164) —
a bigger relative gain than the AUC, and it stays well-calibrated (predicted prob ≈
observed change rate in every decile, monotone 4%→85%).

**Recommended fold band for real folding (route P(fold change) inside the band):**

| route-band | % routed → fold | acc *in* band | % called | acc called |
|------------|:---:|:---:|:---:|:---:|
| 0.45–0.55 (screening) | 8% | 52% | 92% | 79% |
| **0.40–0.60 (default)** | **16%** | **54%** | **84%** | **81%** |
| 0.30–0.70 (high-stakes) | 32% | 58% | 68% | 85% |

Default = **0.40–0.60**: the rationale is principled — the majority-class baseline is
66% accuracy, and *within this band the model only reaches 54%* (worse than a naive
guess), so those cases genuinely need structure; 16% is a manageable fold budget and
84% of events are still called at 81%. In-band accuracy stays below the 66% baseline
out to ~0.30–0.70, so wider bands are defensible for higher stakes. Pair the band
with the high-id+exposed widened trigger to also catch the confidently-wrong
high-id/low-TM outliers (which sit *outside* the band).

**Hyperparameter sweep — no meaningful headroom.** Grid over PCA dim {20,50,100,150}
× logistic C {0.05…1.0} × Ridge α {1…50}, gene-grouped CV. The classifier AUC is
**dead flat at 0.822** regardless of regularization or PCA dim; regression R² is
**insensitive to α** and gains only from *more* PCA components (0.461 at PCA-50 →
**0.472 at PCA-150**, +0.011, adopted). The total insensitivity to regularization
means the model is neither over- nor under-fitting — we're at the information limit
again, so the remaining lever is a **richer embedding** (ESM-2 650M), not tuning.

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
