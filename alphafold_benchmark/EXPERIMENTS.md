# Experiments and conclusions

Chronological log of every analysis run in this investigation, with method and
outcome. The through-line: does DOMAS's impact score capture the structural /
functional effect of alternative splicing? Scripts referenced live in this folder.

---

## Phase 1 — pin down an exact, reproducible benchmark

**E1. Extract exact isoform ground truth.** Parsed the Genome Biology 2025
supplementary (Tables S2/S4) for 7 case-study genes → exact UniProt isoform IDs +
published TM-scores + sequence identity.
→ *Ground truth secured; earlier heuristic isoform matching was replaced by exact IDs.*

**E2. Map UniProt isoforms to DoChaP transcripts by sequence.** Fetched each
isoform + canonical FASTA and matched against DoChaP `ensembl_sequence`.
→ *All 7 pairs matched a DoChaP transcript at 100% identity — the benchmark is
fully reproducible in DoChaP.*

**E3. Whole-protein DOMAS impact on the 7 exact pairs vs TM verdict.**
→ *Sensitivity 3/3 on real structural changes (CFHR3, MDH2, BAX). Specificity gap:
SRP9 (fold preserved, TM 0.76) called HIGH (false positive); CXCR3's N-terminal
insertion missed (NONE). Good recall, poor precision — on 7 cases.*

**E4. Gene PDFs for the 4 divergent genes** (CDK1, CXCR3, SRP9, EEF1AKMT3) with
HMM (Pfam) enrichment tracks + changed-HMM tables.
→ *Rendered correctly once `hmm_by_transcript` was populated and a canonical
reference forced; documented in the gene-PDF workflow.*

**E5. Assess paper bioRxiv 578053 (SEPTIN9, MAPT).**
→ *No TM-scores, no UniProt accessions, no exact isoform pairs anywhere — its case
studies are qualitative (IDR) / disease-missense, not isoform-vs-isoform. Excluded
from the quantitative benchmark.*

## Phase 2 — scale to the full dataset

**E6. Insertion-aware `changed_region` + insertion branch.** Made the divergent
span return canonical+alt coords and a kind; added a net-added-residue impact branch.
→ *CXCR3's insertion stopped reading as a phantom 3-aa change (NONE → MODERATE);
non-insertion cases unchanged. (Later softened — see E14.)*

**E7. Full 11k benchmark: DOMAS impact vs TM verdict.** 11,068 pairs (all Table S4
isoforms with both sequences), one hmmsearch over 16,948 sequences, impact scored
per pair. `run_tm_benchmark.py`.
→ *At impact ≥ moderate: sensitivity 0.79, **specificity 0.29**. Precision 35% vs
33% base rate (**lift 1.08×**). Impact is only weakly predictive of fold change.*

**E8. Sharper cuts.**
→ *Curated real-changes (Table S2): only 68% flagged (the 3/3 pilot was luck).
Clean negatives (TM ≥ 0.95): only 51% correctly quiet. Errors are not near the
TM=0.5 boundary — 1,279 false positives have TM ≥ 0.85; 37% of pairs are off by
≥2 ordinal levels.*

## Phase 3 — is any of the signal real? (identity control)

**E9. Control for sequence identity (partial correlation).** `analyze_identity_control.py`.
→ *Identity↔TM = +0.46; impact↔TM = −0.17 → **+0.05 after controlling identity**.
Within every identity band, impact→TM ≈ 0. **Impact adds nothing beyond sequence
divergence for predicting fold change.** Mean identity by impact: none 86% → high 63%
— impact is essentially an identity proxy.*

**E10. Numeric correlation + error placement.**
→ *impact rank vs TM Spearman −0.19, max-coverage-loss −0.27, identity +0.46. Boxplots
show near-total overlap; even "high" impact has median TM 0.58 (fold preserved).*

**E11. Cross-validated prediction (classification + regression).** `predict_eval.py`,
5-fold CV, 10,232 pairs.
→ *Predicting CHANGE (TM<0.5): DOMAS impact **AUC 0.568** (= majority-class accuracy,
useless); identity 0.697; identity+burial **0.755**. "identity + impact" (0.696) =
identity alone — impact contributes nothing. Regression R²: identity 0.239 →
identity+burial 0.307.*

## Phase 4 — is AlphaFold exhausted? (structural features)

**E12. Burial (RSA) + DSSP of the changed region.** Downloaded canonical AF
structures; computed RSA (Shrake–Rupley) + secondary structure (pydssp) over the
changed region. Pilot (7) → sample (2,085) → **full 11k** (`extract_burial_dssp_full.py`).
→ *Burial **survives identity control**: partial corr **+0.30**, non-zero in every
identity band; incremental R² 0.24 → **0.31**; not a region-length artifact. First
identity-orthogonal signal found. **Counter-intuitive direction**: exposed changed
regions → lower TM (big alt exons sit in exposed loops/termini); buried changes
superpose. DSSP secondary structure adds little (+0.17 partial).*

**E13. Contact modularity.** Cβ<8 Å contacts from the changed region to the rest.
`extract_contacts.py`.
→ *Correlates with TM alone (~+0.47) but **redundant with burial**: partial drops to
+0.10 after burial; incremental R² only +0.006. Burial is the compact summary of the
canonical-structure geometry; contacts (and, by extension, likely PAE) add little.*

## Phase 5 — calibrate the scorer

**E14. Soften the insertion branch.** Analysis showed insertions (7.5% of pairs) have
a 20% true-change rate but were flagged ≥moderate 53% of the time. Reworked
`insertion_impact` to escalate on domain placement, not size.
→ *Insertion flagging 53% → **20%** (matches truth); overall **specificity 0.295 →
0.335** for −0.02 sensitivity. CXCR3's exposed +47 insertion: MODERATE → LOW;
non-insertion cases unchanged.*

**E15. Effect of excluding insertions from the evaluation.**
→ *Negligible — every metric moves ≤0.02. DOMAS's weakness is systemic across
substitutions/deletions, not insertion-specific. (But confirmed insertions are an
over-called subclass, motivating E14.)*

## Phase 6 — the functional axis (the payoff)

**E16. Pathogenic-variant overlap enrichment.** Label = changed region overlaps a
UniProt humsavar LP/P variant (canonical coords). `clinvar_enrichment.py`. Restricted
to the 998 disease proteins (2,023 pairs; 907 positives).
→ *Impact tracks pathogenic overlap: none 8.5% → high 63.5%. **Survives length
control** (within-band 6×). Multivariate over length+identity: impact **+0.016 AUC**,
burial **+0.022**. Unlike on TM, impact adds real functional signal.*

**E17. Benign contrast (specificity check).** Same test with LB/B variants.
→ *Decisive: at fixed length impact separates pathogenic 6× (0.05→0.31) but is
**flat for benign** (0.17→0.17). Multivariate: impact +0.016 (patho) vs **−0.001**
(benign); burial +0.022 vs **+0.000**. Both features are **pathogenicity-specific** —
DOMAS finds functionally important regions, not merely variant-dense ones.*

**E18. Per-pair contingency table.** `variant_impact_contingency.py`.
→ *75% of pathogenic-overlapping splicing events are high-impact (vs 3% none); the
pathogenic share within a column rises none 7% → high 28%, benign stays diffuse.*

## Phase 7 — fold burial into the scorer, and its limits

**E19. Incorporate burial into the impact score.** Added `afdb_rsa` (per-residue
RSA) + `Enricher.rsa`; `hmm_change_impact` weights by burial of the changed region
(buried < 0.30 RSA → raise, exposed > 0.50 → lower). Tried (A) burial as a fallback
below AlphaMissense, then (B) a fall-through ladder where an inconclusive AM yields
to burial.
→ *Direction is function-correct (buried regions overlap pathogenic variants more:
buried_frac 0.38 vs 0.21). Effect on 11k: exposed changes downgraded
(moderate→none). On disease proteins, B cleans the high bucket (benign high −47,
pathogenic +3). **But the predictive gain is marginal (+0.001 AUC)** because
AlphaMissense — top of the ladder, 92% coverage — is a correlated
functional-constraint signal that already carries burial's information (burial is
to AM what contacts were to burial: largely redundant). Kept B; burial's
non-redundant value is the ~8% of pairs lacking AM, plus being non-circular.*

**E20. Enrichment signal vs. classifier (50:50 prior).** Asked the harder question:
given a region overlaps *some* variant, can impact call pathogenic vs benign?
→ *No — near chance. Both classes pile into moderate/high (P(mod-or-high|patho)=95%,
|benign=84%), specificity ~16%. Under a balanced prior, "impact ≥ moderate →
pathogenic" has precision 53% (wrong 47%); "= high" 55%. The length-controlled
enrichment (E16–E17) was largely region-size; conditioning on variant-presence and
balancing classes removes it. **DOMAS impact is a coarse enrichment/prioritisation
signal, never a per-event pathogenic-vs-benign classifier.***

**E21. Peak AM + gnomAD prototype.** Tried the two proposed features on
patho-vs-benign discrimination. (`improve_proto.py`.)
→ *Peak AM **rejected** — for domain-scale regions (median 169 aa) mean AM (0.754)
beats peak (0.703) and every mean-of-top-K (0.681–0.703) monotonically;
pathogenicity is about overall regional constraint, not a hotspot. gnomAD **LOEUF**
is a real, non-circular add (0.754 → 0.770). Biggest finding: the categorical
impact (AUC 0.585) discards most of its own signal vs the continuous mean AM
(0.754) — **quantisation is the bug.***

**E22. Build the calibrated continuous score.** Logistic over region_am + LOEUF +
max_cov_loss + buried_frac → `impact_prob`; provisioned gnomAD (`gene_constraint`,
`Enricher.loeuf`, `build_gnomad`), shipped constants in `utils.impact_probability`.
(`fit_calibrated.py`, `calib_model.json`.)
→ *5-fold CV AUC **0.769**, well-calibrated, deciles monotonic with top decile 58%
pathogenic (vs categorical 38%). Wired into add_scores as a companion column.
region_am dominates (+0.80 std coef), LOEUF the non-circular add (−0.42), burial
+0.24, coverage-loss ~0.*

**E23. Rigorous evaluation — gene-grouped CV.** Checked whether pair-level CV leaks
(a protein's isoforms split across train/test; LOEUF is gene-level).
→ *Gene-grouped AUC **0.765** vs pair-level 0.768 (−0.003): negligible leakage
(~1.7 isoforms/protein). The ~0.77 is honest. Reported number switched to the
gene-grouped one.*

**E24. Second calibrated score — `fold_change_prob` (structural axis).** Same
approach applied to the TM target: logistic over identity + burial + region_am +
LOEUF + max_cov_loss → P(TM<0.5). (`fit_foldchange.py`, `foldchange_model.json`,
`utils.fold_change_probability`.) Runtime identity from the trimmed changed region.
→ *Gene-grouped AUC **0.777**, accuracy 74% (@0.5), well-calibrated. AM helps only
**in combination** (suppressor: 0.54 alone, 0.70 with identity, 0.78 with
identity+burial). **burial flips sign** vs impact_prob (−0.93 structural vs +0.24
functional) — the two scores answer different questions. Wired into add_scores and
shown on the gene PDF with a sign-flip footnote.*

**E25. Model-family comparison.** logistic vs +interactions vs random forest vs
gradient boosting vs hist-GB vs MLP neural net, gene-grouped CV, both targets.
(`ml_compare.py`.)
→ *Logistic **wins** functional (0.761; trees/NN 0.72–0.75) and loses only ~0.01
structural (RF 0.798 vs logistic 0.787). With 4–6 monotone tabular features, tree
ensembles and NNs add nothing; logistic kept for interpretability + native
calibration. The ceiling is set by features/labels, not the algorithm. (Bayesian
would add per-prediction uncertainty, not accuracy.)*

**E26. Exhaust the feature space (length + all pairwise + squares, L1-pruned).**
Expanded both models to degree-2 (main effects + length + squares + all pairwise
products), L1-pruned, gene-grouped CV.
→ *Length adds ~0. Full expansion: no gain functional (0.763→0.759), +0.007
structural (0.788→0.795). **L1 pruned almost nothing** (19/20, 33/35 terms survive
with tiny coefficients) — no clean sparse improvement hiding, no dead parameters.
phyloP conservation (the one untested cheap feature) assessed as **very likely
redundant** (AlphaMissense is trained on conservation; feature space saturated;
constraint axis already covered by region_am + LOEUF) and **not built** — the
protein→genome→bigWig pipeline (Ensembl REST down, remote phyloP slow) was not
justified. **Feature question closed: the ceiling ~0.76 functional / ~0.79
structural is set by the information in the features, not the model or feature
engineering.** The only way up is orthogonal data — folding + better labels.*

## Phase 8 — uncertainty routing (b) and its limits

**E27. Uncertainty routing / selective prediction.** `utils.fold_change_call` maps
`fold_change_prob` → change / preserved / **uncertain**; abstain on uncertain and
route to folding. (`selective_prediction.py`.)
→ *Abstaining lifts accuracy on the called set 75%→86% (structural, 50% coverage).
The routed middle is a **coin flip** (acc ~53-66% in-band, observed change rate
0.47) — genuinely ambiguous, not error. Uncertain cases (17% at 0.4-0.6) have the
SRP9 signature: divergent, exposed, high cov-loss. Band tables for both models
recorded.*

**E28. Widened trigger for the confidently-wrong outliers.** The high-id/low-TM
outliers (identity≥80 AND TM<0.5; 1,035 = 10%) are *confidently* mis-called
`preserved`, so the confidence gate catches only 17%.
→ *Domain-loss widening barely helps (~30%); **"high identity AND exposed"
(mean_rsa≥0.5) catches ~73%** (the outliers are exposed fold-flipping changes — the
burial signal), but doubles the fold budget (17%→39%), precision ~19%. A fix for
high-stakes use only.*

**E29. Specialist model on the high-id regime — redundant.** Tried a second model
trained only on high-identity pairs (identity dropped, since ~constant there) to
crack the outlier regime.
→ *No gain: specialist AUC **0.730** vs global **0.731** on the same subset. Within
high-id, burial is the signal (single-feature ~0.65) but the **global logistic
already uses it** — a near-constant feature doesn't corrupt the others, so the
global model *is* effectively the specialist, and E26 already showed no
identity×burial interaction to exploit. The high-id regime tops out at AUC ~0.73
because the **signal is weak, not the model**. Only folding breaks it.*

## Phase 9 — protein language model embeddings (the first ceiling break)

**E30. ESM-2 embeddings of the changed region.** After hand-crafted features,
feature engineering, model families, cascades and specialists all stalled at ~0.79,
tested the one orthogonal data type: pLM embeddings. Mean-pooled ESM-2 (150M) over
the changed region (canonical) and the canonical−alt difference; PCA-50; gene-grouped
CV. (`esm_embedding_test.py`, `esm_embedding_full.py`.)
→ *Full 11k: baseline 0.787 → **+region-embed 0.816** → **+region+difference 0.822**
(+0.035) — the **first feature to exceed the ceiling**. Nuances: ESM *alone* (0.730)
is weaker than our features (they're good); raw 640-dim overfits (0.782) — must
PCA-compress; the fold-change-targeted canonical−alt difference helps on top of the
region embedding. Confirms the literature: pLM embeddings carry information our
scalars don't. The ceiling was information, and pLMs supply more of it.*

**E31. VEP comparison — we are a region-level meta-model.** Compared to the three
leading variant-effect predictors (AlphaMissense, ESM1v, EVE).
→ *DOMAS operates one level up: VEPs score single missense variants; we score whole
splicing events (regions), on two axes (function + fold change). We **consume the #1
VEP** — AlphaMissense is our top feature (`region_am`). "Running our model on their
data" = already done for AM; done for ESM (E30, embeddings add +0.035); EVE declined
(scalar conservation VEP, worse than & correlated with AM, ~half coverage, gated
behind a JS app). Scalar VEP scores are redundant with AM; only the rich ESM
**embedding** adds.*

**E32. Hyperparameter sweep on the combined model — no headroom.** Grid PCA dim
{20,50,100,150} × logistic C {0.05–1.0} × Ridge α {1–50}, gene-grouped CV.
→ *Classifier AUC **flat at 0.822** across everything; regression R² insensitive to
α, gains only from more PCA components (0.461@50 → **0.472@150**, +0.011). Total
insensitivity to regularization = neither over- nor under-fitting = information
limit again. Remaining lever is a richer embedding (ESM-2 650M), not tuning.
**ESM-2 650M NOT tested — the 2.5 GB weight download from `dl.fbaipublicfiles.com`
stalled dead (~423 MB, 0 KB/s).** Retry in a fresh session via `esm650.py` (3k-subset
150M-vs-650M comparison; needs the model download to succeed).*

---

## Overall conclusion

- **DOMAS impact is NOT a fold-change predictor.** Against TM-score it adds nothing
  beyond sequence identity (partial ≈ 0; AUC 0.568). The earlier "sensitivity" was
  identity leaking through.
- **DOMAS impact IS a functional-disruption signal.** Against disease-variant overlap
  it adds real, **pathogenicity-specific** signal (enriches pathogenic, not benign).
  It was being graded on the wrong exam (fold) when its job is function.
- **AlphaFold burial is the one structural feature worth adding** — identity-orthogonal
  for fold change (partial +0.30; AUC 0.697 → 0.755) and non-circular for the
  functional axis (+0.022, pathogenicity-specific).
- **Ceilings.** Cheap features predict fold change to ~AUC 0.79 (calibrated model) /
  R² 0.31; the residual ("does the remainder refold" — the SRP9 class) needs the
  isoform actually folded, which is what the TM label already encodes. Functional
  truth is positive-only, so these are recall/enrichment results, not specificity.
- **Two shipped calibrated scores.** DOMAS now emits `impact_prob` (functional /
  pathogenicity, gene-grouped AUC 0.765) and `fold_change_prob` (structural /
  P(TM<0.5), 0.777) alongside the categorical `impact`, both logistic (chosen
  empirically over trees/NN) and well-calibrated. Burial enters them with opposite
  signs — the core duality (buried = fold-preserving but pathogenic-variant-rich).
- **The improvement lever is features/labels, not the model.** Un-quantising
  recovered the biggest gain (0.585 → 0.75); gnomAD LOEUF and burial add small
  orthogonal signal; peak-AM and fancier ML add ~nothing. Next real gains need
  cross-species conservation and/or better functional labels.

## Open threads (not done)

- Cross-species conservation (phyloP via DoChaP CDS→genome mapping) — scoped, not run.
- Fold burial into the production scorer (mind its counter-intuitive direction / the
  function-vs-fold target distinction).
- Integrate ASpdb as a second, independent functional label.
