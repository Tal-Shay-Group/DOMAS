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
- **Ceilings.** Cheap features predict fold change to ~AUC 0.75 / R² 0.31; the residual
  ("does the remainder refold" — the SRP9 class) needs the isoform actually folded,
  which is what the TM label already encodes. Functional truth is positive-only, so
  these are recall/enrichment results, not specificity.

## Open threads (not done)

- Cross-species conservation (phyloP via DoChaP CDS→genome mapping) — scoped, not run.
- Fold burial into the production scorer (mind its counter-intuitive direction / the
  function-vs-fold target distinction).
- Integrate ASpdb as a second, independent functional label.
