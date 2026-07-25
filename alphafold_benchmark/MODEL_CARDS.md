# DOMAS calibrated scores — model cards

DOMAS emits two **calibrated probability** scores alongside the categorical `impact`
level, one per axis of alternative-splicing consequence:

| score | axis | question it answers | gene-grouped AUC |
|-------|------|---------------------|:----------------:|
| **`fold_change_prob`** | structural | Does the alternative isoform adopt a *different 3-D fold*? | **0.894** |
| **`impact_prob`** | functional | Is the changed region *functionally important* (pathogenic-variant-like)? | **0.766** |

They are deliberately different questions ("fold-vs-function duality"): a change can
preserve the fold yet abolish function, or wreck the fold while being biologically
silent. Each score is a logistic regression over a small set of interpretable features,
**natively calibrated** (a predicted 0.6 really means ~60%), reported throughout under
**gene-grouped 5-fold cross-validation** (every isoform of a gene is confined to one
fold, so no protein appears in both train and test — no leakage).

Everything below is reproduced end-to-end by [`reproduce_models.py`](reproduce_models.py)
(see §5).

---

## 1. Shared benchmark, labels, and methodology

**Benchmark dataset.** Miller et al., *Predicting the structural impact of human
alternative splicing*, **Genome Biology 2025** (doi:10.1186/s13059-025-03744-x). One
row = **one canonical-vs-alternative-isoform splicing event**. After dropping pure
insertions (no canonical changed region to score), **10,227 pairs across 5,557 human
proteins**.

**Label sources — what they are and what they mean:**

| label | source | meaning | used by |
|-------|--------|---------|---------|
| **TM-score** | The paper folded every canonical and alternative isoform with **AlphaFold2** and computed the pairwise TM-score (Table S4). | Structural-similarity of two folds, 0–1. **TM < 0.5 = different fold; ≥ 0.5 = same fold.** The structural ground truth. | `fold_change_prob` |
| **humsavar** | UniProt **humsavar** disease-variant catalogue (aggregates ClinVar/OMIM), single-residue variants in **canonical protein coordinates**, classed **LP/P** (pathogenic) or **LB/B** (benign). | A changed region "overlaps a pathogenic/benign variant" if a LP/P / LB/B position falls inside it. | `impact_prob` |

**Evaluation protocol.** L2-regularised logistic regression; features standardised
(z-scored); missing features imputed to the training median. Metrics:
- **AUC** — ROC-AUC (rank/threshold-independent; 0.5 = chance).
- **accuracy** — at decision threshold 0.5.
- **R²** — for `fold_change_prob`, variance of the *continuous* TM-score explained by a
  ridge regressor on the same features; for `impact_prob` (a binary target), the
  variance of the 0/1 label explained by the predicted probability (label-R², lower by
  construction).

**Feature provenance** (each feature, and where its value comes from):

| feature | what it is | source |
|---------|-----------|--------|
| `pae_global` | mean **Predicted Aligned Error** over the whole *canonical* AlphaFold structure — how rigidly/confidently the fold is assembled | AlphaFold DB `predicted_aligned_error_v6.json` → `afdb_pae` |
| `identity` | % sequence identity between canonical and alternative protein, from the trimmed changed region | UniProt reference proteome + varsplic sequences |
| `max_cov_loss` | maximum **Pfam** domain coverage (%) lost in the changed region | HMMER `hmmsearch --cut_ga` vs Pfam-A → `pfam_match` |
| `protL` | canonical protein length (aa) | UniProt sequence |
| `region_am` | mean **AlphaMissense** pathogenicity over the changed region | AlphaMissense per-residue → `am_pathogenicity` |
| `loeuf` | gnomAD **LOEUF** gene loss-intolerance (lower = more constrained) | gnomAD v2.1.1 by-gene → `gene_constraint` |
| `buried_frac` | fraction of the changed region **buried** (RSA < 0.25) in the canonical structure | AlphaFold DB model → Shrake–Rupley RSA → `afdb_rsa` |

---

## 2. `fold_change_prob` — structural axis

**What it measures.** Calibrated P(the alternative isoform adopts a **different fold**,
TM-score < 0.5).

**How it's measured / label.** Trained and evaluated against the paper's AlphaFold2
TM-score (§1). TM < 0.5 is the positive class.

**Features & standardized coefficients** (`utils._FOLD_CHANGE_MODEL`; larger |coef| =
more influence):

| feature | std. coef | direction |
|---------|:---------:|-----------|
| **`pae_global`** | **+1.96** | floppier / multi-domain canonical structure → more fold change (dominant) |
| `identity` | −0.68 | more sequence conserved → fold preserved |
| `max_cov_loss` | +0.65 | more Pfam domain lost → more fold change |
| `protL` | −0.54 | longer protein → fold preserved (larger scaffold) |
| *intercept* | −1.26 | |

*(raw feature medians `[17.13, 75.0, 11.0, 388.0]`, means `[16.59, 66.36, 33.22, 379.18]`,
std `[7.12, 28.02, 39.28, 134.36]` for the four features in order.)*

**Performance.** Base rate P(TM<0.5) = 0.338, so the majority-class accuracy baseline is
0.662.

| operating point | AUC | accuracy | R² | coverage |
|-----------------|:---:|:--------:|:--:|:--------:|
| **all events (no band)** | **0.894** | **0.815** | **0.659** | 100% |
| suggested band applied (route 0.40–0.60 to folding) | **0.913** | **0.854** | **0.686** | 88% called |

**Suggested query band: route `fold_change_prob` ∈ [0.40, 0.60] to real folding**
(ESMFold/ColabFold) instead of forcing a call. Rationale, visible in the decile table
and the right-hand graph below: the middle deciles (7–8, prob ≈ 0.37–0.69) are where the
model degrades to a **coin flip** — in-bin accuracy dips to 0.54, *below* the 0.662
majority baseline — because these are the "does the remainder re-fold?" cases (the SRP9
class) that cheap features provably cannot resolve. Abstaining on that ~12% lifts the
called-set accuracy 0.815 → 0.854 and AUC 0.894 → 0.913. Widen to 0.30–0.70 for
higher-stakes screening.

**Decile calibration** (events sorted by predicted probability into 10 equal bins):

| decile | prob range | mean pred | observed rate | n | in-bin acc@0.5 |
|:------:|:----------:|:---------:|:-------------:|:---:|:--------------:|
| 1 | 0.001–0.013 | 0.007 | 0.004 | 1023 | 0.996 |
| 2 | 0.013–0.032 | 0.022 | 0.008 | 1023 | 0.992 |
| 3 | 0.032–0.067 | 0.048 | 0.040 | 1023 | 0.960 |
| 4 | 0.067–0.133 | 0.097 | 0.073 | 1023 | 0.927 |
| 5 | 0.133–0.230 | 0.178 | 0.185 | 1023 | 0.815 |
| 6 | 0.230–0.367 | 0.296 | 0.325 | 1023 | 0.675 |
| 7 | 0.367–0.524 | 0.443 | 0.490 | 1023 | **0.536** |
| 8 | 0.524–0.685 | 0.606 | 0.605 | 1022 | 0.605 |
| 9 | 0.685–0.847 | 0.763 | 0.758 | 1022 | 0.758 |
| 10 | 0.848–0.991 | 0.917 | 0.890 | 1022 | 0.890 |

Predicted ≈ observed in every decile (well-calibrated). The U-shaped accuracy — high at
both confident tails, coin-flip in the middle — is the signature that motivates the band.

![fold_change_prob calibration and per-decile accuracy](model_card_foldchange.png)

**Caveats.**
- **Confounded mechanism.** Part of `pae_global`'s power is that a TM-score between two
  *low-confidence* AlphaFold models is low almost by construction. So a high canonical
  PAE partly flags "this TM is unreliable" as much as "a real fold change occurred." It
  is still non-circular (canonical-only, computed independently of the label) and
  prospectively available, but read it as "confidence the fold is well-defined enough to
  compare," not a pure fold-change detector.
- **Per-protein prior.** `pae_global` is a whole-canonical-structure property, identical
  for every isoform of a gene — it sets a per-protein baseline; only `identity`,
  `max_cov_loss`, `protL` differentiate isoforms of the same gene.
- **Known blind spot.** High-identity / low-TM outliers (≈10% of pairs: identity ≥ 80%
  yet TM < 0.5 — the SRP9 class) are confidently mis-called *preserved*; these sit
  *outside* the routing band, so pair the band with a "high-identity AND exposed region"
  trigger for high-stakes use.
- **Scope.** Pure insertions are excluded (no canonical changed region). The ceiling for
  cheap features is ~AUC 0.90; the residual genuinely needs the isoform folded.

---

## 3. `impact_prob` — functional axis

**What it measures.** Calibrated P(the changed region is **pathogenic-variant-like**).
Concretely, among changed regions that overlap *some* curated variant, P(it is the
**pathogenic** one rather than the benign one). It is an **enrichment / prioritisation**
signal ("look here first"), **not** a per-event clinical classifier.

**How it's measured / label.** Restricted to regions overlapping a humsavar LP/P or LB/B
variant; label = 1 if pathogenic-overlapping, 0 if only benign-overlapping. 3,276
regions. (The benign class is the specificity control — a real functional score must
enrich pathogenic, not merely variant-dense, regions.)

**Features & standardized coefficients** (`utils` impact model):

| feature | std. coef | direction |
|---------|:---------:|-----------|
| **`region_am`** | **+0.80** | higher AlphaMissense constraint → pathogenic (dominant) |
| `loeuf` | −0.42 | more loss-intolerant gene (lower LOEUF) → pathogenic |
| `buried_frac` | +0.24 | buried / structured region → pathogenic |
| `max_cov_loss` | −0.06 | ~no effect |
| *intercept* | −1.19 | |

*(raw medians `[0.482, 0.983, 41.0, 0.340]`, means `[0.482, 0.979, 46.07, 0.309]`,
std `[0.155, 0.443, 41.01, 0.219]`.)*

**Performance.** Base rate P(pathogenic | overlaps a variant) = 0.277 → majority baseline
0.723.

| operating point | AUC | accuracy | R² | coverage |
|-----------------|:---:|:--------:|:--:|:--------:|
| **all events (no band)** | **0.766** | **0.750** | **0.175** | 100% |
| abstain on 0.40–0.60 | **0.770** | **0.793** | **0.181** | 85% called |

**Suggested use / band.** Treat `impact_prob` as a **triage/prioritisation** score, not a
verdict: the **top decile (prob ≳ 0.58)** is ~66% pathogenic — "review these first"; the
**bottom deciles (prob ≲ 0.11)** are ~5–7% pathogenic — reliably de-prioritised; the
**0.40–0.60** middle is genuinely ambiguous and should go to manual curation. Abstaining
on the middle raises called-set accuracy 0.750 → 0.793 (AUC barely moves — the tails are
where the signal is). Do **not** use it to call an individual variant pathogenic.

**Decile calibration:**

| decile | prob range | mean pred | observed rate | n | in-bin acc@0.5 |
|:------:|:----------:|:---------:|:-------------:|:---:|:--------------:|
| 1 | 0.019–0.068 | 0.047 | 0.043 | 328 | 0.957 |
| 2 | 0.068–0.107 | 0.087 | 0.067 | 328 | 0.933 |
| 3 | 0.107–0.143 | 0.124 | 0.113 | 328 | 0.887 |
| 4 | 0.143–0.184 | 0.163 | 0.174 | 328 | 0.826 |
| 5 | 0.184–0.227 | 0.205 | 0.177 | 328 | 0.823 |
| 6 | 0.227–0.281 | 0.254 | 0.293 | 328 | 0.707 |
| 7 | 0.281–0.345 | 0.311 | 0.339 | 327 | 0.661 |
| 8 | 0.345–0.436 | 0.390 | 0.468 | 327 | 0.532 |
| 9 | 0.438–0.576 | 0.505 | 0.437 | 327 | 0.508 |
| 10 | 0.576–0.911 | 0.683 | 0.661 | 327 | 0.661 |

![impact_prob calibration and per-decile accuracy](model_card_impact.png)

**Caveats.**
- **Semi-circular (the biggest one).** `region_am` (AlphaMissense) is itself a
  pathogenicity predictor and dominates the model (+0.80), so part of the signal is built
  in. The **non-circular** contributors are `loeuf` (gnomAD constraint) and `buried_frac`
  (structure). Read the score as "concentrates known functional signals," not independent
  evidence.
- **Positive-only labels.** "No known pathogenic variant" is not proof of neutrality — so
  these are recall/enrichment results, not true specificity.
- **Region-length confound.** On *raw* variant-overlap presence, region length dominates
  (bigger change → more residues → more overlap); `impact_prob` is trained
  length-controlled (patho-vs-benign among overlapping regions) precisely to avoid
  measuring size instead of function. Do not add sequence-identity/length as a feature —
  it would inflate AUC while measuring region size.
- **Not a classifier.** Both classes pile into moderate/high probabilities; specificity
  is limited. It ranks/prioritises; it does not adjudicate.

---

## 4. Why these features and not others (audit trail)

Both feature sets are the survivors of an exhaustive search (EXPERIMENTS.md, E1–E44).
Key negatives worth knowing for review:
- **PAE dominates fold change** and was the single biggest lever (AUC 0.78 → 0.90); it was
  *assumed* redundant for years and only helped when finally measured.
- **`buried_frac` is functional-only; `pae_global` is structural-only** — each was tested
  on the *other* axis and rejected (drop-column ≈ 0 or fails the benign control). The two
  models are cleanly separated.
- **Rejected / redundant:** ESM-2 embeddings (help fold change but not shipped — heavy),
  650M vs 150M ESM (≈equal), region pLDDT (0.79-correlated with burial), contacts
  (redundant with burial), phyloP (only a small structural add, small sample), peak-AM,
  tree/NN model families, multi-task. None beat the shipped logistic on its own axis.

---

## 5. Reproduction

[`reproduce_models.py`](reproduce_models.py) reproduces every number and both graphs in
this document. It:
1. **Downloads** (into `./repro_data/`, ~30 MB, no DoChaP DB needed): UniProt human
   reference proteome + varsplic (sequences → `identity`, `protL`, changed region),
   UniProt humsavar (functional labels), gnomAD v2.1.1 by-gene constraint (`loeuf`).
2. **Reuses** the committed per-pair feature CSVs in this folder (their raw provenance —
   AlphaFold DB, AlphaMissense, Pfam — is documented in `DATA_SOURCES.md`):
   `full_pairs.csv` (TM label + changed region), `pae_features_full.csv` (`pae_global`),
   `bench_rich_results.csv` (`region_am`, `max_cov_loss`), `full_features.csv`
   (`buried_frac`).
3. Refits both models under gene-grouped 5-fold CV (coefficients match the shipped
   `utils` models), prints overall + banded metrics and decile tables, and writes
   `model_card_foldchange.png` / `model_card_impact.png`.

```
cd alphafold_benchmark && python3 reproduce_models.py
```

---

## 6. Review checklist / essentials

- **Intended use.** Screening/prioritising alternative-splicing events at scale:
  `fold_change_prob` for "does the fold change," `impact_prob` for "is the changed region
  functionally important." Both are triage signals feeding human/structural review.
- **Out of scope.** Not a clinical variant classifier; not validated on non-human
  proteins; not for pure insertions; not a substitute for folding the isoform when the
  score is in the uncertain band.
- **Leakage control.** Gene-grouped CV throughout (all isoforms of a gene in one fold);
  `loeuf` is gene-level so pair-level CV would leak — the gene-grouped numbers are the
  honest ones.
- **Class balance.** fold change: 33.8% positive (majority acc 0.662). functional: 27.7%
  pathogenic among overlapping regions (majority acc 0.723). Always compare model accuracy
  to these baselines, not to 0.5.
- **Calibration.** Both are logistic (natively calibrated); the decile tables show
  predicted ≈ observed in every bin.
- **Known failure modes.** fold change: high-identity/low-TM re-folders (SRP9 class),
  outside the band. functional: both variant classes pile high (limited specificity);
  semi-circular via AlphaMissense.
- **Versioning.** Coefficients live in `code/utils.py` (`_FOLD_CHANGE_MODEL` and the
  impact model), regenerated by `fit_foldchange_pae.py` / `fit_calibrated.py`. Card
  reflects the current shipped models; base benchmark = Genome Biology 2025.
- **Provisioning dependency.** `fold_change_prob` needs the `afdb_pae` table
  (`build.py -source afdb_pae`); if absent, `pae_global` imputes to the training median and
  the score degrades gracefully toward an identity/coverage-only prediction.
