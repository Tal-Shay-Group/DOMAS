# DOMAS calibrated scores — model cards

DOMAS emits two **calibrated probability** scores alongside the categorical `impact`
level, one per axis of alternative-splicing consequence:

| score | axis | question it answers | gene-grouped AUC |
|-------|------|---------------------|:----------------:|
| **`fold_change_prob`** | structural | Is the alternative isoform *structurally divergent* (TM < 0.5)? | **0.865** |
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

## Revision — audit of both scores (E44 / E45)

Both models were re-audited from the committed CSVs; both headline AUCs reproduced
exactly (0.894 all-pairs fold, 0.766 impact), so the numbers below were sound *as
computed*. Two findings changed what ships. Scripts:
[`fit_foldchange_hq.py`](fit_foldchange_hq.py), [`am_alone.py`](am_alone.py).

**1. `fold_change_prob` was refit on AlphaFold-confident structures only, and its band
re-derived.** 32% of the 10,227 training pairs are isoform structures Song et al. class
`low`/`unstructured` and exclude from their own metric analyses. Their TM<0.5 rate is
**0.74** against **0.15** in the confident rows, and isoform pLDDT correlates with TM at
**ρ +0.755** — those rows largely label AlphaFold's failure to model the isoform, not a
splicing-induced fold change. Including them inflated the intercept and broke
calibration where it matters: the old model predicted a mean **0.219** on confident
pairs whose observed rate is **0.150**, its top quintile there predicted **0.68 vs an
observed 0.46**, and its `change` calls were only **54.8%** correct — a coin flip
labelled "confident".

| | all-pairs fit (superseded) | confident-only fit (**shipped**) |
|---|:---:|:---:|
| training rows | 10,227 | 6,973 |
| base rate P(TM<0.5) | 0.338 | 0.150 |
| gene-grouped AUC | 0.894 | **0.865** |
| accuracy / R² | 0.815 / 0.659 | **0.871 / 0.578** |
| `pae_global` coef | +1.96 | **+1.32** |
| intercept | −1.26 | **−2.51** |
| top-quintile calibration | 0.68 pred vs 0.46 obs | **0.48 pred vs 0.45 obs** |
| band | 0.40 / 0.60 | **0.20 / 0.50** |
| `change` precision / recall | 0.548 / — | **0.642 / 0.316** |
| `preserved` NPV | 0.906 | **0.938** |

The band is now asymmetric because the calibrated base rate is 0.150, not 0.5. It gives
called-set accuracy **0.911 at 83% coverage**; 380 of the old model's 894 `change` calls
correctly become `uncertain`. Superseded coefficients: `foldchange_model_allpairs.json`.

*Regime guard (implemented).* The filter uses the **isoform** structure's quality, which
DOMAS never has — it does not fold isoforms, and AlphaFold DB holds one model per
canonical accession (no `P12345-2` entries; verified: 0 of 20,504 `afdb_plddt` rows are
isoform-suffixed). Canonical mean pLDDT is the proxy, separating the regimes at AUC
**0.873**, and is now applied in `utils.fold_change_call` / `fold_change_note`: below
mean pLDDT **70** (AlphaFold's own confidence cutoff) a `preserved` verdict is downgraded
to `uncertain`, and every such row carries `fold_change_note` = "canonical structure
poorly resolved (mean pLDDT NN)".

The guard is **one-sided**, which the numbers justify: `preserved` asserts nothing
happened and is only as good as the structure beneath it (39.7% wrong below the floor vs
6.1% above), while `change` rests on sequence-level facts and stays 0.946 precise there.
Guarding both sides gives the same `preserved` gain but cuts `change` calls from 1,664
to 492 and their precision from 0.851 to 0.744 — pure loss.

| | no guard | **one-sided guard** | both sides |
|---|:---:|:---:|:---:|
| `preserved` calls / error | 6,008 / 0.100 | **5,531 / 0.084** | 5,531 / 0.084 |
| `change` calls / precision | 1,664 / 0.851 | **1,664 / 0.851** | 492 / 0.744 |
| abstention | 25% | **30%** | 41% |

It fires on 2,885 rows — 9% of group A, 70% of group B. **It reduces exposure, it does
not fix the regime:** the group-B `preserved` calls that survive are no more accurate
than before (39.7%). Rows carrying the note need a real isoform structure (ESMFold, E42)
before any structural claim is made about them.

**2. `impact_prob` adds +0.013 AUC over raw AlphaMissense.** Ranking the same 3,282
regions by mean `region_am` with *no model at all* gives **0.753** against the model's
**0.766**; gene-clustered bootstrap 95% CI **[+0.001, +0.024]**. Real but marginal — the
score is region-averaged AlphaMissense with a small correction. Nothing leaks
(AlphaMissense is not trained on clinical labels), but its thresholds were calibrated on
ClinVar, which humsavar aggregates. **The AlphaMissense-free combination
`[loeuf, buried_frac, max_cov_loss]` reaches 0.715** — that is the figure to quote for
signal DOMAS contributes independently. Note `max_cov_loss` alone scores **0.523**, i.e.
chance: Pfam coverage change, the one thing DOMAS measures natively, carries no
functional signal against this label.

**3. What the TM<0.5 label actually is.** It fires on three distinct things the model
cannot separate: genuine refolding; rigid-body **re-orientation** of domains that each
keep their fold (the paper's own explanation for most high-identity/low-TM outliers);
and AlphaFold failure on the isoform. TM is also averaged over both length
normalisations, so deletions depress it mechanically — length ratio alone predicts
TM<0.5 at AUC **0.676**. Read the score as *structurally divergent*, not *refolded*.

**4. The two axes are genuinely independent** — ρ = **−0.05** between the scores, and
−0.07 between `impact_prob` and TM. The fold-vs-function duality is empirical, not
rhetorical. Within the *functional* columns of `results.csv`, however, `impact`,
`impact_prob`, `region_am_mean` and `functional_sites` are largely one measurement
displayed four ways (see `enrichment.add_scores`) — their agreement is not
corroboration.

> Sections 1.1, 2 and 3 have been regenerated against the shipped models by
> [`reproduce_models.py`](reproduce_models.py). **Section 4 and Appendix A still quote
> all-pairs-era AUC deltas** (e.g. "PAE lifted 0.78 → 0.90") — those are the numbers
> that drove the feature-selection decisions, and are kept as the historical audit
> trail; the relative conclusions hold, the absolute values are pre-refit.

**Contents**

1. Shared benchmark, labels, and methodology — *incl. 1.1 the two training/evaluation methods*
2. `fold_change_prob` — structural axis
3. `impact_prob` — functional axis
4. Why these features and not others (audit trail)
5. Reproduction
6. Review checklist / essentials
- Appendix A — Every feature evaluated, and where it comes from
- Appendix B — Comparable published methods

---

## 1. Shared benchmark, labels, and methodology

**Benchmark dataset.** Song et al., *Predicting the structural impact of human
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

### 1.1 Two training / evaluation methods (and their agreement)

The script evaluates each model two ways so the honesty of the headline numbers can be
checked:

- **Method 1 — gene-grouped 5-fold CV** (the reported/shipped method). All isoforms of
  a gene are confined to a single fold, so no protein appears in both train and test.
  This is the leakage-free estimate and matters because a gene-level feature (`loeuf`)
  and correlated same-gene isoforms could otherwise leak.
- **Method 2 — random 80/20 hold-out** (pair-level). A single random 80% train / 20%
  test split, ignoring gene membership. This is the naive method; if it scored much
  higher than Method 1 it would reveal same-gene leakage.

**They agree closely — the small gap confirms leakage is minimal:**

| model | metric | Method 1 (gene-grouped CV) | Method 2 (random 80/20, mean of 5 seeds) | Δ (holdout − grouped) |
|-------|--------|:--------------------------:|:-----------------------:|:---------------------:|
| `fold_change_prob` | AUC | 0.865 | 0.863 | −0.002 |
| | accuracy | 0.871 | 0.870 | −0.001 |
| | R² | 0.578 | 0.581 | +0.003 |
| `impact_prob` | AUC | 0.766 | 0.766 | 0.000 |
| | accuracy | 0.748 | 0.752 | +0.004 |
| | R² | 0.170 | 0.169 | −0.001 |

Method 2 is averaged over seeds 0–4 because a single 80/20 split is noisy at this size:
`impact_prob`'s hold-out AUC ranges 0.740–0.814 across those seeds (n_test = 656), and
`fold_change_prob`'s 0.850–0.885. The earlier single-seed figures in this table sat at
one end of that spread; the averaged gaps are ≤ 0.004 throughout.

The random split is *slightly* optimistic on AUC/R² (≤ +0.018), the expected mild effect
of same-gene pairs landing in both train and test; accuracy even drops for `impact_prob`
(single-split noise on 655 test rows). All gaps are ≤ 0.02, so the two methods give
essentially the same answer — **the gene-grouped numbers are honest, not inflated by
leakage.** They are reported as the headline because they are the conservative estimate;
Method 2 is a single-split point estimate (higher variance) shown only for contrast.

---

## 2. `fold_change_prob` — structural axis

**What it measures.** Calibrated P(the alternative isoform is **structurally divergent**
from the canonical, operationalised as TM-score < 0.5). Not "the domain refolds" — see
§Revision point 3 for the three things that label conflates.

**How it's measured / label.** Trained and evaluated against the paper's AlphaFold2
TM-score (§1), restricted to the 6,973 pairs whose isoform structure the paper classes
`high` or `confident` — the same quality filter it applies to its own analyses. TM < 0.5
is the positive class.

**Features & standardized coefficients** (`utils._FOLD_CHANGE_MODEL`; larger |coef| =
more influence):

| feature | std. coef | direction |
|---------|:---------:|-----------|
| **`pae_global`** | **+1.32** | floppier / multi-domain canonical structure → more fold change (dominant) |
| `identity` | −0.58 | more sequence conserved → fold preserved |
| `max_cov_loss` | +0.52 | more Pfam domain lost → more fold change |
| `protL` | −0.51 | longer protein → fold preserved (larger scaffold) |
| *intercept* | −2.51 | |

*(raw feature medians `[13.81, 78.16, 10.0, 388.0]`, means `[14.00, 69.34, 31.03, 381.62]`,
std `[6.22, 26.47, 37.82, 128.61]` for the four features in order.)*

**Performance.** Fit on the 6,973 AlphaFold-confident pairs (see §Revision). Base rate
P(TM<0.5) = 0.150, so the majority-class accuracy baseline is 0.850 — accuracy is a weak
metric at this imbalance and AUC is the one to read.

| operating point | AUC | accuracy | R² | coverage |
|-----------------|:---:|:--------:|:--:|:--------:|
| **all events (no band)** | **0.865** | **0.871** | **0.578** | 100% |
| suggested band applied (route 0.20–0.50 to folding) | **0.891** | **0.911** | **0.582** | 83% called |

**Suggested query band: route `fold_change_prob` ∈ [0.20, 0.50] to real folding**
(ESMFold/ColabFold) instead of forcing a call. The band is **asymmetric** because the
calibrated base rate is 0.150, not 0.5 — a 0.5 reading is already ~3× the average risk,
so it belongs on the `change` side, and the ambiguous region sits below it. Rationale,
visible in the decile table below: the top two deciles (prob ≳ 0.25) are where in-bin
accuracy falls to 0.63–0.67, because these are the "does the remainder re-fold?" cases
(the SRP9 class) that cheap features provably cannot resolve. Abstaining on that 17%
lifts called-set accuracy 0.871 → 0.911 and AUC 0.865 → 0.891, and gives `preserved`
an NPV of 0.938 with `change` precision 0.642 at recall 0.316. Widen to 0.15–0.45 for a
stricter called set (same accuracy, 78% coverage).

For comparison, the superseded all-pairs model with its 0.40/0.60 band made 894 `change`
calls on these same confident rows at **0.548 precision** — a coin flip labelled
"confident"; 380 of those now correctly land in `uncertain`.

**Decile calibration** (events sorted by predicted probability into 10 equal bins):

| decile | prob range | mean pred | observed rate | n | in-bin acc@0.5 |
|:------:|:----------:|:---------:|:-------------:|:---:|:--------------:|
| 1 | 0.001–0.010 | 0.006 | 0.001 | 698 | 0.999 |
| 2 | 0.010–0.018 | 0.014 | 0.006 | 698 | 0.994 |
| 3 | 0.018–0.029 | 0.023 | 0.010 | 698 | 0.990 |
| 4 | 0.029–0.044 | 0.037 | 0.022 | 697 | 0.978 |
| 5 | 0.044–0.067 | 0.055 | 0.042 | 697 | 0.958 |
| 6 | 0.067–0.103 | 0.084 | 0.088 | 697 | 0.912 |
| 7 | 0.103–0.158 | 0.129 | 0.168 | 697 | 0.832 |
| 8 | 0.158–0.248 | 0.198 | 0.253 | 697 | 0.747 |
| 9 | 0.248–0.440 | 0.335 | 0.329 | 697 | 0.671 |
| 10 | 0.440–0.940 | 0.617 | 0.581 | 697 | **0.628** |

Predicted ≈ observed in every decile (well-calibrated: the largest gap is decile 8, 0.198
predicted vs 0.253 observed). Accuracy now falls **monotonically** with the predicted
probability rather than dipping in the middle as it did pre-refit — at a 0.150 base rate
the low deciles are near-certain `preserved` while the doubt concentrates at the top,
where decile 10 spans 0.44–0.94 and is only 0.628 accurate. That is the shape the band
follows, and why it is asymmetric: the ambiguity to route out sits *below* the `change`
threshold, not symmetrically around 0.5.

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
variant; label = 1 if pathogenic-overlapping, 0 if only benign-overlapping. 3,282
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

**Performance.** Base rate P(pathogenic | overlaps a variant) = 0.276 → majority baseline
0.724.

| operating point | AUC | accuracy | R² | coverage |
|-----------------|:---:|:--------:|:--:|:--------:|
| **all events (no band)** | **0.766** | **0.748** | **0.170** | 100% |
| abstain on 0.40–0.60 | **0.774** | **0.788** | **0.177** | 85% called |

**Suggested use / band.** Treat `impact_prob` as a **triage/prioritisation** score, not a
verdict: the **top decile (prob ≳ 0.58)** is ~66% pathogenic — "review these first"; the
**bottom deciles (prob ≲ 0.11)** are ~5–7% pathogenic — reliably de-prioritised; the
**0.40–0.60** middle is genuinely ambiguous and should go to manual curation. Abstaining
on the middle raises called-set accuracy 0.750 → 0.793 (AUC barely moves — the tails are
where the signal is). Do **not** use it to call an individual variant pathogenic.

**Decile calibration:**

| decile | prob range | mean pred | observed rate | n | in-bin acc@0.5 |
|:------:|:----------:|:---------:|:-------------:|:---:|:--------------:|
| 1 | 0.018–0.069 | 0.047 | 0.040 | 329 | 0.960 |
| 2 | 0.069–0.105 | 0.087 | 0.073 | 329 | 0.927 |
| 3 | 0.105–0.143 | 0.123 | 0.088 | 328 | 0.912 |
| 4 | 0.143–0.185 | 0.163 | 0.192 | 328 | 0.808 |
| 5 | 0.185–0.228 | 0.207 | 0.183 | 328 | 0.817 |
| 6 | 0.228–0.280 | 0.254 | 0.284 | 328 | 0.716 |
| 7 | 0.280–0.344 | 0.311 | 0.338 | 328 | 0.662 |
| 8 | 0.345–0.440 | 0.388 | 0.460 | 328 | 0.540 |
| 9 | 0.440–0.577 | 0.504 | 0.451 | 328 | **0.485** |
| 10 | 0.577–0.885 | 0.681 | 0.655 | 328 | 0.655 |

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
3. **Trains + tests** both models by both methods of §1.1 (gene-grouped CV and random
   80/20), prints coefficients (which match the shipped `utils` models), overall +
   banded metrics, the method comparison, and the decile tables, and writes
   `model_card_foldchange.png` / `model_card_impact.png`.

```
cd alphafold_benchmark && python3 reproduce_models.py            # everything, both models, both methods
python3 reproduce_models.py --download-only                      # just fetch inputs into repro_data/ and exit
python3 reproduce_models.py --models fold_change --methods grouped   # one model, one method
python3 reproduce_models.py --methods holdout --seed 7            # only the random 80/20 split, seed 7
python3 reproduce_models.py --no-graphs                           # skip the PNGs
```

Flags: `--models {fold_change,impact,both}`, `--methods {grouped,holdout,both}`,
`--download-only`, `--no-graphs`, `--seed N` (hold-out split seed).

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

---

## Appendix A — Every feature evaluated, and where it comes from

The two shipped feature sets are the survivors of the full investigation (EXPERIMENTS.md,
E1–E44). This appendix catalogues **all** features that were built and tested — shipped
and rejected — grouped by data source, with what each is and why it did or didn't make
the cut. "S" = in `fold_change_prob` (structural); "F" = in `impact_prob` (functional).

**Sequence-derived** (UniProt reference proteome + `uniprot_sprot_varsplic`):

| feature | what it is | status |
|---------|-----------|--------|
| `identity` | % sequence identity between canonical and alternative protein, from the trimmed changed region | **shipped (S)** — #2 structural feature |
| `protL` | canonical protein length (aa) | **shipped (S)** — small real add |
| `reglen` | length (aa) of the changed region | dropped — it is the region-length confound; deliberately excluded (would measure size, not function) |

**Canonical-AlphaFold-structure-derived** (AlphaFold DB model + its confidence files):

| feature | what it is | source | status |
|---------|-----------|--------|--------|
| `pae_global` | mean Predicted Aligned Error over the whole canonical structure (fold rigidity/confidence) | AFDB `predicted_aligned_error_v6.json` → `afdb_pae` | **shipped (S)** — dominant feature (the biggest lever, AUC 0.78→0.90) |
| `buried_frac` | fraction of the changed region buried (RSA < 0.25) | AFDB model → Shrake–Rupley RSA → `afdb_rsa` | **shipped (F)** — functional (buried → pathogenic-rich) |
| `mean_rsa` | mean relative solvent accessibility of the changed region | same as above | dropped — near-duplicate of `buried_frac` |
| `pae_reg2rest` | mean PAE between the changed region and the rest of the fold (region-specific rigidity) | AFDB PAE | dropped — strong alone (AUC 0.73) but redundant with `pae_global` in the full model |
| `region_plddt` | mean per-residue pLDDT of the changed region (local order/confidence) | AFDB model b-factors → `afdb_plddt` | dropped — 0.79-correlated with `buried_frac`; +0.001 functional / +0.002 structural |
| `helix_frac` / `strand_frac` / `structured_frac` | DSSP secondary-structure fractions of the changed region | AFDB model → pydssp | dropped — small, subsumed by burial/PAE |
| contact density | Cβ<8 Å contacts from the changed region to the rest of the fold | AFDB model coordinates | dropped — redundant with burial (partial corr collapses after burial) |
| `dist_term` / `rel_center` / `frac_before` / `is_terminal` | position of the changed region relative to the protein termini | changed-region coords + `protL` | dropped — +0.017 structural *before* PAE, but PAE absorbs it (+0.002 after) |

**Isoform-AlphaFold-structure-derived** (from Song et al. Genome Biology 2025 Table S4 —
computed from the *isoform's own* AF2 fold, hence **circular** for predicting TM):

| feature | what it is | status |
|---------|-----------|--------|
| iso helix/sheet/loop, radius of gyration, surface charge, IDR % | global structural descriptors of the alternative isoform's predicted fold | dropped — predict TM at AUC 0.92 but **circular** (need the isoform folded, i.e. already have TM) |
| PTM-change counts (`ptm_added`, buried↔exposed, `ptm_missed`) | how the isoform gains/loses/relocates canonical PTM sites | dropped — non-circular but +0.002 functional (redundant with `region_am`/length) |
| charge / radius deltas | canonical-vs-isoform change in charge / radius of gyration | dropped — ~0 contribution |

**Protein language model** (ESM-2, run over the sequences):

| feature | what it is | status |
|---------|-----------|--------|
| ESM-2 150M / 650M region + (canonical−alt) difference embeddings (mean-pooled, PCA-50) | learned sequence representation of the changed region | tested — *breaks the structural ceiling* (AUC 0.79→0.82) but heavy to provision; **not shipped**. 650M ≈ 150M. |
| ESM multi-pool (mean+max+std+first+last token) | positional pooling of per-residue embeddings | dropped — equal to mean-pooling |

**Variant / constraint databases:**

| feature | what it is | source | status |
|---------|-----------|--------|--------|
| `region_am` | mean AlphaMissense pathogenicity over the changed region | AlphaMissense → `am_pathogenicity` | **shipped (F)** — dominant functional feature |
| `loeuf` | gnomAD gene loss-intolerance (lower = more constrained) | gnomAD v2.1.1 by-gene → `gene_constraint` | **shipped (F)** — non-circular constraint signal |
| `pLI` | gnomAD probability of LoF-intolerance | gnomAD | dropped — weaker than LOEUF |
| peak-AM / mean-of-top-K AM | hotspot AlphaMissense instead of regional mean | AlphaMissense | dropped — regional mean beats every peak (domain-scale regions) |
| `func_site` | whether the region overlaps a curated UniProt functional residue | UniProt features | used in the *categorical* `impact` only (semi-circular for the probability) |

**Domain / coverage:**

| feature | what it is | source | status |
|---------|-----------|--------|--------|
| `max_cov_loss` | maximum Pfam domain coverage (%) lost in the changed region | HMMER `hmmsearch --cut_ga` vs Pfam-A → `pfam_match` | **shipped (S and F)** — the one feature in both models (opposite signs) |
| `fold_state` | folded / disordered, from UniProt HELIX/STRAND/disorder annotation | UniProt features | used upstream in the categorical `impact` |

**Cross-species conservation:**

| feature | what it is | source | status |
|---------|-----------|--------|--------|
| phyloP (region / whole-protein mean) | evolutionary conservation of the changed region | UCSC hg38 100-way phyloP bigWig + EBI coordinates API | tested — a small structural add (+0.037 AUC on a 363-pair sample); not shipped pending scale-up |

**Not available:** **NMD (nonsense-mediated decay)** — DoChaP filters NMD-targeted
transcripts at DB-build time, so degraded isoforms never reach the benchmark; NMD is
neither a feature nor a confound here.

---

## Appendix B — Comparable published methods

We looked for published methods attempting each task and compare them below. A recurring
theme: for the **structural** axis the prior work *folds every isoform* (expensive) and
*catalogues* the result, whereas DOMAS *predicts* the fold-change label from cheap
canonical-only features; for the **functional** axis the closest analogues score exon /
isoform importance from protein + conservation + expression features, on labels that are
related but not identical to ours (so AUCs are not head-to-head).

### B.1 Structural axis — `fold_change_prob`

| method | approach | reported result | relation to ours |
|--------|----------|-----------------|-------------------|
| **Song et al., *Genome Biology* 2025** (our label source) | Fold all 11,161 human isoforms with **AlphaFold2**; compare each to its canonical with TM-score + 5 other metrics (secondary structure, surface charge, radius of gyration, PTM SASA, IDR). | Structural similarity largely tracks sequence identity; flags **328 high-identity / low-TM** outliers (53 with clear domain alterations). A *characterisation*, not a predictor. | This is exactly the TM ground truth we predict. We reproduce their label from **4 cheap canonical features (no isoform folding)** at **AUC 0.90 / R² 0.66**, and quantify the high-id/low-TM blind spot they identified. |
| **Yang et al., bioRxiv 2024.01.30.578053** | AF2 structures of AS isoforms; statistical characterisation of where AS falls. | AS regions are **enriched in loops / exposed / disordered residues**; qualitative case studies (Septin-9, Tau). No TM predictor, no isoform-vs-isoform structural score. | Supports our finding that exposed/peripheral changes drive low TM; provides no quantitative predictor to compare against. |
| **ASpdb** (Nucleic Acids Research, Database issue 2025) | Knowledgebase of >3,400 canonical + >7,200 alternative-isoform AF2 models with **comparative structural-alteration calls** and visualisation. | A resource, not a predictive model; provides a second, **non-TM** structural-alteration label. | Complementary independent label source (flagged in our EXPERIMENTS E40 as a future cross-check; its calls are per-entry web content, not bulk-downloadable). |
| **DIGGER / DIGGER 2.0** (Louadi et al., NAR 2021 / 2025) | Augment the protein–protein interaction network with **domain–domain interactions** + splicing; identify which interactions/domains an AS event disrupts (interaction rewiring). | Network/rule-based (no AUC); maps exon → lost domain → lost interactions. | Different structural readout — *interaction* consequence, not global fold. Closest in spirit to our `max_cov_loss` (domain loss) feature; complementary to a fold-change probability. |

**Bottom line (structural).** No prior method offers a cheap *predictor* of the
fold-change (TM) label — they fold the isoform or catalogue outcomes. DOMAS's contribution
is a calibrated ~AUC-0.90 predictor from canonical-only features, with **`pae_global`
(PAE)** as the dominant signal — to our knowledge not previously used for this task.

### B.2 Functional axis — `impact_prob`

| method | approach | reported result | relation to ours |
|--------|----------|-----------------|-------------------|
| **ExonImpact** (2017, [PMC5390777](https://pmc.ncbi.nlm.nih.gov/articles/PMC5390777)) | **Random forest** on protein-level features of the spliced exon: predicted disorder (SPINE-D), ASA + secondary structure (SPINE-X), Pfam-domain overlap, PTM sites (dbPTM), conservation, exon length. | **10-fold CV AUC 0.83–0.835.** Positives = 4,211 HGMD splice-disrupting exons; negatives = 2,664 1000-Genomes in-frame indel exons. | The closest analogue. Higher AUC, but on an **easier, different label** (disease-exon vs benign-indel-exon, with exon length as a feature) and using *predicted* disorder/ASA where we use AlphaFold-derived burial. Our label (pathogenic-vs-benign among variant-overlapping regions, length-controlled, gene-grouped) is harder — AUCs are not head-to-head. |
| **TRIFID** (Pozo et al., *NAR Genomics & Bioinformatics* 2021) | Gradient-boosted trees scoring the **functional relevance of each splice isoform** (0–1) from cross-species conservation, expression, and APPRIS/Pfam-integrity features. | Conservation is the top predictor; high-scoring isoforms' exons are under purifying selection, low-scoring ones neutral. Bimodal scores. Integrated into APPRIS. | Isoform-level (not region-level) and no single AUC against a variant label. Its top signal (conservation) is the axis we tested as phyloP (small add). Complementary; different granularity. |
| **pext** (Cummings et al., *Nature* 2020) | Not ML: **GTEx transcript-expression**-based per-base annotation (proportion expressed across transcripts) — is the region actually expressed? | Low-pext filtering removes **22.8%** of pLoF in haploinsufficient genes but only **3.8%** of ClinVar pathogenic variants. | Orthogonal axis (expression) to ours (constraint + structure). We *considered but did not use* GTEx (E41); pext is the strong published version of that idea and is complementary to `impact_prob`. |

**Not direct comparators.** SpliceAI / MMSplice / AbSplice predict *whether a variant
alters splicing* (the upstream event), not the functional consequence of the resulting
spliced region — a different task one level up. AlphaMissense scores single missense
variants; we *consume* it as `region_am` and operate at the region level.

**Bottom line (functional).** DOMAS's `impact_prob` sits in the same family as ExonImpact
and TRIFID (protein/constraint features → functional-importance score) but is
**region-level, calibrated, gene-grouped, and length-controlled**, and it is explicitly
paired with the orthogonal structural score. Its AUC (0.766) is lower than ExonImpact's
headline 0.83, but on a deliberately harder, confound-controlled label; the honest framing
throughout is that it is an **enrichment/prioritisation** signal, not a classifier.

**Sources.**
[Song et al., Genome Biology 2025](https://link.springer.com/article/10.1186/s13059-025-03744-x) ·
[Yang et al., bioRxiv 578053](https://www.biorxiv.org/content/10.1101/2024.01.30.578053v2) ·
[ASpdb, NAR 2025](https://academic.oup.com/nar/article/53/D1/D331/7893317) ·
[DIGGER 2.0, NAR 2025](https://academic.oup.com/nar/article/53/W1/W245/8126897) ·
[ExonImpact, PMC5390777](https://pmc.ncbi.nlm.nih.gov/articles/PMC5390777) ·
[TRIFID, NAR Genom. Bioinform. 2021](https://academic.oup.com/nargab/article/3/2/lqab044/6281449) ·
[pext / Cummings et al., Nature 2020](https://www.nature.com/articles/s41586-020-2329-2)
