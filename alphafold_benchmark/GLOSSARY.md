# Glossary

Terms used across the AlphaFold-isoform / DOMAS-validation investigation,
grouped by area. Definitions are scoped to how the term was used here.

---

## Molecular biology

- **Alternative splicing (AS)** — one gene producing multiple mRNAs (and proteins)
  by including/excluding different exons.
- **Isoform** — one of the alternative protein products of a gene. UniProt names
  them `ACC-1` (canonical), `ACC-2`, `ACC-3`, …
- **Canonical isoform / sequence** — UniProt's chosen reference isoform for a
  protein (the `-1`, or the bare accession).
- **Transcript** — an mRNA product of a gene; in DoChaP identified by an Ensembl
  ID (`ENST…`) with a version suffix.
- **Exon / exon skipping** — a segment of a transcript; skipping removes one.
- **CDS (coding sequence)** — the protein-coding portion of a transcript; DoChaP
  stores per-exon CDS genomic coordinates.
- **Splice junction** — the boundary between two joined exons; DOMAS maps these to
  domain structures.
- **Domain** — a structural/functional protein module (e.g. a kinase domain).
- **Motif** — a short functional sequence element (e.g. the BH3 motif in BAX).
- **Post-translational modification (PTM)** — chemical change to a residue after
  translation (phosphorylation, acetylation, etc.).
- **Functional residue / functional site** — a curated position with a role:
  active site (ACT_SITE), ligand binding (BINDING), modified residue (MOD_RES),
  general site (SITE), motif (MOTIF), disulfide bond (DISULFID), glycosylation
  (CARBOHYD).
- **Intrinsically disordered region (IDR)** — a segment with no stable fold;
  low AlphaFold confidence by nature.
- **N-terminus / C-terminus** — the start / end of a protein chain.
- **Ortholog / homolog** — corresponding genes across species / by common descent;
  the basis of cross-species conservation.
- **Conservation (evolutionary)** — how unchanged a position is across species;
  a proxy for functional importance.
- **Missense variant** — a single-amino-acid-changing mutation.
- **Pathogenic / benign variant** — disease-causing vs harmless; here from
  UniProt humsavar (LP/P vs LB/B).
- **Case-study genes** — CFHR3 (Sushi/CCP domain), MDH2 (malate dehydrogenase),
  BAX (apoptosis, BH3), CDK1 (kinase, Thr15), CXCR3 (7-transmembrane chemokine
  receptor), SRP9 (signal-recognition particle), EEF1AKMT3 (SAM methyltransferase),
  SEPTIN9 (GTPase), MAPT/Tau (3R/4R microtubule-binding repeats).

## Structural biology / AlphaFold

- **AlphaFold / AlphaFold2 (AF2)** — deep-learning protein-structure predictor.
- **AlphaFold DB (AFDB)** — repository of AF-predicted structures, one per UniProt
  entry (`AF-<acc>-F1-model_v6.pdb`).
- **pLDDT** — AF's per-residue confidence (0–100); high = well-modeled/ordered,
  low = disordered/uncertain.
- **PAE (predicted aligned error)** — AF's confidence in the *relative* position of
  residue pairs; captures inter-domain rigidity. **The single strongest fold-change
  feature (E38):** canonical-structure PAE lifts structural AUC 0.78 → 0.90. Two
  summaries used: `pae_global` (whole-structure mean, a protein-level fold-stability
  prior) and `pae_reg2rest` (mean PAE between the changed region and the rest,
  region-specific). Downloaded as `predicted_aligned_error_v6.json` from AFDB.
- **TM-score** — structural-similarity metric between two structures, 0–1; **>0.5 =
  same fold, <0.5 = different fold**. The ground-truth label in this work.
- **RMSD / GDT-TS** — other structure-comparison metrics (seen in a literature search).
- **SASA (solvent-accessible surface area)** — surface area a water-sized probe can
  reach on a residue; high = surface, ~0 = buried.
- **RSA (relative solvent accessibility)** — SASA normalized by the residue's max
  (0 = buried core, 1 = fully exposed); our "burial" feature.
- **Buried fraction** — share of a region's residues with RSA < 0.25.
- **Secondary structure** — local backbone conformation: α-helix (H), β-strand (E),
  coil/loop (C/-).
- **Contact / contact density** — residues whose Cβ atoms are within ~8 Å; used to
  test structural "embedding" of the changed region.
- **Radius of gyration / surface charge** — global shape/charge descriptors (columns
  in the Genome Biology Table S4).
- **ColabFold** — accessible AF2 (MSA-based) re-implementation; the recommended folder
  for the hard high-id/low-TM band that ESMFold over-preserves (E42).
- **ESMFold** — single-sequence (MSA-free) structure predictor built on ESM-2, run via
  HuggingFace `transformers` `EsmForProteinFolding` (E42). Its TM reproduces AF2-TM
  (Pearson 0.90) but is systematically optimistic on subtle fold changes.
- **tmtools** — Python wrapper for US-align/TM-align; used to compute TM-score between
  two folded structures (E42).
- **EBI Proteins coordinates API** — `https://www.ebi.ac.uk/proteins/api/coordinates/<acc>`;
  maps UniProt protein positions to genomic exon coordinates (protein→genome), used to
  place the changed region on the genome for phyloP lookup (E39).

## Statistics / machine learning

- **Pearson correlation** — linear association between two continuous variables.
- **Spearman correlation** — rank-based (monotonic) association; robust to nonlinearity.
- **Confound** — a variable that drives both predictor and outcome (here: sequence
  identity, and region length).
- **Partial correlation / "controlling for X"** — the association that remains after
  removing X's linear contribution from both variables. The key test in this work.
- **R² (coefficient of determination)** — fraction of variance in the target
  explained by a model.
- **Incremental R²** — the R² a feature adds on top of a baseline model.
- **RMSE** — root-mean-square error of a regression's predictions.
- **AUC (ROC-AUC)** — probability a random positive scores above a random negative;
  0.5 = chance, 1.0 = perfect. Computed here via the Mann–Whitney rank formula.
- **Sensitivity (recall) / specificity** — true-positive rate / true-negative rate.
- **Precision** — fraction of positive calls that are correct.
- **Base rate** — prevalence of the positive class; **lift** = precision ÷ base rate.
- **Majority-class baseline** — accuracy of always predicting the commonest class.
- **Confusion matrix** — TP / FP / TN / FN counts at a decision threshold.
- **Cross-validation (k-fold, out-of-fold)** — train on k−1 folds, predict the held-out
  fold; gives honest (non-overfit) performance. Used 5-fold throughout.
- **Logistic regression** — linear model for a binary outcome (fit here by
  numpy gradient descent).
- **Least squares / linear regression** — fit continuous outcome (via `np.linalg.lstsq`).
- **Enrichment / dose-response** — outcome rate rising across ordered predictor levels.
- **Youden's J** — sensitivity + specificity − 1; picks a balanced threshold.
- **Stratified sampling** — sampling to balance a variable (here: equal per identity band).
- **Standardization (z-score)** — rescaling features to mean 0, sd 1 before fitting.
- **Ground truth / label** — the target being predicted (TM verdict; variant overlap).
- **Positive-only labels** — labels where absence ≠ negative (e.g. "no known
  pathogenic variant" doesn't prove neutrality); limits specificity claims.
- **Contingency table (row% / col%)** — cross-tab of two categoricals; row% within a
  row, col% within a column.
- **Calibration** — whether a predicted probability matches the observed rate (a "0.6"
  really means 60%). Logistic is naturally calibrated; trees/NNs usually need
  post-hoc Platt/isotonic calibration.
- **Gene-grouped (grouped) CV** — cross-validation where all of a gene's isoforms sit
  in one fold, so no gene is in both train and test. The leakage-free evaluation when
  a feature is gene-level (LOEUF) — used for the shipped scores.
- **Suppressor / interaction effect** — a feature useless *alone* that becomes useful
  *in combination* (here AM: ~0 for TM by itself, but +0.03 AUC on top of
  identity+burial). Logistic captures it only via explicit interaction terms; trees
  capture it natively.
- **Random forest / gradient boosting / hist-gradient-boosting** — tree-ensemble
  classifiers (sklearn); tested and *not* better here (few, monotone features).
- **Neural network / MLP (multi-layer perceptron)** — feed-forward net; tested (32,16)
  and no better — wrong regime for 6 tabular features.
- **Bayesian inference** — gives posterior *uncertainty* on predictions/coefficients;
  same point accuracy as frequentist logistic, so not pursued for AUC.

## Tools / software

- **HMMER** — profile-HMM search suite. **hmmsearch** (profiles → sequences, fast for
  bulk) vs **hmmscan** (sequences → profiles). **hmmpress** compiles the HMM library.
- **`--cut_ga`** — HMMER gathering-threshold cutoffs (Pfam's curated per-family cutoffs).
- **domtbl (`--domtblout`)** — HMMER per-domain hit table (coords, bitscore, coverage).
- **Coverage %** — `round(100·(hmm_to−hmm_from+1)/hmm_len)`: fraction of a Pfam model matched.
- **Biopython** — Python structural-biology toolkit; used `PDBParser`, `SASA.ShrakeRupley`.
- **Shrake–Rupley** — rolling-probe SASA algorithm (Biopython implementation).
- **pydssp** — pure-Python DSSP (secondary-structure assignment from coordinates).
- **DSSP / mkdssp** — canonical secondary-structure + accessibility program (binary
  unavailable here; pydssp used instead).
- **freesasa** — fast C SASA library (failed to build here; Biopython used instead).
- **pandas / numpy** — dataframes / numerical arrays.
- **matplotlib** — plotting (hexbin, boxplots).
- **scikit-learn** — the ML library used for the model-family comparison
  (`GroupKFold`, `LogisticRegression`, tree ensembles, `MLPClassifier`, `roc_auc_score`).
- **PyMuPDF (`fitz`)** — read/rasterize PDFs (extract text, render benchmark PDFs to PNG).
- **pyBigWig** — remote/local bigWig reader (proposed for phyloP; not run).
- **curl / xargs -P** — parallel bulk download of AF structures.

## Databases / data resources

- **UniProt (SwissProt)** — curated protein knowledgebase; source of canonical
  sequences, features, isoforms, disease variants.
- **UniProt reference proteome** — per-species canonical set (`UP000005640_9606`
  = human).
- **UniProt varsplic** — all curated splice-variant isoform sequences
  (`uniprot_sprot_varsplic.fasta`).
- **UniProt humsavar** — disease-associated single-AA variants in canonical
  coordinates, classified LP/P (pathogenic), LB/B (benign), US (uncertain).
- **Ensembl** — genome/transcript resource; peptide FASTA + transcript IDs used by DoChaP.
- **AlphaFold DB** — see above.
- **AlphaMissense** — per-variant pathogenicity predictions (mean per residue used as
  a constraint proxy, `region_am`).
- **Pfam** — protein-family HMM library (CC0); the HMM domains DOMAS scans.
- **InterPro** — integrative domain database; supplies representative-domain `type`.
- **ClinVar** — clinical variant archive (humsavar aggregates its interpretations).
- **OMIM** — Mendelian-disease/gene catalog (referenced by humsavar disease names).
- **gnomAD** — population variation database; source of per-gene loss-of-function
  constraint. Provisioned into the `gene_constraint` table.
- **LOEUF** (`oe_lof_upper`) — gnomAD's observed/expected LoF upper-bound; **lower =
  more loss-intolerant** (dosage-sensitive gene). The non-circular gene-level feature.
- **pLI** — gnomAD probability a gene is LoF-intolerant (near 1 = intolerant); weaker
  than LOEUF here.
- **ASpdb** — isoform-structure knowledgebase (>7,200 alt-isoform AF2 models +
  structural-alteration calls); identified as an independent cross-check, not integrated.
- **APPRIS / TRIFID** — principal-isoform annotation; **SPADE** is its Pfam
  domain-integrity score.
- **GTEx** — tissue expression (proposed as an isoform-relevance signal; not used).
- **ConSurf / phyloP / phastCons** — per-residue / per-base conservation resources.
  **phyloP now run (E39)** via the UCSC hg38 100-way bigWig (remote `pyBigWig`): a
  modest structural add (+0.037 AUC on 363 pairs), partly beyond AlphaMissense.
- **pyBigWig** — reads bigWig tracks, including *remote* range queries over HTTP (used
  to pull phyloP for genomic intervals without downloading the ~9 GB file).
- **IsoformMapper / AlphaSync** — related isoform tools/DBs found in the literature search.
- **ECOD** — structural-domain classification (deliberately excluded from enrichment).

## DOMAS / DoChaP internals

- **DOMAS** — the tool being validated: maps splicing changes to domain structures
  and scores their impact.
- **DoChaP** — domain-charting database; `DB_merged.sqlite` holds Genes / Transcripts
  / Proteins / **RepresentativeDomains** (InterPro-keyed) built by
  `RepresentativeDomainsBuilder`.
- **enrichment.sqlite** — DOMAS's local evidence DB: `afdb_plddt`, **`afdb_rsa`**
  (per-residue relative solvent accessibility / burial), `am_pathogenicity`,
  `uniprot_feature`, `uniprot_alias`, `ensembl_sequence`, **`gene_constraint`**
  (gnomAD LOEUF/pLI per accession), `pfam_match`, `source_meta`.
- **impact level** — DOMAS's categorical change severity: `none` / `low` / `moderate` /
  `high` (plus `gain`, `n/a`), from `hmm_change_impact`.
- **`impact_prob`** — calibrated **functional** probability (pathogenicity-relevance),
  logistic over region_am + LOEUF + max_cov_loss + buried_frac
  (`utils.impact_probability`, gene-grouped AUC 0.765).
- **`fold_change_prob`** — calibrated **structural** probability, P(TM-score < 0.5),
  logistic over identity + burial + region_am + LOEUF + max_cov_loss
  (`utils.fold_change_probability`, gene-grouped AUC 0.777). Burial enters with the
  **opposite sign** to `impact_prob`.
- **`buried_frac` / `region_rsa`** — fraction of the changed region with RSA<0.25 /
  mean RSA, from `afdb_rsa` (AlphaFold burial).
- **`changed_region`** — the divergent span between canonical and isoform; now returns
  canonical + alt spans and a **kind** (`substitution` / `insertion` / `deletion`).
- **`insertion_impact`** — scores alt-only added residues (softened to escalate on
  domain placement, not size).
- **`fold_state`** — 'folded' / 'disordered' / None from UniProt HELIX/STRAND/disorder.
- **`region_am`** — mean AlphaMissense over the changed region.
- **`func_site`** — whether the changed region overlaps a UniProt functional residue.
- **RepresentativeDomains `type`** — InterPro entry type: Domain / Repeat / Family /
  Homologous_superfamily / Conserved_site / Active_site / Binding_site / PTM.
- **NON_COMPARISON_EVENTS** — result rows with no domain comparison (gene not in DB,
  no canonical junctions, etc.).

## File formats

- **FASTA** — plain sequence format (`>header` then residues).
- **PDB / mmCIF** — 3D atomic-coordinate formats (AF structures).
- **OOXML `.xlsx`** — Excel; the Genome Biology supplementary used a non-standard
  namespace (`purl.oclc.org/ooxml/...`) requiring a custom parser.
- **domtbl** — HMMER domain-hit table.
- **bigWig** — indexed genome-track format (phyloP/phastCons).
- **`.gz`** — gzip compression.
