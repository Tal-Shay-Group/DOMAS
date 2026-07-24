# Data sources and exactly what was retrieved

Every external source touched in this investigation, the precise file/query, what
it gave us, and its license/access notes. Large downloads live outside the repo
(scratchpad or `DoChaP-db/data/`); the derived CSVs in this folder are committed.

---

## 1. Genome Biology 2025 — the ground-truth dataset

**Paper:** Miller et al., *Predicting the structural impact of human alternative
splicing*, Genome Biology 2025, doi:10.1186/s13059-025-03744-x.
- PMC: <https://pmc.ncbi.nlm.nih.gov/articles/PMC12442299/>
- Springer: <https://link.springer.com/article/10.1186/s13059-025-03744-x>
- Local PDF: `DOMAS/papers/Predicting the structural impact of human alternative splicing.pdf`

**Supplementary retrieved** (static-content Springer pattern
`https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-025-03744-x/MediaObjects/13059_2025_3744_MOESM{N}_ESM.xlsx`):
- **MOESM3 = Table S2** (`moesm3.xlsx`): **328** high-identity/low-TM outlier
  isoforms — columns: isoform ID, identity, TM-score, gene, protein name.
  → extracted to `table_s2_changes.csv`.
- **MOESM5 = Table S4** (`moesm5.xlsx`): **11,161** isoforms × 28 columns — isoform
  ID, length, pLDDT, class, sequence identity, **TM-score**, helix/sheet/loop
  fractions, surface charge, radius of gyration, PTM-change and IDR-percentage
  fields. Originally only 7 of the 28 columns were kept (`table_s4_all.csv`); the
  full 28 were later recovered to **`table_s4_full.csv`** (E35). The 21 extra columns
  turned out either **circular** for TM (isoform SS/Rg/IDR come from the isoform's own
  fold) or non-additive (PTM-change).
- *Note:* the xlsx used a non-standard OOXML namespace (`purl.oclc.org/ooxml/...`);
  parsed by reading `sheet1.xml` + `sharedStrings.xml` directly.

This is the **TM-score ground truth** for the whole benchmark: 11,161 isoforms,
5,923 genes.

## 2. UniProt — sequences and disease variants

Base FTP: `https://ftp.uniprot.org/pub/databases/uniprot/current_release/`

- **Per-isoform FASTA (REST)** — `https://rest.uniprot.org/uniprotkb/<isoform>.fasta`
  for the 7 pilot isoforms (Q02985-2, P40926-2, Q07812-4, P06493-2, P49682-2,
  P49458-2, Q96AZ1-2) + their canonicals. Used for the exact DoChaP mapping.
- **varsplic** — `.../knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz`
  (8.6 MB): **41,343** curated splice-variant isoform sequences; covered **99%** of
  Table S4 isoforms. → alt-isoform sequences for the 11k benchmark.
- **Human reference proteome** —
  `.../reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz` (7.8 MB):
  **20,652** canonical sequences; covered **100%** of Table S4 accessions. → canonical
  sequences for the 11k benchmark.
- **humsavar** — `.../variants/humsavar.txt` (9.1 MB): **84,845** disease-associated
  single-AA variants in canonical coordinates — **33,233 LP/P** (pathogenic),
  **39,962 LB/B** (benign), 11,650 US. → the functional-axis labels (E16–E18).
- License: UniProt CC BY 4.0.

## 3. AlphaFold DB — 3D structures

- **API** (to resolve the file URL): `https://alphafold.ebi.ac.uk/api/prediction/<acc>`
  → revealed current version **v6**.
- **Structures** — `https://alphafold.ebi.ac.uk/files/AF-<acc>-F1-model_v6.pdb`,
  downloaded for **5,880** canonical accessions (the benchmark proteins),
  parallel via `curl`/`xargs`. → per-residue RSA + secondary structure + contacts
  (E12–E13). Stored in scratchpad (`afstruct/`), not committed.
- **PAE** — `https://alphafold.ebi.ac.uk/files/AF-<acc>-F1-predicted_aligned_error_v6.json`,
  downloaded for **5,557** canonical accessions (parallel `curl`/`xargs`, scratchpad
  `pae/`, not committed). → `pae_global` + `pae_reg2rest` (E38, the strongest fold-change
  feature) in `pae_features_full.csv`. NxN predicted-aligned-error matrix per protein.
- License: AlphaFold DB CC-BY 4.0.

## 3b. UCSC phyloP + EBI coordinates — cross-species conservation (E39)

- **EBI Proteins coordinates API** — `https://www.ebi.ac.uk/proteins/api/coordinates/<acc>`
  (JSON via `Accept: application/json`): per-exon protein→genome mapping (chromosome,
  strand, exon protein/genome spans). Maps the changed region onto hg38.
- **UCSC hg38 100-way phyloP** — `http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw`
  (~9 GB), queried **remotely** with `pyBigWig` range requests (no full download). → mean
  phyloP over the changed region / whole protein (`phylop_features.csv`). License: UCSC free.

## 3c. ESMFold + tmtools — self-folding the isoform (E42)

- **ESMFold** via HuggingFace `transformers` `EsmForProteinFolding` (`facebook/esmfold_v1`,
  ~2.7 GB weights + wrapped ESM-2 LM, HF cache). CPU folding ~20 s per short pair.
- **tmtools** (pip; wraps US-align) for TM-score between folded structures. → `esmfold_results.csv`.

## 3d. ASpdb — attempted independent structural label (E40, blocked)

- **Bulk file** — `https://biodataai.uth.edu/ASpdb/tables/AS_event_info.txt` (14,645 AS
  events: acc, gene, isoform pair, lengths, event type, changed sequences). Contains **no
  usable structural-alteration call** — the comparative structural verdicts are per-entry
  web content only. Not usable as a bulk label. (Website: `https://biodataai.uth.edu/ASpdb/`.)

## 4. DoChaP — transcript / domain database (local)

- **`DoChaP-web/DB_merged.sqlite`** — Genes, Transcripts, Proteins
  (`protein_interpro_id` = UniProt accession), **RepresentativeDomains** (InterPro-
  keyed, with `type`). Used for isoform→transcript mapping and gene PDFs. 97% of
  Table S4 accessions are present.
- **`DoChaP-db/data/enrichment.sqlite`** — DOMAS's local evidence DB:
  `afdb_plddt` (per-residue pLDDT), `am_pathogenicity` (AlphaMissense per-residue mean),
  `uniprot_feature` (functional sites / fold-state), `uniprot_alias`,
  `ensembl_sequence` (Ensembl peptides for sequence matching), `source_meta`.

## 5. Pfam — domain HMMs (local, already provisioned)

- **`DoChaP-db/data/pfam/Pfam-A.hmm`** (hmmpress-ed): the profile-HMM library
  scanned by `hmmsearch --cut_ga` for the impact score and SPADE. 19k+ families.
- License: Pfam CC0.

## 6. gnomAD — per-gene loss-of-function constraint

- **File** — `gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz` (~5 MB) from
  `https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/`.
  Columns used: `gene`, **`oe_lof_upper` (LOEUF)**, `pLI`. Gene symbol → UniProt
  accession via the reference-proteome headers (17,142 accessions mapped) →
  `gene_constraint` table. The non-circular gene-level feature for both calibrated
  scores. License: gnomAD terms (freely available, CC0).

## 6b. AlphaFold burial (derived from §3 structures)

- Per-residue **RSA** (Shrake–Rupley via Biopython) computed from the 5,880 AF
  structures → `afdb_rsa` table. No new download — same structures as §3.

## 7. Literature search — resources identified (not all used)

Via web search; catalogued for context / future work:
- **ASpdb** — isoform-structure knowledgebase, >7,200 alt-isoform AF2 models with
  comparative structural-alteration calls. NAR 2025,
  <https://academic.oup.com/nar/article/53/D1/D331/7893317>. *Flagged as an
  independent functional cross-check; not integrated.*
- **bioRxiv 2024.01.30.578053** (Yang et al.) — AF2 isoform study; qualitative
  case studies only, no TM/accession ground truth. Local PDF in `DOMAS/papers/`.
  *Excluded.*
- **IsoformMapper** — bioRxiv 2025.03.05.641708 (splice-variant community analysis).
- **AlphaSync** — bioRxiv 2025.03.12.642845 (AFDB synced to UniProt).
- **AF2/ESMFold/OmegaFold benchmark** — bioRxiv 2025.06.20.660709 (not isoform-specific).

## 8. Considered but NOT retrieved

- **AlphaMissense** raw file — already built into `am_pathogenicity`; not re-downloaded.
- **phyloP / phastCons** (UCSC bigWig) — proposed for cross-species conservation; the
  plan (DoChaP CDS→genome mapping + `pyBigWig` interval queries) was scoped, not run.
- **ClinVar** VCF — not needed; UniProt humsavar already gives variants in canonical
  coordinates, avoiding genomic mapping.
- **PAE** files from AFDB — deemed likely-redundant with burial and very large; skipped.

---

## Derived data committed in this folder

| file | what it is |
|------|-----------|
| `table_s2_changes.csv` / `table_s4_all.csv` | extracted Genome Biology supplementary |
| `full_pairs.csv` | 11,068 runnable pairs with canonical changed-region spans |
| `bench_full_results.csv` | per-pair DOMAS impact + kind vs TM |
| `bench_rich_results.csv` | per-pair impact + drivers (coverage loss, region_am, func_site, …) |
| `full_features.csv` / `sample_features.csv` | per-pair burial (RSA), buried_frac, DSSP secondary-structure |
| `contact_features.csv` | per-pair contact-modularity features |
| `predict_eval_results.txt` | cross-validated AUC / R² comparison |
| `clinvar_enrichment_results.txt` | pathogenic vs benign enrichment |
| `calib_model.json` / `foldchange_model.json` | shipped logistic coefficients for `impact_prob` / `fold_change_prob` |
| `fit_calibrated.py` / `fit_foldchange.py` | fit + gene-grouped CV for the two scores |
| `improve_proto.py` | peak-AM / gnomAD feature prototype (peak rejected) |
| `ml_compare.py` | model-family comparison (logistic chosen over trees/NN) |
| `variant_impact_contingency.py` / `variant_impact_by_score_decile.py` | pathogenic/benign × impact tables |
| `domas_vs_tm.png` / `pred_vs_tm.png` | the two scatter figures |
| `esm650_compare.py` / `esm650_accuracy.py` | ESM-2 650M vs 150M (E33) |
| `exp_common.py` | Phase-10 shared harness (master table + gene-grouped CV, all 3 metrics) |
| `exp1_position.py` | position of changed region (E34) |
| `exp2_tables4.py` | recovered Table S4 columns (E35) |
| `exp3_multitask.py` | multi-task joint model (E36) |
| `exp4_esmhead.py` | ESM multi-pool vs mean-pool (E37) |
| `exp5_pae.py` / `exp5b_pae_rigor.py` / `exp5_pae_full.py` | **PAE — the ceiling-break (E38)** |
| `exp6_conservation.py` | phyloP conservation pipeline (E39) |
| `esmfold_probe.py` / `esmfold_run.py` | ESMFold feasibility + self-folding test (E42) |
