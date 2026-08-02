# Review cases — one example per decision the algorithm makes

> **Generated file — do not edit by hand.** Rebuild with
> `python3 tests/generate_review_index.py`.
>
> The picks are properties of the output, so they stop being true when the
> outputs move — which has happened both when a fixture was refreshed and when
> DoChaP was rebuilt. Deriving them is what keeps this file honest.

- reference outputs as of commit `f2adf69`
- DoChaP database built `2026-07-31 09:29`
- generated `2026-08-02 08:28`

The PDFs are not in git (`.gitignore` excludes `*.pdf`). Rebuild them with
`python3 tests/generate_reference_outputs.py` and
`python3 tests/generate_review_pdfs.py`, then collect them into one directory
with `python3 tests/collect_review_pdfs.py`.

---

## A. Choosing the comparable transcript

| # | What it shows | Where | Detail |
| :-- | :--- | :--- | :--- |
| A1 | the outside-exon gate picks a DIFFERENT transcript than longest-CDS - the one case where is_most_like_canonical changes the answer | `*_TMED8_*.pdf` in `review_cases/leafcutter_all_transcripts` | TMED8 `chr14:clu_9428` · 5 tx, 2 comparable · longest CDS `XM_017021224.2`, most like canonical `NM_001346133.2` · dropped domain |
| A2 | is_most_like_canonical is UNSET - no candidate passes the gate, so longest-CDS decides among coding candidates alone | `*_GLIPR1_*.pdf` in `review_cases/leafcutter_all_transcripts` | GLIPR1 `chr12:clu_6316` · 5 tx, 3 comparable · longest CDS `ENST00000456650.7` · dropped domain, shorter |
| A3 | a non-coding candidate holds the longest CDS - step 1 excludes it from selection, so the tag goes to a coding transcript instead | `*_FAM209A_*.pdf` in `review_cases/leafcutter_all_transcripts` | FAM209A `chr20:clu_4032:FAM209A` · 4 tx, 3 comparable · 1 non-coding · longest CDS `XM_047439965.1` · dropped domain |
| A4 | several comparable transcripts, so the flags have to choose between them | `*_PRRG3_*.pdf` in `tests/reference_outputs/category_examples__restrict_False__representative_True` | PRRG3 `AF_ENSG00000130032_X_265953` · 11 tx, 2 comparable · longest CDS `ENST00000884468.1`, most like canonical `ENST00000370353.3` · merged domain, same_domains |

## B. Domain resolution (InterPro entry-type ladder)

| # | What it shows | Where | Detail |
| :-- | :--- | :--- | :--- |
| B1 | the tier ladder discards a redundant annotation (2 -> 1 annotations) | `*_CACNG3_*.pdf` in `tests/reference_outputs/ioe_csv__restrict_True__representative_True` | CACNG3 `ENST00000005284.4` · 254 such transcripts |
| B2 | demote, don't delete - a member-DB hit kept as the only evidence (2 -> 1 annotations) | `*_CACNG3_*.pdf` in `tests/reference_outputs/ioe_csv__restrict_True__representative_True` | CACNG3 `ENST00000005284.4` · 384 such transcripts |
| B3 | all three tiers present in one transcript (6 -> 4 annotations) | `*_DIAPH2_*.pdf` in `tests/reference_outputs/ioe_csv__restrict_True__representative_True` | DIAPH2 `ENST00000324765.13` · 113 such transcripts |
| B4 | Site/PTM entries removed - unconditional, unlike the tier rule (11 -> 8 annotations) | `*_KMT2D_*.pdf` in `tests/reference_outputs/category_examples__restrict_True__representative_True` | KMT2D `ENST00000301067.12` · 5 such transcripts |
| B5 | the same accession more than once (3 -> 3 annotations) | `*_BPIFB2_*.pdf` in `tests/reference_outputs/ioe_csv__restrict_True__representative_True` | BPIFB2 `ENST00000170150.4` · 366 such transcripts |
| B6a | an instance of a repeat DISCARDED for overlapping its neighbour of the same accession - the domain count drops though both instances are real (11 -> 9 annotations) | `*_TOPBP1_*.pdf` in `review_cases/leafcutter` | TOPBP1 `ENST00000260810.10` · IPR036420 102-199 dropped, 199-302 kept, overlap 1 aa (1% of the shorter) · 43 such transcripts |
| B6b | an instance of a repeat DISCARDED for overlapping its neighbour of the same accession - the domain count drops though both instances are real (10 -> 8 annotations) | `*_ZNF317_*.pdf` in `review_cases/leafcutter` | ZNF317 `ENST00000247956.11` · IPR036236 175-231 dropped, 230-287 kept, overlap 2 aa (4% of the shorter) · 43 such transcripts |
| B6c | an instance of a repeat DISCARDED for overlapping its neighbour of the same accession - the domain count drops though both instances are real (9 -> 5 annotations) | `*_ZNF680_*.pdf` in `tests/reference_outputs/category_examples__restrict_True__representative_True` | ZNF680 `ENST00000309683.11` · IPR036236 246-303 dropped, 190-247 kept, overlap 2 aa (3% of the shorter) · 43 such transcripts |
| B6d | an instance of a repeat DISCARDED for overlapping its neighbour of the same accession - the domain count drops though both instances are real (5 -> 4 annotations) | `*_CSMD2_*.pdf` in `tests/reference_outputs/category_examples__restrict_True__representative_True` | CSMD2 `ENST00000465819.1` · IPR035976 95-158 dropped, 154-222 kept, overlap 5 aa (8% of the shorter) · 43 such transcripts |
| B6e | an instance of a repeat DISCARDED for overlapping its neighbour of the same accession - the domain count drops though both instances are real (20 -> 18 annotations) | `*_SORL1_*.pdf` in `tests/reference_outputs/category_examples__restrict_True__representative_True` | SORL1 `ENST00000260197.12` · IPR003961 2025-2112 dropped, 1934-2029 kept, overlap 5 aa (6% of the shorter) · 43 such transcripts |
| B7a | the same rule applied where it IS right - two calls of one region, near-identical spans, one discarded (2 -> 1 annotations) | `*_SLC6A19_*.pdf` in `tests/reference_outputs/ioe_csv__restrict_True__representative_True` | SLC6A19 `ENST00000304460.11` · IPR000175 31-610 dropped, 25-608 kept, overlap 578 aa (100% of the shorter) · 20 such transcripts |
| B7b | the same rule applied where it IS right - two calls of one region, near-identical spans, one discarded (2 -> 1 annotations) | `*_WNT16_*.pdf` in `tests/reference_outputs/ioe_csv__restrict_True__representative_True` | WNT16 `ENST00000222462.3` · IPR005817 51-365 dropped, 49-365 kept, overlap 315 aa (100% of the shorter) · 20 such transcripts |
| B7c | the same rule applied where it IS right - two calls of one region, near-identical spans, one discarded (2 -> 1 annotations) | `*_DNASE1L1_*.pdf` in `review_cases/leafcutter` | DNASE1L1 `ENST00000862494.1` · IPR016202 3-294 dropped, 2-293 kept, overlap 291 aa (100% of the shorter) · 20 such transcripts |

## C. Intron retention

| # | What it shows | Where | Detail |
| :-- | :--- | :--- | :--- |
| C1 | a retained intron producing real domain calls | `*_ZNF593_*.pdf` in `tests/reference_outputs/ioe_csv__restrict_True__representative_True` | ZNF593 `RI_ENSG00000142684_1_269` · 2 tx, 1 comparable · longest CDS `ENST00000270812.6` · added_domain, shorter |
| C2 | a retained intron in a cluster with a non-coding transcript | `*_MSLNL_*.pdf` in `tests/reference_outputs/ioe_csv__restrict_True__representative_True` | MSLNL `RI_ENSG00000162006_16_14761` · 2 tx, 1 comparable · 1 non-coding · longest CDS `ENST00000537221.1` · dropped domain |

## D. Input format × AS event type

| # | What it shows | Where | Detail |
| :-- | :--- | :--- | :--- |
| D1a | plus-strand A3SS | `*_A3SS_COL18A1_*.pdf` in `review_cases/rmats` | COL18A1 `A3SS_ENSG00000182871_21_63` · 33 tx, 4 comparable · added_domain, dropped domain |
| D1b | minus-strand A3SS - the boundary that varies with strand | `*_A3SS_SBF1_*.pdf` in `review_cases/rmats` | SBF1 `A3SS_ENSG00000100241_22_78` · 80 tx, 10 comparable · added_domain, dropped domain, longer |
| D2a | plus-strand A5SS | `*_A5SS_GGA1_*.pdf` in `review_cases/rmats` | GGA1 `A5SS_ENSG00000100083_22_37` · 49 tx, 7 comparable · dropped domain, longer, split domain |
| D2b | minus-strand A5SS - the boundary that varies with strand | `*_A5SS_POFUT2_*.pdf` in `review_cases/rmats` | POFUT2 `A5SS_ENSG00000186866_21_62` · 33 tx, 1 comparable · added_domain, dropped domain |
| D3a | plus-strand SE - NO comparable transcript, canonical only | `*_SE_DIP2A_*.pdf` in `review_cases/rmats` | DIP2A `SE_ENSG00000160305_21_1` · 38 tx, 0 comparable |
| D3b | minus-strand SE | `*_SE_POFUT2_*.pdf` in `review_cases/rmats` | POFUT2 `SE_ENSG00000186866_21_29` · 33 tx, 5 comparable · dropped domain |
| D4a | plus-strand MXE | `*_MXE_PDXK_*.pdf` in `review_cases/rmats` | PDXK `MXE_ENSG00000160209_21_94` · 36 tx, 2 comparable · added_domain, dropped domain, shorter |
| D4b | minus-strand MXE | `*_MXE_ATP5PO_*.pdf` in `review_cases/rmats` | ATP5PO `MXE_ENSG00000241837_21_114` · 26 tx, 4 comparable · dropped domain, shorter |
| D5 | MAJIQ per-LSV classification | `*.pdf` in `review_cases/majiq` | — |
| D6a | one LeafCutter event, first of the two genes it names | `*_HYPK_*.pdf` in `review_cases/leafcutter` | HYPK `chr15:clu_17411:HYPK` · 3 comparable |
| D6b | the same event, second of the two genes it names | `*_SERF2_*.pdf` in `review_cases/leafcutter` | SERF2 `chr15:clu_17411:SERF2` · 4 comparable |

## E. Classification outcomes

| # | What it shows | Where | Detail |
| :-- | :--- | :--- | :--- |
| E1 | longer_domains - equal counts >1 on both sides, compared side longer | `*_DHX29_*.pdf` in `review_cases/leafcutter` | — |

---

## Table S4 outcome coverage

| Outcome | Clusters in category_examples |
| :--- | ---: |
| dropped domain | 12 |
| added_domain | 10 |
| longer | 7 |
| same | 7 |
| shorter | 3 |
| same_domains | 4 |
| merged domain | 5 |
| split domain | 2 |
| increased_domain_number | 2 |
| reduced_domain_number | 2 |
| shorter_domains | 1 |
| longer_domains | **0 — see E1** |

---

## Gaps

- None: every case above resolved against the current outputs.
- `no_gene_specified` cannot have a PDF: the figure is per-gene, and a cluster
  naming no gene has no transcripts, exons or domains to draw.
