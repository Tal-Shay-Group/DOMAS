# Review cases — one example per decision the algorithm makes

Picks re-derived against the reference outputs as of commit `85bf0b6` (the refreshed
SUPPA fixture). Paths are relative to the repository root.

`REF = tests/reference_outputs/<case>` — use the `restrict_False__representative_True`
variant of each case unless a row says otherwise; the four flag variants differ only
in PDF scope and representative-domain source.

> The ioe fixture was refreshed in `85bf0b6` and now samples 2-transcript genes, so
> multi-transcript examples come from `category_examples`. Every pick below was
> re-derived after that change; earlier drafts of this file named clusters that no
> longer exist.

---

## A. Choosing the comparable transcript

Ladder: **protein-coding → most-like-canonical → longest CDS**. The first flag is
left unset when nothing passes the outside-exon gate.

| # | What it shows | Gene | Cluster | Case |
| :-- | :--- | :--- | :--- | :--- |
| A1 | `is_most_like_canonical` **set** (14 of 59 clusters) | RIMS4 | `A3_ENSG00000101098_20_67494` | ioe_csv |
| A2 | `is_most_like_canonical` **unset** — no candidate passes the gate, only the longest-CDS flag is awarded (45 of 59) | IBSP | `A3_ENSG00000029559_4_18475` | ioe_csv |
| A3 | Cluster containing a **non-coding** transcript — excluded from selection by step 1, still compared and reported | FGF4 | `SE_ENSG00000075388_11_83184` | ioe_csv |
| A4 | **Several comparable transcripts**, 6 transcripts of which 2 non-coding — the most readable multi-candidate case | POLQ | `A3_ENSG00000051341_3_16572` | category_examples |
| A5 | Two comparable transcripts out of 10 | PRRG3 | `AF_ENSG00000130032_X_265953` | category_examples |
| A6 | Five transcripts, one non-coding, outcome `split domain` | PEX3 | `SE_ENSG00000034693_6_53833` | category_examples |

The transcript header carries the tags — `Protein: ENSP…  [longest CDS, most like
canonical]`. `Protein: N/A` marks a non-coding transcript.

---

## B. Domain resolution (InterPro entry-type ladder)

**Domain/Repeat** always kept → **Family/Homologous_superfamily** dropped only where a
higher tier covers >50 % of it → **member-database hits** dropped on the same rule.
Site/PTM removed outright.

| # | What it shows | Gene / transcript | Count |
| :-- | :--- | :--- | :--- |
| B1 | Reduction happens, 2 annotations → 1 | CACNG3 `ENST00000005284.4` | 28 transcripts |
| B2 | **Demote, don't delete** — a member-DB hit is kept because no InterPro Domain covers the region | CACNG3 `ENST00000005284.4` | 25 |
| B3 | **All three tiers** in one transcript, 6 annotations → 4 | DIAPH2 `ENST00000324765.13` | 2 |
| B4 | Site/PTM entries removed | — | **0 — see gaps** |
| B5 | **Repeated accession**: same domain id more than once | BPIFB2 `ENST00000170150.4` | 12 |

---

## C. Intron retention (new in this work)

The retained interval is drawn as a **dotted span inside the exon**, not as a raised
bracket — it is the opposite of an excision. The first-page junction table carries a
`feature_type` column, present only for events that contain a retained intron.

| # | What it shows | Gene | Cluster | Case |
| :-- | :--- | :--- | :--- | :--- |
| C1 | RI producing real domain calls (`added_domain`, `shorter`) | ZNF593 | `RI_ENSG00000142684_1_269` | ioe_csv |
| C2 | RI with a non-coding transcript, outcome `dropped domain` | MSLNL | `RI_ENSG00000162006_16_14761` | ioe_csv |
| C3 | RI from rMATS, both strands | PFKL (+), MCM3AP (−) | — | `review_cases/rmats/` |

Before `85bf0b6` these five ioe RI clusters produced **zero** analysed rows — the
fixture had frozen the single-junction representation, so the golden test asserted
the broken behaviour. They now yield 7 domain calls.

---

## D. Input format × AS event type

`ioe_csv` covers all seven SUPPA types. The other three formats had **no PDFs at all**
— their tests run `create_pdf=False` — so a set was generated here.

| Format | Coverage | Location |
| :--- | :--- | :--- |
| SUPPA (`ioe`) | A3, A5, AF, AL, MX, RI, SE — 5 clusters each | `REF/ioe_csv__…/` |
| category_examples | A3, A5, AF, AL, MX, SE | `REF/category_examples__…/` |
| hadas | human/mouse comparison | `REF/hadas_xlsx__…/` |
| rMATS | SE, A5SS, A3SS, MXE, RI — **one + and one − strand each** | `review_cases/rmats/` (10) |
| MAJIQ | per-LSV classification | `review_cases/majiq/` (4) |
| LeafCutter | incl. the multi-gene pair | `review_cases/leafcutter/` (6) |

The rMATS ten, paired by strand so the strand-dependent construction can be compared:

| Event | + strand | − strand |
| :--- | :--- | :--- |
| SE | `human_SE_DIP2A_10_*.pdf` | `human_SE_MCM3AP_9_*.pdf` |
| A5SS | `human_A5SS_PDXK_3_*.pdf` | `human_A5SS_POFUT2_4_*.pdf` |
| A3SS | `human_A3SS_COL18A1_2_*.pdf` | `human_A3SS_MCM3AP_1_*.pdf` |
| MXE | `human_MXE_PDXK_5_*.pdf` | `human_MXE_MCM3AP_6_*.pdf` |
| RI | `human_RI_PFKL_7_*.pdf` | `human_RI_MCM3AP_8_*.pdf` |

**Highest value to look at first:**

1. **A minus-strand A3SS or A5SS** — `human_A3SS_MCM3AP_1_*.pdf`,
   `human_A5SS_POFUT2_4_*.pdf`. Every such event was unanalysable before this work,
   so no one has ever seen one.
2. **A retained intron** — C1 or `review_cases/rmats/human_RI_PFKL_7_*.pdf`. New
   capability and a new drawing convention.
3. **The multi-gene LeafCutter pair** — `human_RBM5_5_*.pdf` and `human_RBM6_6_*.pdf`:
   one event, `chr3:clu_10461`, now analysed once per gene. Adjacent paralogues at
   3p21.3; previously only the first-named gene was reported.

---

## E. Classification outcomes (Table S4)

`category_examples` is built for this — 11 of the 12 outcomes have 1–13 clusters each.

| Outcome | Clusters | Example |
| :--- | ---: | :--- |
| dropped domain | 13 | PRRG3 |
| added_domain | 10 | SMAD9 |
| longer | 7 | ZNF680 |
| same | 7 | KMT2D |
| shorter | 3 | FKBP5 |
| same_domains | 3 | FKBP5 |
| merged domain | 3 | PRRG3 |
| split domain | 2 | PEX3 |
| increased_domain_number | 2 | ZNF141 |
| reduced_domain_number | 2 | DAZ4 |
| shorter_domains | 1 | PCDH7 |
| **longer_domains** | **0** | **missing** |

---

## Gaps

1. **`longer_domains`** — the only S4 outcome with no example anywhere.
2. **Site/PTM removal (B4) is never exercised.** No transcript in the ioe fixture
   carries an `Active_site` / `Binding_site` / `Conserved_site` / `PTM` InterPro entry,
   so that branch of the ladder is not merely unreviewable — it is **untested**.
3. **`no_gene_specified` has no PDF** — the LeafCutter runs that produce it are the
   ones without PDF generation.
4. The generated sets under `review_cases/` are **not** golden references and are
   compared against nothing; they exist to be looked at.
