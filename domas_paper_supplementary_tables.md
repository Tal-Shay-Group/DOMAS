# DOMAS — Supplementary Tables

Companion to the Materials and Methods section. Table numbers match the citations
in the manuscript.

---

## Table S1. Genome assemblies and annotation releases

The assembly and annotation behind the DoChaP database each species was built from.
Taken from the GFF3 headers of the Ensembl annotation files under `data/<species>/ensembl/`.

| Species | Common name | Assembly | Accession | Annotation date |
| :--- | :--- | :--- | :--- | :--- |
| *Homo sapiens* | human | GRCh38 | GCA_000001405.29 | 2025-05 |
| *Mus musculus* | mouse | GRCm39 | GCA_000001635.9 | 2025-05 |
| *Rattus norvegicus* | rat | GRCr8 | GCA_036323735.1 | 2024-03 |
| *Danio rerio* | zebrafish | GRCz11 | GCA_000002035.4 | 2018-04 |
| *Xenopus tropicalis* | frog | UCB_Xtro_10.0 | GCA_000004195.4 | 2021-02 |

> **To confirm before submission:** these are the builds present in the database
> directory. They should be checked against the builds actually used for the
> reported results, and the Ensembl and RefSeq release numbers added.

---

## Table S2. Per-tool mapping onto DOMAS event features

Every reader produces the same internal record: an event feature — an intron given
by two genomic coordinates, to be found either **spliced out** (`junction`) or
**retained** (`retained_intron`) — carrying its gene, chromosome, cluster and species.

### S2a. Where each field comes from

| Field | LeafCutter | rMATS-turbo | MAJIQ | SUPPA | hadas |
| :--- | :--- | :--- | :--- | :--- | :--- |
| gene | `genes` (significance), resolved via symbol | `GeneID` | `Gene ID` | gene prefix of `event_id` | `ensembl_h` / `genes` |
| gene symbol | `genes` | `geneSymbol` | `#Gene Name` | — (looked up) | `symbol_h` / `genes` |
| chromosome | `intron` field 1 | `chr` (`chr` prefix stripped) | `chr` (stripped) | `seqname` | from junction string |
| coordinates | `intron` fields 2–3 | per event type, see S2b | `Junctions coords`, `IR coords` | coordinate tokens of `event_id` | junction string |
| cluster | `intron` field 4, strand suffix dropped | `{type}_{gene}_{chr}_{n}` | `LSV ID` | `{type}_{gene}_{chr}_{n}` | `cluster` |
| species | `-specie` | `-specie` | `-specie` | `-specie` | per row (human + mouse) |

Notes:

- **LeafCutter** cluster ids drop the strand suffix (`chr:clu_<n>_<strand>` →
  `chr:clu_<n>`) so the significance and effect-size files join. Where a cluster
  names several genes it is analysed once per gene and the symbol is appended to
  the cluster id (`chr1:clu_15645:CDK11B`). A cluster naming no gene is kept and
  labelled `no_gene_specified`.
- **rMATS** counts clusters across all five files, not per file, so `n` is unique
  within a run.
- The species is stated by the user; only hadas carries it per row, being a
  human/mouse comparison.

### S2b. rMATS event types → features

Coordinates are genomic. `A5SS`/`A3SS` name the long and short forms of one exon
plus its flanking exon, so which boundary distinguishes the two forms depends on
the strand; `SE` and `MXE` are reported in genomic order on both strands.

| Event | Features emitted | Type |
| :--- | :--- | :--- |
| SE | (upstreamEE, exonStart_0base), (exonEnd, downstreamES), (upstreamEE, downstreamES) | junction ×3 |
| A5SS, + strand | (longExonEnd, flankingES), (shortEE, flankingES) | junction ×2 |
| A5SS, − strand | (flankingEE, longExonStart_0base), (flankingEE, shortES) | junction ×2 |
| A3SS, + strand | (flankingEE, longExonStart_0base), (flankingEE, shortES) | junction ×2 |
| A3SS, − strand | (longExonEnd, flankingES), (shortEE, flankingES) | junction ×2 |
| MXE | (upstreamEE, 1stExonStart_0base), (1stExonEnd, downstreamES), (upstreamEE, 2ndExonStart_0base), (2ndExonEnd, downstreamES) | junction ×4 |
| RI | (upstreamEE, downstreamES) **twice** | junction + retained_intron |

Retained introns are the one case where a single interval is emitted twice. The
retaining isoform is defined by the *absence* of the junction, so a junction test
alone cannot find it; the two types ask opposite questions of the same coordinates.
SUPPA `RI` events are handled the same way. MAJIQ marks the retained intron through
`IR coords`, which names one of the intervals already listed in `Junctions coords`;
that interval is emitted once, as `retained_intron`.

---

## Table S3. Labels for events and transcripts that cannot be analysed

Every input event reaches the output. Nine labels record why an event or transcript
produced no domain comparison, so coverage is auditable.

| Label | Level | Meaning |
| :--- | :--- | :--- |
| `no_gene_specified` | event | The event names no gene at all. Routine for annotation-free clustering; no database lookup was possible. |
| `gene_not_in_db` | event | A gene **was** named but its identifier is absent from the database. |
| `no_canonical_transcript` | event | The gene is present but has no canonical transcript to compare against. |
| `only_one_transcript` | event | The gene has a single transcript, so there is nothing to compare. |
| `feature_not_mapped` | event | A feature of the event matches no transcript of the gene. |
| `no_canonical_features` | event | No feature of the event is present in the canonical transcript. |
| `transcript_doesnt_have_features` | transcript | The transcript carries no feature of the event. |
| `no_unique_features` | transcript | The transcript carries only features the canonical transcript also has, so it does not differ within the event. |
| `no_domains_in_region` | transcript | The transcript was compared, but neither it nor the canonical transcript has an annotated domain in the region. A meaningful negative result. |

`no_gene_specified` and `gene_not_in_db` are deliberately distinct: the first is an
event with no gene to look up, the second a lookup that was possible and failed.

---

## Table S4. Classification of domain changes

Each identity group is classified from the number of member domains it contributes
on each side and the amino-acid length each side covers. **C** is the canonical
transcript, **T** the compared one. Lengths are the union of the member intervals,
so overlapping annotations of the same residues are not double-counted.

| C count | T count | Length comparison | Classification |
| :---: | :---: | :--- | :--- |
| 0 | ≥ 1 | — | `added_domain` |
| ≥ 1 | 0 | — | `dropped domain` |
| 1 | 1 | T = C | `same` |
| 1 | 1 | T > C | `longer` |
| 1 | 1 | T < C | `shorter` |
| 1 | > 1 | — | `split domain` |
| > 1 | 1 | — | `merged domain` |
| n | n (n > 1) | T = C | `same_domains` |
| n | n (n > 1) | T > C | `longer_domains` |
| n | n (n > 1) | T < C | `shorter_domains` |
| c | t (c < t, both > 1) | — | `increased_domain_number` |
| c | t (c > t, both > 1) | — | `reduced_domain_number` |

The `_domains` suffix marks the same length outcome reached over a group of equal,
multiple count rather than a single domain pair; the two forms are normalised
together for reporting.
