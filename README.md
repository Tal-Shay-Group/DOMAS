# DOMAS: Domain Oriented Mapping of Alternative Splicing

**Domain Oriented Mapping of Alternative Splicing (DOMAS)** is a computational framework
designed to bridge the gap between differential splicing events and their protein-level
effect. DOMAS accepts as input a list of differential splicing events and performs a
coordinate-aware mapping of these events onto protein domain architectures. It then
annotates each event with the affected domain(s) and classifies the effect as causing
domain loss, gain, alteration, truncation or elongation, compared to the protein encoded
by the canonical transcript.

## Key features

- **Input formats.** Output of several differential splicing tools: LeafCutter
  (Li et al., 2016), rMATS-turbo (Wang et al., 2024), MAJIQ (Vaquero-Garcia et al., 2023)
  and SUPPA `.ioe` files, plus a hadas-style comparative-splicing Excel file. To
  accommodate additional input formats, please write to us.
- **Output format.** A CSV listing the domains affected by each alternative splicing
  event and the class of the effect compared to the protein encoded by the canonical
  transcript — loss, gain, alteration, elongation or truncation. Each event can be
  linked to its visualization in DoChaP format (Gal-Oz et al., 2021).
- **Protein domain databases.** DOMAS builds upon the DoChaP database
  (Gal-Oz et al., 2021), which integrates domain annotation from several databases.
- Stay tuned — DOMAS will also be available as a web server soon.

## Installation and requirements

- Python 3.x
- Dependencies: `pandas`, `numpy`, `sqlite3`, `openpyxl` (for Excel input), `matplotlib`
- **Database:** requires access to a local instance of the **DoChaP DB**.
  See Gal-Oz et al. (2021) for installation instructions.

## Usage

Run the utility from the command line. `-format` selects the input reader, and the
remaining input arguments depend on it.

LeafCutter — takes a pair of `leafcutter_ds` output files and no `-input`:

```bash
python3 code/domas.py -format leafcutter -lc_sig <leafcutter_ds_cluster_significance.txt> -lc_effect <leafcutter_ds_effect_sizes.txt> -dochap <dochap.sqlite> -output_csv results.csv
```

rMATS-turbo — `-input` is the output *directory*:

```bash
python3 code/domas.py -format rmats -input <rmats_output_dir> -dochap <dochap.sqlite> -output_csv results.csv
```

MAJIQ — `-input` is a `voila tsv` file:

```bash
python3 code/domas.py -format majiq -input <voila.tsv> -dochap <dochap.sqlite> -output_csv results.csv
```

SUPPA — `-input` is a single `.ioe` file, or a directory of them (see `-ioe_pattern`):

```bash
python3 code/domas.py -format ioe -input <events.ioe> -dochap <dochap.sqlite> -output_csv results.csv
```

hadas-format Excel:

```bash
python3 code/domas.py -format hadas -input <comparative_splicing.xlsx> -dochap <dochap.sqlite> -output_csv results.csv
```

### Input parameters

| Parameter | Description |
| :---- | :---- |
| `-format` | Input format: `hadas`, `ioe`, `leafcutter`, `rmats` or `majiq`. Required. |
| `-input` | Input path for `-format hadas/ioe/rmats/majiq`. Not used with `-format leafcutter`. |
| `-lc_sig` | Path to `leafcutter_ds_cluster_significance.txt` (only with `-format leafcutter`). |
| `-lc_effect` | Path to `leafcutter_ds_effect_sizes.txt` (only with `-format leafcutter`). |
| `-dochap` | Path to the DoChaP sqlite database. Required. |
| `-output_csv` | Destination path for the generated results (default `junctions_analysis.csv`). |
| `-ioe_pattern` | Filename regex used to find `.ioe` files when `-input` is a directory. |
| `-examples_per_event` | Per event type, keep only this many example clusters (0 keeps every cluster); only when `-input` is a directory with `-format ioe`. |
| `-gene_ids` | Comma-separated gene symbols to generate PDFs for (only with `-format hadas` and `-pdf`). |
| `-num_workers` | Number of parallel worker processes (default: the machine's CPU count). |
| `-max_clusters` | If > 0, analyze only the first N clusters (sorted). 0 (default) means no limit. |
| `-stats_out_dir` | Directory for the stats report (default: alongside `-output_csv`). |

### Optional behaviour flags

The defaults are what a normal run wants; each flag turns one of them off.

| Flag | Effect |
| :---- | :---- |
| `-pdf` | Generate a per-gene PDF. Off by default — a full-scale run would otherwise produce one PDF per gene. Only honored with `-format hadas`. |
| `-no_representative_domains` | Use `DomainEvent`/`DomainType` only. By default domains come from the `RepresentativeDomains` table where available, falling back per protein. |
| `-no_stats` | Skip the statistics report (event distribution, domain frequency, etc.) that is otherwise generated for `-output_csv` after the run. |
| `-keep_non_comparable` | Keep rows for non-comparable transcripts (`gene_not_in_db`, `junction_not_mapped`, `no_unique_junctions`, …). By default the output CSV holds only transcripts actually compared to the canonical transcript. |

### A note on rMATS event types

DOMAS reads the `SE`, `A5SS`, `A3SS` and `MXE` `[EventType].MATS.JC.txt` files.
`RI.MATS.JC.txt` is deliberately not read: a retained-intron event yields a single
junction, and a single-junction event can never produce a transcript that differs from
the canonical one within the event, so such events are structurally unanalysable.

## Output format

The resulting CSV file provides:

- **Domain and change columns** — the affected domain identifier and the predicted change
  (e.g. gain / loss), with the canonical and compared domain lengths and counts.
- **Domain descriptions** — functional description of the affected domain, from DoChaP.
- **Identification columns taken from the input** — cluster, gene symbol, gene ID,
  species and the coordinates of the alternative splicing event.
- **`is_longest_cds` / `is_most_like_canonical`** — per-transcript flags marking the
  transcript with the longest annotated CDS and the transcript most like the canonical
  one. Every comparable transcript is reported; filter on a flag for a single-isoform
  view. A transcript may carry both flags, one, or neither.

## References

Gal-Oz, S. T., Haiat, N., Eliyahu, D., Shani, G., & Shay, T. (2021). DoChaP: the domain
change presenter. *Nucleic Acids Research, 49*(W1), W162–W168.

Li, Y. I., Knowles, D. A., & Pritchard, J. K. (2016). LeafCutter: annotation-free
quantification of RNA splicing. *bioRxiv*, 044107.

Vaquero-Garcia, J., Aicher, J. K., Jewell, S., Gazzara, M. R., Radens, C. M., Jha, A.,
Norton, S. S., Lahens, N. F., Grant, G. R., & Barash, Y. (2023). RNA splicing analysis
using heterogeneous and large RNA-seq datasets. *Nature Communications, 14*(1), 1230.

Wang, Y., Xie, Z., Kutschera, E., Adams, J. I., Kadash-Edmondson, K. E., & Xing, Y.
(2024). rMATS-turbo: an efficient and flexible computational tool for alternative
splicing analysis of large-scale RNA-seq data. *Nature Protocols, 19*(4), 1083–1104.
