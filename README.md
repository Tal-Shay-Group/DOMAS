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
| `-stats_out_dir` | Directory for the `-stats` report (default: alongside `-output_csv`). |

### Optional behaviour flags

The defaults are what a normal run wants; each flag turns one of them off.

| Flag | Effect |
| :---- | :---- |
| `-pdf` | Generate a per-gene PDF. Off by default — a full-scale run would otherwise produce one PDF per gene. Only honored with `-format hadas`. |
| `-no_representative_domains` | Use `DomainEvent`/`DomainType` only. By default domains come from the `RepresentativeDomains` table where available, falling back per protein. |
| `-stats` | Also generate the statistics report (event distribution, domain frequency, etc.) for `-output_csv` after the run. Off by default. |
| `-keep_non_comparable` | Keep rows for non-comparable transcripts (`gene_not_in_db`, `junction_not_mapped`, `no_unique_junctions`, …). By default the output CSV holds only transcripts actually compared to the canonical transcript. |

### A note on rMATS event types

DOMAS reads all five `[EventType].MATS.JC.txt` files: `SE`, `A5SS`, `A3SS`, `MXE` and
`RI`. A retained-intron record names one interval, the spliced form, so it is emitted
as two features at the same coordinates — one matched by adjacency (the intron is
spliced out) and one by containment (a single exon holds it). That is what lets a
retaining transcript carry a feature the canonical one lacks.

## Example runs

`run_examples.sh` runs DOMAS once per input format against the fixtures in `tests/`
and compares each result against the reference stored in `tests/run_examples/`. The
DoChaP database is not in this repository, so pass its path:

```bash
./run_examples.sh /path/to/DB_merged.sqlite [output_dir]
```

`output_dir` defaults to `./run_examples_output`, and each format writes `<format>.csv`
there. DOMAS creates the directory if it does not exist. The whole set takes about two
minutes; add `-stats` to a command to get the statistics report as well.

Each run also writes `domas.log` beside its `-output_csv` and prints nothing to the
console. The four examples share one output directory, so only the last run's log
survives.

The script is four `domas.py` command lines, one per format — copy the one you need and
swap in your own input.

| Example | Input | Rows |
| :--- | :--- | ---: |
| `ioe` | `tests/ioe/` — 20 SUPPA events per event type | 265 |
| `rmats` | `tests/rmats/` — the five `[EventType].MATS.JC.txt` files | 344 |
| `majiq` | `tests/majiq/NveB_Mono_voila.txt` — a `voila tsv` | 372 |
| `leafcutter` | `tests/leafcutter/` — the `leafcutter_ds` output pair | 316 |

Each invocation uses the minimal arguments its format needs: `-format`, the input, and
`-specie`, which is required because none of these files carries a species field. The
three large fixtures also pass `-max_clusters 200` so an example finishes in seconds —
drop it to analyse the whole fixture.

Each reference is the exact CSV its run produces, so a mismatch means the analysis
changed. The script checks them at the end with `tests/compare_run_example.py`, which
sorts both sides first — row order is not part of the result, since clusters are
analysed in parallel and returned as they finish.

## Output format

The resulting CSV file provides:

- **Domain and change columns** — the affected domain identifier and the predicted change
  (e.g. gain / loss), with the canonical and compared domain lengths and counts.
- **Domain descriptions** — functional description of the affected domain, from DoChaP.
- **Identification columns taken from the input** — cluster, gene symbol, gene ID,
  species and the coordinates of the alternative splicing event.
- **`is_longest_cds` / `is_most_like_canonical`** — written only with
  `-write_all_comparable`, which compares every comparable transcript and keeps a row
  for each. The flags mark the transcript with the longest annotated CDS and the one
  most like the canonical; a transcript may carry both, one, or neither. By default
  only the transcript those rules select is compared, so each cluster holds one
  comparison and the two columns are omitted.

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
