# DOMAS: Domain Oriented Mapping of Alternative Splicing

**Domain Oriented Mapping of Alternative Splicing (DOMAS)** is a computational framework
designed to bridge the gap between differential splicing events and their protein-level
effect. DOMAS accepts as input a list of differential splicing events and performs a
coordinate-aware mapping of these events onto protein domain architectures. It then
annotates each event with the affected domain(s) and classifies the effect as causing
domain loss, gain, truncation or elongation, compared to the protein encoded
by the canonical transcript.

## Key features

- **Input formats.** Output of several differential splicing tools: LeafCutter
  (Li et al., 2016), rMATS-turbo (Wang et al., 2024), MAJIQ (Vaquero-Garcia et al., 2023)
  and SUPPA (Trincado et al., 2018). To accommodate additional input formats, please write to us.
- **Output format.** Two Excel files, listing the domains affected by each alternative splicing
  event and the class of the effect compared to the protein encoded by the canonical
  transcript — loss, gain, elongation or truncation. Each event can be
  linked to its visualization in DoChaP format (Gal-Oz et al., 2021). The other file list
  the events that cannot be compared.
- **Protein domains.** DOMAS uses **Pfam** domains only. It builds upon the DoChaP
  database (Gal-Oz et al., 2021), which integrates annotation from several signature
  databases, and keeps the Pfam signatures: one region matched by Pfam, SMART and CDD
  would otherwise be counted as three overlapping domains, and there is no entry type
  to rank them by. Two instances of the same Pfam accession overlapping by at least
  half of the shorter one are treated as a single domain; below that they are separate
  instances, so tandem repeats are counted individually.
- **Web server.** The same analysis in the browser, with nothing to install and no
  local database — see [Web server](#web-server) below.

## Web server

**[dochap.bgu.ac.il/#!/domas](https://dochap.bgu.ac.il/#!/domas)** runs the whole
analysis in the browser — no installation, and no local copy of the DoChaP database.
Each upload is capped at the 100 most significant events the input names, so filter
longer lists first, or use the stand-alone version below, which has no cap.

**Running a bundled example.** Under *Try it with sample data*, pick a format
(LeafCutter, rMATS or MAJIQ) and press **Run example**. These are real human data and
need no files of your own. *Download example input* hands you the same files as a zip,
named the way the upload form expects them, which is the quickest way to see what your
own input should look like.

**Running your own files.** Choose an **Input format**, pick the **Species** the events
were called against (human, mouse or rat — none of the input formats carries a species
field, so it has to be stated), select the file(s) and press **Process**:

| Input format | What to upload |
| --- | --- |
| LeafCutter | both files: `..._cluster_significance.txt` and `..._effect_sizes.txt` |
| rMATS | the `[EventType].MATS.JC.txt` files (SE, A3SS, A5SS, MXE, RI) |
| MAJIQ | the `voila tsv` output |
| SUPPA | a single events `.ioe` file |

A run takes a minute or two. The results appear as a table below the form, where each
gene symbol links to that gene's DoChaP visualization, with only the canonical
transcript and the compared alternative transcript shown (that is, not hidden) and 
**Download results** gives a zip file holding the compared events, the events that could not be compared, 
and the run summary. The summary is also shown on the page.

## Installation and requirements

- Python 3.9 or higher
- Dependencies: see `requirements.txt`

  ```bash
  pip install -r requirements.txt
  ```
- **Database:** requires access to a local instance of the **DoChaP DB**.
  See Gal-Oz et al. (2021) for installation instructions.

## Usage

Run the utility from the command line. `-format` selects the input reader, and the
remaining input arguments depend on it. `-species` is required throughout: none of
these formats carries a species field, and DOMAS stops if the gene ids turn out to
belong to a different species than the one stated.

The file `run_examples.sh` contains examples of calling DOMAS using the `domas.py` 
command line for each of the input formats — read it
for the invocation you need and swap in your own input file names. It runs DOMAS once per format
against the fixtures in `tests/` and checks each result against the reference stored
in `tests/run_examples/`. The DoChaP database is not in this repository, so pass its
path:

```bash
./run_examples.sh /path/to/DB_merged.sqlite [output_dir]
```

The whole set takes about two minutes.

### Reproducing the paper's examples

`paper_examples/` holds the two use cases reported in the paper — the immune and the
placenta comparison — as the input files, a run script each, and the results those
scripts produce. Running them needs only the database:

```bash
cd paper_examples
./run_immune.sh   /path/to/DB_merged.sqlite
./run_placenta.sh /path/to/DB_merged.sqlite
```

Each takes well under a minute and writes into `results/immune/` or
`results/placenta/`, overwriting what is committed there. Pass an output directory as
a second argument to write somewhere else and keep the committed results to compare
against.
Note that the input files can also be run on the web server with no installation, to
receive the results for the 100 most significant events, as described in the paper.

Read the output in this order:

1. `compared_summary.txt` — start here. It gives the headline (how many event-gene
   pairs changed at least one domain), then the input as DOMAS read it, how many of
   the events were comparable and why the rest were not, and the count of every
   outcome. Every figure the paper reports for these two runs is in this file or
   derived from the next one.
2. `compared.csv` — one row per alternative-transcript group:domain_combination, 
   ranked by expected effect: non-coding products first, then domain gains, losses,
   length changes by descending percentage, then the unchanged, then the rows whose
   window held no domain. Reading from the top is reading the candidates in priority
   order.
3. `non_compared.csv` — the events that could not be compared, with the reason, plus
   the per-transcript and per-junction notes saying why a particular transcript or
   junction was not used. Almost all of its rows are those notes, not failures.

`paper_examples/README.md` tabulates the figures the two runs should produce, so a
rerun can be checked line by line against the paper.

**A difference from the paper worth knowing.** Both scripts pass `-no_excel`, which
the paper does not mention because it is not part of the method: DOMAS writes
`compared.xlsx` and `non_compared.xlsx` by default, and `-no_excel` asks for the CSVs
instead. The results are the same rows either way. CSV is committed here because it
can be read and diffed as text in a pull request, where a workbook cannot. Drop the
flag from either script to get exactly the file set the paper describes — the
workbooks, with each gene symbol hyperlinked to that gene's DoChaP page, which the
CSV cannot carry.

Both scripts also pass `-num_workers 1`. With several workers each event's rows are
written as its worker finishes, so the row order varies between runs even though the
result does not; one worker makes the output reproducible. A rerun therefore matches
the committed results exactly, apart from the `Date` line of the summary. The cost is
negligible here - the analysis is a fraction of a second either way, and the rest of
the run is reading the database, which is not parallelised.

Every parameter and flag, with its default, is documented in the CLI's own help:

```bash
python3 code/domas.py -h
```

The defaults are what a normal run wants; each flag either turns one of them off or
asks for something extra, such as the statistics report.

## Output format

Each results table is saved as an Excel workbook (`.xlsx`), so it opens in a spreadsheet
without an import step. `-no_excel` leaves it as CSV instead, and so does a table with more
rows than an Excel sheet holds (1,048,575) — that case is reported in the log and the CSV
stands, so a result is never made unreadable for being too large. 

The analysis writes CSV throughout and the command converts it to Excel once at the end. Anything that
reads a run's output back — `results_stats`, `compare_results_csv`, `analyze_results.py`,
the test suite — reads a CSV, so pass `-no_excel` when a run feeds one of those. `-stats`
needs nothing: it runs before the conversion.

The resulting table provides:

- **Domain and change columns** — the affected domain identifier and the predicted change
  (e.g. gain / loss), with the canonical and compared domain lengths and counts.
- **`length_change_pct`** — how much the domain's length changed, as a percentage of its
  length in the canonical transcript: `(|alternative − canonical|) / canonical × 100`, to two
  decimals. A magnitude, so a domain cut in half and one grown by half both read `50.0`
  and the outcome column says which happened. Empty where the change cannot be measured:
  rows carrying no domain lengths at all (`non_coding_alternative`, `gained_protein`), and
  any row whose canonical domain length is zero — the `domain_gain` rows where the
  canonical side holds no instance of the domain at all. A `domain_gain` that goes from one
  instance to two does have a canonical length, and so does get a percentage.
  This is the number the default row order ranks the `longer_domain` and `shorter_domain`
  blocks by, so sorting on the column reproduces the order the file is already in.
- **Domain descriptions** — functional description of the affected domain, from DoChaP.
- **Identification columns taken from the input** — cluster, gene symbol, gene ID,
  species and the coordinates of the alternative splicing event.
- **`is_longest_cds` / `is_most_like_canonical`** — written only with
  `-write_all_comparable`, which compares every comparable transcript and keeps a row
  for each. The flags mark the transcript with the longest annotated CDS and the one
  most like the canonical; a transcript may carry both, one, or neither. By default
  only the transcript those rules select is compared, so each cluster holds one
  comparison and the two columns are omitted.
- **`rank` / `canonical_junction_in_cds` / `alternative_junction_in_cds`** — written
  only with `-extra_columns`, and available for every input format. `rank` names the exons of
  the canonical transcript that the event's junctions join — `E2_E4` where the event
  skips exon 3, `E11_E13Last` where it reaches the final exon, and `*` for a splice
  site that is no exon edge of the canonical, which is what an alternative-splice-site
  event has. The other two say whether those junctions fall inside the coding sequence
  of the canonical and of the compared transcript respectively: `yes` when every one of
  them does, `no` when none does (the event is in a UTR), `partial` when a junction
  straddles the start or stop codon or the group is mixed, and `no_cds` for a
  transcript with no annotated protein.

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
