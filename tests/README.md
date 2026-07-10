# Tests

Tests for the junction/domain comparison pipeline (`junction_analisys.py`,
`generate_gene_pdf.py`). `test_flags.py` runs against the real DoChaP
database; `test_junction_analysis.py` is DB-free unit tests for
`junction_analisys.py`'s pure functions.

## Fixtures

- `ioe_example_junctions.csv` - plain CSV junctions file (human genes,
  `ENSG*` ids), read with `alternative_splicing.read_junctions_csv()`.
- `short_H_vs_M_HN6.xlsx` - human-vs-mouse comparative splicing table, read
  with `alternative_splicing.hadas_read_input_file()` (needs a DB connection
  to resolve mouse gene symbols to ensembl ids).

Both are read into a DataFrame first, then passed as `df_junctions` to
`JunctionsAnalysis.analyze_junctions()` - reading junctions from a file is
alternative_splicing.py's job, not junction_analisys.py's.

## Running the tests

```bash
cd tests
python3 -m pytest -v
```

Runs both `test_flags.py` and `test_junction_analysis.py`. Drop `-v` for
terse output, or run a subset by name:

```bash
python3 -m pytest test_flags.py -v -k "ioe_csv and True"
```

By default the tests run against
`/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-web/DB_merged.sqlite`.
Point them at a different database with `--db-path`:

```bash
python3 -m pytest test_flags.py -v --db-path /path/to/other/DB_merged.sqlite
```

For local runs against the usual database, just omit the flag - the default
already points at it.

### What's covered

Every comparable transcript is always compared to canonical now - there's no
tie-break flag that narrows a cluster down to one transcript. Each row is
tagged `is_longest_cds`/`is_most_like_canonical` instead. `FLAG_COMBINATIONS`
varies `(restrict_pdf_to_comparable, use_representative_domains)`.

- **`test_*_all_flag_combinations`** (8 tests): runs
  `analyze_junctions()` end-to-end against both fixture files, for every
  combination of `restrict_pdf_to_comparable` and `use_representative_domains`.
  Verifies the run completes without error, produces a results CSV and at
  least one PDF, and that each tie-break rule tags at most one transcript per
  cluster. Output goes to a pytest-managed temp dir (`tmp_path`) and is not kept.
- **`test_compare_against_reference_outputs`** (8 tests, automatic golden
  comparison): for each of the same 4 cases x 2 fixtures, runs `analyze_junctions()` and
  writes its CSV/PDF output to a persistent directory,
  `tests/generated_outputs/<case_name>/`, then compares the generated
  `results.csv` row-for-row against the golden reference at
  `tests/reference_outputs/<case_name>/results.csv` (order-independent), and
  compares each generated PDF's extracted text against the golden manifest at
  `tests/reference_outputs/<case_name>/pdf_text_reference.json` (see
  "PDF text reference" below).
  Skips with a clear message if no reference exists yet for that case.
  The generated output directory is deleted **only after a successful
  comparison** - on failure it's left in place for inspection. Pass
  `--keep-test-output` to keep it even on success, e.g. to look at the PDFs
  by hand:

  ```bash
  python3 -m pytest test_flags.py -v --keep-test-output
  ```
- **`test_comparable_transcript_ids_*`** (2 tests): unit tests for
  `JunctionsAnalysis._comparable_transcript_ids`, which derives "canonical +
  actually-compared transcripts" from a cluster's results.
- **`test_transcript_matches_ids_checks_both_id_columns`**: unit test for
  `GeneVisualization._transcript_matches_ids`.
- **`test_create_pdf_transcript_ids_filters_transcripts`** /
  **`test_create_pdf_transcript_ids_filters_real_gene`**: confirm
  `create_pdf(transcript_ids=...)` actually restricts which transcripts get
  drawn (the second uses a real 12-transcript gene and checks the resulting
  PDF page count).

## Two output directories - don't confuse them

- **`reference_outputs/<case_name>/`** - the **golden** results, checked into
  git intentionally. This is what `test_compare_against_reference_outputs`
  compares against. Only update these on purpose, after manually confirming a
  change is correct (see below).
- **`generated_outputs/<case_name>/`** - **scratch** output produced by
  `test_compare_against_reference_outputs` on every run. Deleted
  automatically after a successful comparison; kept (for that run) on
  failure, or always with `--keep-test-output`. Not meant to be committed.

## Reference outputs (golden, manual inspection / regression reference)

`generate_reference_outputs.py` is **not** part of the pytest run - it's a
standalone script that runs every `(restrict_pdf_to_comparable,
use_representative_domains)` combination against both fixture files and
writes the full, persistent output (results.csv + every generated PDF) into
`reference_outputs/<case_name>/`, e.g.:

```
reference_outputs/
├── ioe_csv__restrict_False__representative_False/
├── ioe_csv__restrict_True__representative_False/
├── ioe_csv__restrict_False__representative_True/
├── ioe_csv__restrict_True__representative_True/
├── hadas_xlsx__restrict_False__representative_False/
├── hadas_xlsx__restrict_True__representative_False/
├── hadas_xlsx__restrict_False__representative_True/
└── hadas_xlsx__restrict_True__representative_True/
```

Regenerate them with:

```bash
python3 generate_reference_outputs.py
```

This wipes and rebuilds the entire `reference_outputs/` directory from
scratch each time it's run. After a deliberate algorithm change, use it to
manually eyeball the new results and diff them against the previously
committed golden case before re-committing
(e.g. `diff <(sort old/results.csv) <(sort new/results.csv)`) - don't
overwrite a golden reference without checking the diff first.

`results.csv` and `pdf_text_reference.json` are committed for all 8 cases.
The full PDFs themselves are only committed for
`hadas_xlsx__restrict_True__representative_False/`, as a reference point for
manual/visual inspection - the other 7 cases' PDFs exist locally once
generated but aren't tracked (regenerating writes them, but don't `git add`
them). `test_compare_against_reference_outputs` skips (rather than fails)
any case whose `results.csv` or `pdf_text_reference.json` isn't present.

## PDF text reference

Raw PDF bytes aren't diffable as a golden reference: matplotlib embeds a
CreationDate (and similar metadata) in every PDF it writes, so re-running the
exact same analysis produces a byte-different file even when the visible
content is identical - confirmed by regenerating a case and diffing old vs
new: same file sizes, only a handful of differing bytes per PDF, no content
change.

Instead, `tests/pdf_text_utils.py` extracts each PDF page's text with
PyPDF2 and stores it as `pdf_text_reference.json` (one file per case,
mapping PDF filename to a list of per-page text) - gene/transcript/protein
ids, event names, and domain labels are all real text elements, so this
catches real content regressions (wrong transcript compared, dropped
domain, mislabeled event, ...) without false-failing on timestamp noise.
It won't catch a purely visual/layout bug that doesn't change the extracted
text (e.g. a domain drawn at the wrong x-position but still listed in the
per-page text) - an image-diff approach would be needed for that, at the
cost of being flakier across machines (font hinting/anti-aliasing differ by
OS).

`generate_reference_outputs.py` writes `pdf_text_reference.json` alongside
`results.csv` automatically. If you regenerate references by hand instead,
rebuild the manifest with:

```python
from pdf_text_utils import write_pdf_text_manifest
write_pdf_text_manifest('reference_outputs/<case_name>')
```

## Other files

- `conftest.py` - registers the `--keep-test-output` and `--db-path` pytest
  options, and the `keep_test_output` / `db_path` fixtures built from them.
- `pdf_text_utils.py` - builds/writes the `pdf_text_reference.json` manifest
  used by both `generate_reference_outputs.py` and
  `test_compare_against_reference_outputs` (see "PDF text reference" above).
