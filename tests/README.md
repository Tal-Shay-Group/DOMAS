# Tests

Tests for the junction/domain comparison pipeline (`junction_analisys.py`,
`generate_gene_pdf.py`), run against the real DoChaP database.

## Fixtures

- `ioe_example_junctions.csv` - plain CSV junctions file (human genes,
  `ENSG*` ids), read with `hadas_format=False`.
- `short_H_vs_M_HN6.xlsx` - human-vs-mouse comparative splicing table, read
  with `hadas_format=True` (parsed via `domas.hadas_read_input_file`, which
  needs a DB connection to resolve mouse gene symbols to ensembl ids).

Both are passed as `junctions_csv` to `JunctionsAnalysis.analyze_junctions()`.

## Running the tests

```bash
cd tests
python3 -m pytest test_flags.py -v
```

Drop `-v` for terse output, or run a subset by name:

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

- **`test_*_all_flag_combinations`** (8 tests): runs
  `analyze_junctions()` end-to-end against both fixture files, for every
  combination of `use_longest_cds` and `restrict_pdf_to_comparable`.
  Verifies the run completes without error, produces a results CSV and at
  least one PDF, and - when `use_longest_cds=True` - that no cluster ever
  compares more than one non-canonical transcript to the canonical one.
  Output goes to a pytest-managed temp dir (`tmp_path`) and is not kept.
- **`test_compare_against_reference_outputs`** (8 tests, automatic golden
  comparison): for each of the same 8 cases, runs `analyze_junctions()` and
  writes its CSV/PDF output to a persistent directory,
  `tests/generated_outputs/<case_name>/`, then compares the generated
  `results.csv` row-for-row against the golden reference at
  `tests/reference_outputs/<case_name>/results.csv` (order-independent).
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
standalone script that runs all 8 flag combinations against both fixture
files and writes the full, persistent output (results.csv + every generated
PDF) into `reference_outputs/<case_name>/`, e.g.:

```
reference_outputs/
├── ioe_csv__longestcds_False__restrict_False/
├── ioe_csv__longestcds_False__restrict_True/
├── ioe_csv__longestcds_True__restrict_False/
├── ioe_csv__longestcds_True__restrict_True/
├── hadas_xlsx__longestcds_False__restrict_False/
├── hadas_xlsx__longestcds_False__restrict_True/
├── hadas_xlsx__longestcds_True__restrict_False/
└── hadas_xlsx__longestcds_True__restrict_True/
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

Only `hadas_xlsx__longestcds_True__restrict_True/` is currently committed to
git, as a golden reference point. The other 7 cases exist locally once
generated but aren't tracked yet; `test_compare_against_reference_outputs`
skips (rather than fails) any case whose reference isn't present.

## Other files

- `conftest.py` - registers the `--keep-test-output` and `--db-path` pytest
  options, and the `keep_test_output` / `db_path` fixtures built from them.
