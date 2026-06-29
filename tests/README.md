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

`test_flags.py` points at a hardcoded DB path
(`/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-web/DB_merged.sqlite`).
Update `DB_PATH` near the top of the file if the database moves.

### What's covered

- **`test_*_all_flag_combinations`** (8 tests): runs
  `analyze_junctions()` end-to-end against both fixture files, for every
  combination of `use_longest_cds` and `restrict_pdf_to_comparable`.
  Verifies the run completes without error, produces a results CSV and at
  least one PDF, and - when `use_longest_cds=True` - that no cluster ever
  compares more than one non-canonical transcript to the canonical one.
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

## Reference outputs (manual inspection / regression reference)

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
scratch each time it's run. Use it to manually eyeball results after a code
change, or to diff against a previously-committed "golden" case to confirm a
fix only changed what you expect (e.g. `diff <(sort old/results.csv) <(sort
new/results.csv)`).

Only `hadas_xlsx__longestcds_True__restrict_True/` is currently committed to
git, as a golden reference point.
