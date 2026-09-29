# Paper examples

The two use cases reported in the DOMAS paper, as input, a script and the results.
Everything here regenerates from the files in `input/`, so the numbers in the paper
can be checked without assembling anything.

| | Use case 1 — immune | Use case 2 — placenta |
| --- | --- | --- |
| comparison | human monocytes vs CD4 T cells | preeclamptic vs preterm birth placentas |
| dataset | GSE60424 (Linsley et al., 2014) | GSE114691 |
| caller | LeafCutter | LeafCutter |
| script | `run_immune.sh` | `run_placenta.sh` |
| results | `results/immune/` | `results/placenta/` |

## Running them

The DoChaP database is not in this repository — download it and pass its path:

```bash
./run_immune.sh   /path/to/DB_merged.sqlite
./run_placenta.sh /path/to/DB_merged.sqlite
```

Each takes well under a minute and overwrites `results/<case>/`. An output directory
can be given as a second argument to write somewhere else instead.

Both scripts pass `-max_clusters 100`, which is what the paper reports and what the
web server applies: the 100 events with the lowest `p.adjust`, not the first 100 in
the file. They also pass `-no_excel`, so the results stay CSV and can be read and
diffed as text — drop it and DOMAS writes `.xlsx` instead, its default, with every
gene symbol linked to that gene's DoChaP page.

## What the results hold

- `compared.csv` — one row per alternative-transcript group per domain, ranked by
  expected effect.
- `non_compared.csv` — the events that could not be compared, and the per-transcript
  and per-junction notes explaining why a transcript or junction was not used.
- `compared_summary.txt` — how much of the input was analysed, how much of it reached
  a comparison, and the outcome counts.

## The figures these reproduce

| | immune | placenta |
| --- | --- | --- |
| input events | 100 | 100 |
| event–gene pairs | 108 | 119 |
| comparable | 98 | 70 |
| with at least one domain change | 58 | 37 |
| group × domain rows | 185 | 168 |
| non-coding alternative | 15 | 13 |
| domain gain | 0 | 1 |
| domain loss | 60 | 46 |
| longer domain | 4 | 6 |
| shorter domain | 12 | 12 |
| no domain change | 28 | 28 |
| no domains in region | 66 | 62 |

The scripts pass `-num_workers 1`, so a rerun reproduces these files exactly - only
the `Date` line of the summary changes. (With several workers each event's rows are
written as its worker finishes, which reorders the file without changing the result.)

Produced with DOMAS 1.0.0 against the DoChaP build listed in Table S1 of the paper.
Domain annotation changes with a database rebuild, so a different build can move
these counts.
