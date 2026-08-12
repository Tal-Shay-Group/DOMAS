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
  and SUPPA `.ioe` files. To accommodate additional input formats, please write to us.
- **Output format.** A CSV listing the domains affected by each alternative splicing
  event and the class of the effect compared to the protein encoded by the canonical
  transcript — loss, gain, alteration, elongation or truncation. Each event can be
  linked to its visualization in DoChaP format (Gal-Oz et al., 2021).
- **Protein domain databases.** DOMAS builds upon the DoChaP database
  (Gal-Oz et al., 2021), which integrates domain annotation from several databases.
- Also available as a web server in dochap.bgu.ac.il. Look for the DOMAS page.

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

`run_examples.sh` is one working `domas.py` command line per input format — read it
for the invocation you need and swap in your own input. It runs DOMAS once per format
against the fixtures in `tests/` and checks each result against the reference stored
in `tests/run_examples/`. The DoChaP database is not in this repository, so pass its
path:

```bash
./run_examples.sh /path/to/DB_merged.sqlite [output_dir]
```

The whole set takes about two minutes.

Every parameter and flag, with its default, is documented in the CLI's own help:

```bash
python3 code/domas.py -h
```

The defaults are what a normal run wants; each flag either turns one of them off or
asks for something extra, such as the statistics report.

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
- **`rank` / `c_junction_in_cds` / `t_junction_in_cds`** — written only with
  `-extra_columns`, and available for every input format. `rank` names the exons of
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
