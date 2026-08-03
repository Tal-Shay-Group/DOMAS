"""
Download the raw source files for the DOMAS enrichment database.

This module only fetches files onto disk (no parsing) so it can be re-run
cheaply: anything already present is skipped unless -force_download is given.
build.py turns what this downloads into enrichment.sqlite.

Sources (ECOD deliberately excluded - not consumed by the analysis, and large):
  - uniprot : reviewed (SwissProt) reference-proteome flat file per species,
              parsed later into uniprot_feature / uniprot_alias
  - afdb    : AlphaFold DB per-species proteome tar, parsed into afdb_plddt
  - ensembl : Ensembl peptide FASTA per species, parsed into ensembl_sequence
  - pfam    : Pfam-A HMM library; stays a file on disk (hmmscan reads it),
              build.py only runs hmmpress and records provenance

Nothing here is licensed NC/ND: UniProt = CC BY 4.0, AlphaFold DB = CC-BY 4.0,
Ensembl = open, Pfam = CC0. (APPRIS, which is CC BY-NC-ND, is not used.)
"""
import os
import re
import sys
import urllib.request

# DoChaP species label -> the per-source identifiers needed to build URLs.
#   up      : UniProt reference-proteome id "<proteomeID>_<taxon>"
#   afdb    : AlphaFold DB proteome tar basename (without .tar)
#   ensembl : Ensembl species directory name under current_fasta/
# AlphaFold DB re-versioned its proteome archives to v6; the v4 basenames used here
# until now 404 (the /latest/ directory serves *_v6.tar only). Sizes at v6: human
# 4.8 GB, mouse 3.5 GB, rat 3.5 GB.
SPECIES = {
    "H_sapiens":    {"up": "UP000005640_9606",  "afdb": "UP000005640_9606_HUMAN_v6",  "ensembl": "homo_sapiens"},
    "M_musculus":   {"up": "UP000000589_10090", "afdb": "UP000000589_10090_MOUSE_v6", "ensembl": "mus_musculus"},
    "R_norvegicus": {"up": "UP000002494_10116", "afdb": "UP000002494_10116_RAT_v6",   "ensembl": "rattus_norvegicus"},
}

# DoChaP also carries D_rerio and X_tropicalis. AlphaFold DB publishes a zebrafish
# proteome archive (UP000000437_7955_DANRE_v6.tar, 4.6 GB) but NONE for X. tropicalis,
# so structural enrichment is reachable for zebrafish and unreachable in bulk for frog.
# Neither is wired into the download/build paths yet (uniprot/ensembl entries would be
# needed too).
_AFDB_EXTRA = {"D_rerio": "UP000000437_7955_DANRE_v6"}

_UNIPROT_BASE = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
                 "knowledgebase/reference_proteomes/Eukaryota")
_AFDB_BASE = "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest"
_ENSEMBL_BASE = "https://ftp.ensembl.org/pub/current_fasta"
_PFAM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

# Per-source metadata surfaced in enrichment.sqlite's source_meta table.
LICENSES = {
    "uniprot": "CC BY 4.0",
    "afdb":    "CC-BY 4.0",
    "ensembl": "Ensembl terms (open, no restrictions)",
    "pfam":    "CC0 1.0",
    "alphamissense": "CC BY 4.0",
    "gnomad":  "CC0 1.0 (gnomAD terms, freely available)",
}
ALL_SOURCES = ("uniprot", "afdb", "ensembl", "pfam")
# Opt-in sources: large and/or not part of a default provision run.
OPTIN_SOURCES = ("alphamissense", "gnomad")

# gnomAD v2.1.1 per-gene LoF constraint (LOEUF = oe_lof_upper, pLI). Small (~5 MB
# bgzip); build.py maps gene symbol -> UniProt accession via the reference
# proteome and stores gene_constraint. Gene-level loss-intolerance is an
# orthogonal, non-circular signal for the calibrated impact probability.
_GNOMAD_URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/"
               "2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")

# AlphaMissense per-variant pathogenicity predictions (human canonical UniProt
# isoforms; ~1.2 GB gzipped, 216M variant rows). build.py aggregates them to a
# per-residue mean pathogenicity used as a functional-constraint signal.
#
# LICENSE: the AlphaMissense *predictions* are CC BY 4.0. DeepMind relicensed
# them from the original CC BY-NC-SA 4.0 after release, so they are usable
# (incl. commercially) with attribution (Cheng et al., Science 2023). DeepMind's
# current terms: https://github.com/google-deepmind/alphamissense#readme
# NOTE: the Zenodo record below still shows *stale* CC BY-NC-SA 4.0 metadata -
# the operative license is CC BY 4.0 per DeepMind's own repo. (Code is Apache
# 2.0; model weights are not released - we use only the prediction file.)
_ALPHAMISSENSE_URL = ("https://zenodo.org/records/8208688/files/"
                      "AlphaMissense_aa_substitutions.tsv.gz")


def _subdir(data_dir, source):
    d = os.path.join(data_dir, source)
    os.makedirs(d, exist_ok=True)
    return d


def _download(url, dest, force=False):
    """Stream `url` to `dest`, skipping if already present. Uses a .part temp
    so an interrupted download never leaves a truncated file in place."""
    if os.path.exists(dest) and not force:
        print(f"  skip (exists): {os.path.basename(dest)}")
        return dest
    tmp = dest + ".part"
    print(f"  GET {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "DOMAS-enrichment/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as out:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    os.replace(tmp, dest)
    print(f"  saved -> {dest} ({os.path.getsize(dest)/1e6:.1f} MB)")
    return dest


def _resolve_ensembl_pep_url(ensembl_species):
    """Find the *.pep.all.fa.gz in an Ensembl species pep/ directory. The
    assembly name in the filename changes between releases, so we read the
    directory listing rather than hard-coding it."""
    listing_url = f"{_ENSEMBL_BASE}/{ensembl_species}/pep/"
    req = urllib.request.Request(listing_url, headers={"User-Agent": "DOMAS-enrichment/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        html = r.read().decode("utf-8", "replace")
    names = re.findall(r'href="([^"]+\.pep\.all\.fa\.gz)"', html)
    if not names:
        raise RuntimeError(f"no .pep.all.fa.gz found in {listing_url}")
    return listing_url + names[0], names[0]


def download_uniprot(data_dir, species, force=False):
    out = _subdir(data_dir, "uniprot")
    meta = SPECIES[species]
    proteome = meta["up"]                       # e.g. UP000005640_9606
    pid = proteome.split("_")[0]                # e.g. UP000005640
    url = f"{_UNIPROT_BASE}/{pid}/{proteome}.dat.gz"
    _download(url, os.path.join(out, f"{proteome}.dat.gz"), force)


def download_afdb(data_dir, species, force=False):
    out = _subdir(data_dir, "afdb")
    name = SPECIES[species]["afdb"] + ".tar"
    _download(f"{_AFDB_BASE}/{name}", os.path.join(out, name), force)


def download_ensembl(data_dir, species, force=False):
    out = _subdir(data_dir, "ensembl")
    url, name = _resolve_ensembl_pep_url(SPECIES[species]["ensembl"])
    _download(url, os.path.join(out, name), force)


def download_pfam(data_dir, force=False):
    out = _subdir(data_dir, "pfam")
    _download(_PFAM_URL, os.path.join(out, "Pfam-A.hmm.gz"), force)


def download_alphamissense(data_dir, force=False):
    # See the _ALPHAMISSENSE_URL definition above for the license (CC BY 4.0).
    out = _subdir(data_dir, "alphamissense")
    _download(_ALPHAMISSENSE_URL, os.path.join(out, "AlphaMissense_aa_substitutions.tsv.gz"), force)


def download_gnomad(data_dir, force=False):
    out = _subdir(data_dir, "gnomad")
    _download(_GNOMAD_URL, os.path.join(out, "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"), force)


def download_all(data_dir, species=None, only=None, force=False):
    """Download every requested source for every requested species.
    `species` limits to a subset of SPECIES (default all); `only` limits to a
    subset of ALL_SOURCES (default all). Pfam is species-independent."""
    species = species or list(SPECIES)
    only = only or list(ALL_SOURCES)
    bad_sp = [s for s in species if s not in SPECIES]
    bad_src = [s for s in only if s not in ALL_SOURCES + OPTIN_SOURCES]
    if bad_sp:
        sys.exit(f"Unknown species: {bad_sp}. Known: {list(SPECIES)}")
    if bad_src:
        sys.exit(f"Unknown source: {bad_src}. Known: {list(ALL_SOURCES)}")

    if "pfam" in only:
        print("[pfam] (species-independent)")
        download_pfam(data_dir, force)
    if "alphamissense" in only:
        print("[alphamissense] (human canonical, ~1.2 GB; CC BY 4.0)")
        download_alphamissense(data_dir, force)
    if "gnomad" in only:
        print("[gnomad] (per-gene LoF constraint, ~5 MB)")
        download_gnomad(data_dir, force)
    for sp in species:
        for src in only:
            if src == "pfam":
                continue
            print(f"[{src}] {sp}")
            {"uniprot": download_uniprot, "afdb": download_afdb,
             "ensembl": download_ensembl}[src](data_dir, sp, force)
    print("download: done")
