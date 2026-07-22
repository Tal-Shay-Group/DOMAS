"""
Build enrichment.sqlite from the raw files fetched by download.py.

Parsing only - no network. Run download.py first (domas.py -download), then
this (domas.py -build). Tables use source-prefixed, self-describing names:

    afdb_plddt        accession, length, mean_plddt, plddt
    uniprot_feature   acc, ftype, start, end, note      (+ index on acc)
    uniprot_alias     acc, primacc                       (+ index on acc)
    ensembl_sequence  protein_ensembl_id, transcript_ensembl_id, species, seq
    source_meta       source, dataset, species, version, source_url,
                      license, downloaded_at, built_at, row_count, notes

Pfam is not a table: hmmscan reads Pfam-A.hmm from disk, so build only runs
hmmpress on it and records a source_meta row for provenance.
"""
import datetime
import gzip
import os
import re
import sqlite3
import subprocess
import tarfile

import download  # SPECIES, LICENSES, ALL_SOURCES, _UNIPROT_BASE etc.

_SCHEMA = """
CREATE TABLE IF NOT EXISTS afdb_plddt (
    accession   TEXT PRIMARY KEY,
    length      INTEGER,
    mean_plddt  REAL,
    plddt       TEXT      -- comma-joined per-residue pLDDT
);
CREATE TABLE IF NOT EXISTS uniprot_feature (
    acc   TEXT,
    ftype TEXT,
    start INTEGER,
    end   INTEGER,
    note  TEXT
);
CREATE INDEX IF NOT EXISTS ix_uniprot_feature_acc ON uniprot_feature(acc);
CREATE TABLE IF NOT EXISTS uniprot_alias (
    acc     TEXT,   -- secondary accession
    primacc TEXT    -- current primary accession
);
CREATE INDEX IF NOT EXISTS ix_uniprot_alias_acc ON uniprot_alias(acc);
CREATE TABLE IF NOT EXISTS ensembl_sequence (
    protein_ensembl_id    TEXT PRIMARY KEY,   -- versionless ENSP
    transcript_ensembl_id TEXT,               -- versionless ENST
    species               TEXT,
    seq                   TEXT
);
CREATE INDEX IF NOT EXISTS ix_ensembl_sequence_tx ON ensembl_sequence(transcript_ensembl_id);
CREATE TABLE IF NOT EXISTS source_meta (
    source        TEXT,
    dataset       TEXT,
    species       TEXT,
    version       TEXT,
    source_url    TEXT,
    license       TEXT,
    downloaded_at TEXT,
    built_at      TEXT,
    row_count     INTEGER,
    notes         TEXT
);
"""


def _now():
    return datetime.datetime.now().replace(microsecond=0).isoformat()


def _file_mtime(path):
    if not os.path.exists(path):
        return None
    return datetime.datetime.fromtimestamp(os.path.getmtime(path)).replace(microsecond=0).isoformat()


def _record_meta(con, source, dataset, species, source_url, raw_path, row_count, notes=""):
    con.execute(
        "DELETE FROM source_meta WHERE source=? AND dataset=? AND IFNULL(species,'')=IFNULL(?,'')",
        (source, dataset, species))
    con.execute(
        "INSERT INTO source_meta(source,dataset,species,version,source_url,license,"
        "downloaded_at,built_at,row_count,notes) VALUES(?,?,?,?,?,?,?,?,?,?)",
        (source, dataset, species, "current_release", source_url,
         download.LICENSES.get(source, ""), _file_mtime(raw_path), _now(), row_count, notes))
    con.commit()


# ---------------------------------------------------------------- UniProt ----
_LOC_RE = re.compile(r"[<>?]")


def _parse_loc(tok):
    """'14..74' -> (14,74); '2' -> (2,2); '<1..>300' -> (1,300). Returns
    (None,None) if either coord is missing or implausibly large (>7 digits,
    which in SwissProt flat files means a wrapped /note leaked into the slot)."""
    tok = _LOC_RE.sub("", tok).strip()
    parts = tok.split("..")
    try:
        start = int(parts[0])
        end = int(parts[-1])
    except (ValueError, IndexError):
        return None, None
    if len(str(start)) > 7 or len(str(end)) > 7:
        return None, None
    return start, end


def _iter_uniprot_entries(path):
    """Yield (primary_acc, [secondary_accs], [(ftype,start,end,note)]) per entry."""
    accs, feats = [], []
    cur = None            # (ftype, start, end, [note_fragments])
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            tag = line[:2]
            if tag == "AC":
                for a in line[5:].strip().rstrip(";").split(";"):
                    a = a.strip()
                    if a:
                        accs.append(a)
            elif tag == "FT":
                if line[5:6] != " ":                     # header: "FT   TYPE  loc"
                    if cur:
                        feats.append(cur)
                    body = line[5:].rstrip("\n")
                    bits = body.split(None, 1)
                    ftype = bits[0]
                    start, end = _parse_loc(bits[1]) if len(bits) > 1 else (None, None)
                    cur = [ftype, start, end, []]
                elif cur is not None:                    # continuation (/note=, /evidence= ...)
                    cur[3].append(line[5:].strip())
            elif tag == "//":
                if cur:
                    feats.append(cur)
                    cur = None
                if accs:
                    prim = accs[0]
                    parsed = []
                    for ftype, start, end, frags in feats:
                        m = re.search(r'/note="([^"]*)"', " ".join(frags))
                        parsed.append((ftype, start, end, m.group(1) if m else ""))
                    yield prim, accs[1:], parsed
                accs, feats, cur = [], [], None


def build_uniprot(data_dir, con, species_list):
    total_feat = total_alias = 0
    for sp in species_list:
        proteome = download.SPECIES[sp]["up"]
        raw = os.path.join(data_dir, "uniprot", f"{proteome}.dat.gz")
        if not os.path.exists(raw):
            print(f"  [uniprot] {sp}: missing {raw} - run -download first; skipping")
            continue
        nfeat = nalias = 0
        for prim, secs, feats in _iter_uniprot_entries(raw):
            for s in secs:
                con.execute("INSERT INTO uniprot_alias(acc,primacc) VALUES(?,?)", (s, prim))
                nalias += 1
            for ftype, start, end, note in feats:
                con.execute(
                    "INSERT INTO uniprot_feature(acc,ftype,start,end,note) VALUES(?,?,?,?,?)",
                    (prim, ftype, start, end, note))
                nfeat += 1
        con.commit()
        pid = proteome.split("_")[0]
        url = f"{download._UNIPROT_BASE}/{pid}/{proteome}.dat.gz"
        _record_meta(con, "uniprot", proteome, sp, url, raw, nfeat,
                     notes=f"{nalias} aliases; features parsed from SwissProt FT lines")
        total_feat += nfeat
        total_alias += nalias
        print(f"  [uniprot] {sp}: {nfeat} features, {nalias} aliases")
    return total_feat, total_alias


# ------------------------------------------------------------------- AFDB ----
def _plddt_from_pdb(raw_bytes):
    vals = []
    for line in raw_bytes.decode("latin-1").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                vals.append(round(float(line[60:66]), 1))
            except ValueError:
                pass
    return vals


_AF_NAME = re.compile(r"AF-([A-Z0-9]+)-F1-model_v\d+\.pdb\.gz$")


def build_afdb(data_dir, con, species_list):
    total = 0
    for sp in species_list:
        dataset = download.SPECIES[sp]["afdb"]
        raw = os.path.join(data_dir, "afdb", dataset + ".tar")
        if not os.path.exists(raw):
            print(f"  [afdb] {sp}: missing {raw} - run -download first; skipping")
            continue
        n = 0
        with tarfile.open(raw, "r") as tar:
            for member in tar:
                m = _AF_NAME.search(member.name)
                if not m:
                    continue          # only fragment 1 (covers proteins <=2700 aa)
                acc = m.group(1)
                data = tar.extractfile(member).read()
                vals = _plddt_from_pdb(gzip.decompress(data))
                if not vals:
                    continue
                con.execute(
                    "INSERT OR REPLACE INTO afdb_plddt(accession,length,mean_plddt,plddt) "
                    "VALUES(?,?,?,?)",
                    (acc, len(vals), round(sum(vals) / len(vals), 2),
                     ",".join(str(v) for v in vals)))
                n += 1
                if n % 5000 == 0:
                    con.commit()
        con.commit()
        url = f"{download._AFDB_BASE}/{dataset}.tar"
        _record_meta(con, "afdb", dataset, sp, url, raw, n,
                     notes="pLDDT from CA B-factors; fragment F1 only")
        total += n
        print(f"  [afdb] {sp}: {n} proteins")
    return total


# ---------------------------------------------------------------- Ensembl ----
def _open_maybe_gz(path):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path, errors="replace")


def build_ensembl(data_dir, con, species_list):
    total = 0
    ens_dir = os.path.join(data_dir, "ensembl")
    for sp in species_list:
        prefix = download.SPECIES[sp]["ensembl"]        # e.g. homo_sapiens
        # filename is <Species>.<assembly>.pep.all.fa[.gz]; match on the leading token
        want = prefix.split("_")[0].capitalize()        # 'Homo'
        cands = [f for f in os.listdir(ens_dir)
                 if f.lower().startswith(prefix.split("_")[0]) and ".pep.all.fa" in f] if os.path.isdir(ens_dir) else []
        if not cands:
            print(f"  [ensembl] {sp}: no *.pep.all.fa* for {prefix} - run -download first; skipping")
            continue
        raw = os.path.join(ens_dir, cands[0])
        n = 0
        pid = tx = None
        seq = []

        def _flush():
            nonlocal n
            if pid and seq:
                con.execute(
                    "INSERT OR REPLACE INTO ensembl_sequence"
                    "(protein_ensembl_id,transcript_ensembl_id,species,seq) VALUES(?,?,?,?)",
                    (pid, tx, sp, "".join(seq)))
                n += 1

        with _open_maybe_gz(raw) as fh:
            for line in fh:
                if line.startswith(">"):
                    _flush()
                    seq = []
                    pid = line[1:].split()[0].split(".")[0]
                    mtx = re.search(r"transcript:(\S+)", line)
                    tx = mtx.group(1).split(".")[0] if mtx else None
                else:
                    seq.append(line.strip())
            _flush()
        con.commit()
        _record_meta(con, "ensembl", cands[0], sp,
                     f"{download._ENSEMBL_BASE}/{prefix}/pep/{cands[0]}", raw, n,
                     notes="peptide FASTA indexed by versionless ENSP/ENST")
        total += n
        print(f"  [ensembl] {sp}: {n} proteins")
    return total


# ------------------------------------------------------------------- Pfam ----
def build_pfam(data_dir, con):
    pfam_dir = os.path.join(data_dir, "pfam")
    gz = os.path.join(pfam_dir, "Pfam-A.hmm.gz")
    hmm = os.path.join(pfam_dir, "Pfam-A.hmm")
    if not os.path.exists(hmm) and os.path.exists(gz):
        print("  [pfam] gunzip Pfam-A.hmm.gz")
        with gzip.open(gz, "rb") as fi, open(hmm, "wb") as fo:
            while True:
                b = fi.read(1 << 20)
                if not b:
                    break
                fo.write(b)
    if not os.path.exists(hmm):
        print(f"  [pfam] missing {hmm} - run -download first; skipping")
        return 0
    if not os.path.exists(hmm + ".h3i"):
        print("  [pfam] hmmpress")
        subprocess.run(["hmmpress", hmm], check=True, capture_output=True)
    # count families for the meta row
    n = 0
    with open(hmm, errors="replace") as fh:
        for line in fh:
            if line.startswith("NAME"):
                n += 1
    _record_meta(con, "pfam", "Pfam-A.hmm", None, download._PFAM_URL, hmm, n,
                 notes="on-disk HMM library (hmmpress'd); not stored as a table")
    print(f"  [pfam] {n} families pressed")
    return n


# ------------------------------------------------------------------- driver --
def build_all(data_dir, db_path, species=None, only=None, delete_raw=False):
    species = species or list(download.SPECIES)
    only = only or list(download.ALL_SOURCES)
    os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)
    con = sqlite3.connect(db_path)
    con.executescript(_SCHEMA)
    print(f"building {db_path}")
    if "uniprot" in only:
        build_uniprot(data_dir, con, species)
    if "afdb" in only:
        build_afdb(data_dir, con, species)
    if "ensembl" in only:
        build_ensembl(data_dir, con, species)
    if "pfam" in only:
        build_pfam(data_dir, con)
    con.execute("ANALYZE")
    con.commit()
    con.close()
    if delete_raw:
        _delete_raw(data_dir, species, only)
    print("build: done")


def _delete_raw(data_dir, species, only):
    """Remove the large raw downloads once they are in the DB. Pfam-A.hmm is
    kept (hmmscan needs it); only its .gz is removed."""
    victims = []
    for sp in species:
        if "uniprot" in only:
            victims.append(os.path.join(data_dir, "uniprot", download.SPECIES[sp]["up"] + ".dat.gz"))
        if "afdb" in only:
            victims.append(os.path.join(data_dir, "afdb", download.SPECIES[sp]["afdb"] + ".tar"))
    if "pfam" in only:
        victims.append(os.path.join(data_dir, "pfam", "Pfam-A.hmm.gz"))
    for p in victims:
        if os.path.exists(p):
            os.remove(p)
            print(f"  removed raw: {p}")
