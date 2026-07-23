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
CREATE TABLE IF NOT EXISTS afdb_rsa (
    accession   TEXT PRIMARY KEY,
    length      INTEGER,
    mean_rsa    REAL,
    rsa         TEXT      -- comma-joined per-residue relative solvent accessibility (0=buried..1=exposed)
);
CREATE TABLE IF NOT EXISTS gene_constraint (
    accession   TEXT PRIMARY KEY,   -- UniProt accession
    gene        TEXT,
    loeuf       REAL,               -- gnomAD oe_lof_upper (lower = more loss-intolerant)
    pli         REAL
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

# Tien et al. 2013 theoretical maximum solvent-accessible surface area (A^2),
# used to turn absolute SASA into relative solvent accessibility (RSA).
_MAX_ASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225,
            "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,
            "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285,
            "TYR": 263, "VAL": 174}


def _rsa_from_pdb(raw_text):
    """Per-residue relative solvent accessibility (0=buried..1=exposed) from a PDB
    model, via Bio.PDB Shrake-Rupley. Biopython is imported lazily so the rest of
    build.py works without it. Returns [] on any parse/compute failure."""
    from io import StringIO
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    try:
        st = PDBParser(QUIET=True).get_structure("m", StringIO(raw_text))
        res = [r for r in list(st[0].get_chains())[0] if r.id[0] == " "]
        if not res:
            return []
        ShrakeRupley().compute(st, level="R")
        return [round(min(1.0, r.sasa / _MAX_ASA.get(r.resname, 200)), 3) for r in res]
    except Exception:
        return []


def build_afdb_rsa(data_dir, con, species_list):
    """Per-residue relative solvent accessibility from the AlphaFold DB proteome
    tars (same source as build_afdb), stored in afdb_rsa. Opt-in: Shrake-Rupley is
    ~0.1 s/structure, so this is slower than the pLDDT pass and is run explicitly."""
    total = 0
    for sp in species_list:
        dataset = download.SPECIES[sp]["afdb"]
        raw = os.path.join(data_dir, "afdb", dataset + ".tar")
        if not os.path.exists(raw):
            print(f"  [afdb_rsa] {sp}: missing {raw} - run -download afdb first; skipping")
            continue
        n = 0
        with tarfile.open(raw, "r") as tar:
            for member in tar:
                m = _AF_NAME.search(member.name)
                if not m:
                    continue
                acc = m.group(1)
                text = gzip.decompress(tar.extractfile(member).read()).decode("latin-1")
                vals = _rsa_from_pdb(text)
                if not vals:
                    continue
                con.execute(
                    "INSERT OR REPLACE INTO afdb_rsa(accession,length,mean_rsa,rsa) "
                    "VALUES(?,?,?,?)",
                    (acc, len(vals), round(sum(vals) / len(vals), 3),
                     ",".join(str(v) for v in vals)))
                n += 1
                if n % 2000 == 0:
                    con.commit(); print(f"  [afdb_rsa] {sp}: {n}...")
        con.commit()
        _record_meta(con, "afdb_rsa", dataset, sp, f"{download._AFDB_BASE}/{dataset}.tar",
                     raw, n, notes="RSA from Shrake-Rupley on F1 model; Tien 2013 max-ASA")
        total += n
        print(f"  [afdb_rsa] {sp}: {n} proteins")
    return total


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


# -------------------------------------------------------------- pfam_match ----
_PFAM_MATCH_SCHEMA = """
CREATE TABLE IF NOT EXISTS pfam_match (
    protein_ensembl_id TEXT,     -- versionless ENSP (matches ensembl_sequence)
    pfam_acc           TEXT,     -- PFxxxxx (version stripped)
    pfam_name          TEXT,     -- e.g. zf-C2H2, RRM_1, Pkinase
    ali_start          INTEGER,  -- alignment on the protein
    ali_end            INTEGER,
    env_start          INTEGER,  -- envelope (fuller domain extent)
    env_end            INTEGER,
    hmm_from           INTEGER,  -- match on the Pfam model
    hmm_to             INTEGER,
    hmm_len            INTEGER,  -- model length
    bitscore           REAL,     -- domain bitscore (SPADE-style integrity metric)
    ievalue            REAL,     -- domain independent E-value
    coverage           INTEGER   -- round(100*(hmm_to-hmm_from+1)/hmm_len)
);
CREATE INDEX IF NOT EXISTS ix_pfam_match_prot ON pfam_match(protein_ensembl_id);
"""


def _parse_hmmsearch_domtbl(path):
    """Yield pfam_match rows from an hmmsearch --domtblout file.

    hmmsearch searches profiles->sequences, so vs hmmscan the target/query are
    swapped: target = the sequence (our protein), query = the Pfam profile.
    1-based columns used: 1 target(seq) name, 4 query(profile) name, 5 query acc
    (PFxxxxx), 6 qlen (model length), 13 i-Evalue, 14 domain score, 16-17 hmm
    from/to, 18-19 ali from/to, 20-21 env from/to.
    """
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 23:
                continue
            prot = f[0]
            pfam_name = f[3]
            pfam_acc = f[4].split(".")[0]                 # strip version suffix
            hmm_len = int(f[5])
            ievalue = float(f[12])
            bitscore = float(f[13])
            hmm_from, hmm_to = int(f[15]), int(f[16])
            ali_from, ali_to = int(f[17]), int(f[18])
            env_from, env_to = int(f[19]), int(f[20])
            cov = round(100 * (hmm_to - hmm_from + 1) / hmm_len)
            yield (prot, pfam_acc, pfam_name, ali_from, ali_to, env_from, env_to,
                   hmm_from, hmm_to, hmm_len, bitscore, ievalue, cov)


def build_pfam_match(data_dir, con, species=None, chunk_size=20000):
    """Populate the pfam_match cache: every protein in ensembl_sequence scanned
    against Pfam-A, one row per domain hit (see _PFAM_MATCH_SCHEMA).

    Uses hmmsearch (profiles -> sequences): it reads the pressed Pfam profile DB
    once and streams the proteome, so it is markedly faster than per-protein
    hmmscan for a bulk fill and gives identical hits. Incremental - proteins
    already in pfam_match are skipped, so it can resume or top up new isoforms.

    NOTE: this is the expensive step (order ~0.15 CPU-s/protein). It is NOT part
    of the default build_all/only set; run it explicitly (only=['pfam_match']).
    """
    hmm = os.path.join(data_dir, "pfam", "Pfam-A.hmm")
    if not os.path.exists(hmm + ".h3i"):
        print(f"  [pfam_match] {hmm} not pressed - run build_pfam first; skipping")
        return 0
    con.executescript(_PFAM_MATCH_SCHEMA)
    done = {r[0] for r in con.execute("SELECT DISTINCT protein_ensembl_id FROM pfam_match")}

    q = "SELECT protein_ensembl_id, seq FROM ensembl_sequence WHERE seq IS NOT NULL AND seq!=''"
    params = ()
    if species:
        q += " AND species IN (%s)" % ",".join("?" * len(species))
        params = tuple(species)
    todo = [(pe, s) for pe, s in con.execute(q, params) if pe not in done]
    print(f"  [pfam_match] {len(todo):,} proteins to scan ({len(done):,} already present)")

    insert = ("INSERT INTO pfam_match(protein_ensembl_id,pfam_acc,pfam_name,ali_start,ali_end,"
              "env_start,env_end,hmm_from,hmm_to,hmm_len,bitscore,ievalue,coverage) "
              "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)")
    total = 0
    for i in range(0, len(todo), chunk_size):
        chunk = todo[i:i + chunk_size]
        fa = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
        for pe, s in chunk:
            fa.write(f">{pe}\n{s}\n")
        fa.close()
        out = fa.name + ".domtbl"
        subprocess.run(["hmmsearch", "--cut_ga", "--domtblout", out, hmm, fa.name],
                       check=True, capture_output=True)
        con.executemany(insert, _parse_hmmsearch_domtbl(out))
        con.commit()
        os.remove(fa.name)
        os.remove(out)
        total = con.execute("SELECT COUNT(*) FROM pfam_match").fetchone()[0]
        print(f"  [pfam_match] {min(i + chunk_size, len(todo)):,}/{len(todo):,} proteins, {total:,} matches")
    _record_meta(con, "pfam", "pfam_match", None, download._PFAM_URL, hmm, total,
                 notes="per-protein Pfam-A matches (hmmsearch --cut_ga); cache for HMMER-free lookup")
    print(f"  [pfam_match] done: {total:,} matches")
    return total


# ---------------------------------------------------------- AlphaMissense ----
_AM_SCHEMA = """
CREATE TABLE IF NOT EXISTS am_pathogenicity (
    accession TEXT PRIMARY KEY,   -- UniProt canonical accession
    length    INTEGER,
    mean_am   REAL,               -- whole-protein mean pathogenicity
    am        TEXT                -- comma-joined per-residue mean pathogenicity
);
"""
_AM_VAR = re.compile(r'^[A-Z](\d+)[A-Z]$')   # e.g. V2L -> position 2


def build_alphamissense(data_dir, con):
    """Aggregate AlphaMissense per-variant predictions to a per-residue mean
    pathogenicity table (am_pathogenicity), mirroring afdb_plddt. The raw file
    lists every substitution; we average the ~19 substitutions per position.
    Streaming, grouped by protein -> flat memory. See download.py for license
    (predictions are CC BY 4.0). Explicit-only build step (not in default set).
    """
    raw = os.path.join(data_dir, "alphamissense", "AlphaMissense_aa_substitutions.tsv.gz")
    if not os.path.exists(raw):
        print(f"  [alphamissense] missing {raw} - run -download first; skipping")
        return 0
    con.executescript(_AM_SCHEMA)
    insert = "INSERT OR REPLACE INTO am_pathogenicity(accession,length,mean_am,am) VALUES(?,?,?,?)"
    cur_acc, pos_scores, n = None, {}, 0

    def flush():
        nonlocal n
        if not cur_acc or not pos_scores:
            return
        maxp = max(pos_scores)
        arr, allv = [], []
        for p in range(1, maxp + 1):
            v = pos_scores.get(p)
            arr.append(str(round(sum(v) / len(v), 2)) if v else "")
            if v:
                allv.extend(v)
        mean = round(sum(allv) / len(allv), 3) if allv else None
        con.execute(insert, (cur_acc, maxp, mean, ",".join(arr)))
        n += 1

    with gzip.open(raw, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("uniprot"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            acc, var, score = f[0], f[1], f[2]
            m = _AM_VAR.match(var)
            if not m:
                continue
            if acc != cur_acc:
                flush()
                pos_scores, cur_acc = {}, acc
                if n and n % 2000 == 0:
                    con.commit()
            pos_scores.setdefault(int(m.group(1)), []).append(float(score))
    flush()
    con.commit()
    _record_meta(con, "alphamissense", "AlphaMissense_aa_substitutions", None,
                 download._ALPHAMISSENSE_URL, raw, n,
                 notes="per-residue mean pathogenicity; predictions CC BY 4.0 (see download.py)")
    print(f"  [alphamissense] {n:,} proteins")
    return n


# ------------------------------------------------------------------- driver --
_GN_RE = re.compile(r'Name=([^;{ ]+)')


def _acc_gene_from_dat(path):
    """UniProt accession -> gene symbol from a reference-proteome .dat(.gz)."""
    out = {}; acc = gene = None
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("AC") and acc is None:
                acc = ln.split()[1].rstrip(";")
            elif ln.startswith("GN") and gene is None:
                m = _GN_RE.search(ln)
                if m:
                    gene = m.group(1)
            elif ln.startswith("//"):
                if acc and gene:
                    out[acc] = gene
                acc = gene = None
    return out


def build_gnomad(data_dir, con, species_list):
    """Per-gene gnomAD LoF constraint (LOEUF, pLI) keyed by UniProt accession
    (gene symbol -> accession via the human reference proteome). Opt-in; human
    only (gnomAD is human)."""
    raw = os.path.join(data_dir, "gnomad", "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
    if not os.path.exists(raw):
        print(f"  [gnomad] missing {raw} - run -download gnomad first; skipping")
        return 0
    g2l, g2p = {}, {}
    with gzip.open(raw, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        gi, li, pi = h.index("gene"), h.index("oe_lof_upper"), h.index("pLI")
        for ln in fh:
            c = ln.rstrip("\n").split("\t")
            try:
                if c[li] not in ("", "NA"):
                    g2l[c[gi]] = float(c[li])
                if c[pi] not in ("", "NA"):
                    g2p[c[gi]] = float(c[pi])
            except (ValueError, IndexError):
                pass
    dat = os.path.join(data_dir, "uniprot", download.SPECIES["H_sapiens"]["up"] + ".dat.gz")
    acc2gene = _acc_gene_from_dat(dat) if os.path.exists(dat) else {}
    n = 0
    for acc, gene in acc2gene.items():
        if gene in g2l or gene in g2p:
            con.execute("INSERT OR REPLACE INTO gene_constraint(accession,gene,loeuf,pli) "
                        "VALUES(?,?,?,?)", (acc, gene, g2l.get(gene), g2p.get(gene)))
            n += 1
    con.commit()
    _record_meta(con, "gnomad", "v2.1.1", "H_sapiens", download._GNOMAD_URL, raw, n,
                 notes="LOEUF=oe_lof_upper, pLI; gene->accession via reference proteome .dat")
    print(f"  [gnomad] {n} accessions with constraint")
    return n


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
    # afdb_rsa is the opt-in burial pass (Shrake-Rupley, slower); reuses the afdb tar.
    if "afdb_rsa" in only:
        build_afdb_rsa(data_dir, con, species)
    if "ensembl" in only:
        build_ensembl(data_dir, con, species)
    if "pfam" in only:
        build_pfam(data_dir, con)
    # pfam_match is the expensive per-protein HMMER scan; only when asked for
    # explicitly (never in the default source set).
    if "pfam_match" in only:
        build_pfam_match(data_dir, con, species=species)
    if "alphamissense" in only:
        build_alphamissense(data_dir, con)
    if "gnomad" in only:
        build_gnomad(data_dir, con, species)
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
