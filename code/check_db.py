"""
DB sanity-check script for the DoChaP merged SQLite database.

Canonical field values in Transcripts table:
    0 = not canonical
    1 = canonical from RefSeq only
    2 = canonical from Ensembl only
    3 = canonical from both RefSeq and Ensembl

Two kinds of checks:

1. `run_anomaly_checks(con)` - purely local SQL-based consistency checks
   (no network access). Returns a dict of {check_name: DataFrame of
   offending rows}.

2. `compare_random_genes_with_external_dbs(con, n)` - picks `n` random
   genes and compares their transcripts/exons/proteins/domains against
   Ensembl and NCBI/RefSeq via their public REST APIs. Both APIs are
   rate-limited, so requests are throttled and retried on 429/5xx.

Run as a script:
    python check_db.py anomalies --db path/to/DB.sqlite --out-dir reports/
    python check_db.py compare --db path/to/DB.sqlite -n 20 --out report.csv
"""

import argparse
import json
import logging
import random
import sqlite3
import time
import urllib.error
import urllib.request

import pandas as pd

logger = logging.getLogger(__name__)

SPECIE_TO_ENSEMBL = {
    'H_sapiens': 'homo_sapiens',
    'M_musculus': 'mus_musculus',
    'R_norvegicus': 'rattus_norvegicus',
}


# ---------------------------------------------------------------------------
# Local, SQL-based anomaly checks
# ---------------------------------------------------------------------------

def check_multiple_canonical_transcripts(con):
    """Genes with more than one transcript flagged canonical.

    Canonical field: 0=not canonical, 1=RefSeq only, 2=Ensembl only, 3=both.
    canonical != 0 detects any canonical source.
    """
    query = """
        SELECT gene_GeneID_id, gene_ensembl_id, COUNT(*) AS canonical_count
        FROM Transcripts
        WHERE canonical != 0
        GROUP BY gene_GeneID_id, gene_ensembl_id
        HAVING COUNT(*) > 1
    """
    return pd.read_sql_query(query, con)


def check_genes_without_canonical_transcript(con):
    """Genes that have transcripts but none flagged canonical.

    Canonical field: 0=not canonical, 1=RefSeq only, 2=Ensembl only, 3=both.
    Detects genes where all transcripts have canonical=0 (no canonical source).
    """
    query = """
        SELECT g.gene_GeneID_id, g.gene_ensembl_id, g.gene_symbol, g.specie
        FROM Genes g
        JOIN Transcripts t
          ON t.gene_GeneID_id = g.gene_GeneID_id
        GROUP BY g.gene_GeneID_id, g.gene_ensembl_id
        HAVING SUM(CASE WHEN t.canonical != 0 THEN 1 ELSE 0 END) = 0
    """
    return pd.read_sql_query(query, con)


def check_genes_without_transcripts(con):
    """Genes with no rows at all in Transcripts."""
    query = """
        SELECT g.gene_GeneID_id, g.gene_ensembl_id, g.gene_symbol, g.specie
        FROM Genes g
        LEFT JOIN Transcripts t ON t.gene_GeneID_id = g.gene_GeneID_id
        WHERE t.gene_GeneID_id IS NULL
    """
    return pd.read_sql_query(query, con)


def check_missing_transcript_ids(con):
    """Transcripts with neither a transcript_ensembl_id nor a transcript_refseq_id."""
    query = """
        SELECT *
        FROM Transcripts
        WHERE (transcript_ensembl_id IS NULL OR transcript_ensembl_id = '')
          AND (transcript_refseq_id IS NULL OR transcript_refseq_id = '')
    """
    return pd.read_sql_query(query, con)


def check_transcripts_without_exons(con):
    """
    Transcripts with no rows in Transcript_Exon (matched on whichever of
    transcript_ensembl_id / transcript_refseq_id is set).
    """
    query = """
        SELECT t.gene_GeneID_id, t.gene_ensembl_id,
               t.transcript_ensembl_id, t.transcript_refseq_id, t.canonical
        FROM Transcripts t
        LEFT JOIN Transcript_Exon te
          ON (t.transcript_ensembl_id IS NOT NULL AND te.transcript_ensembl_id = t.transcript_ensembl_id)
          OR (t.transcript_ensembl_id IS NULL AND te.transcript_refseq_id = t.transcript_refseq_id)
        WHERE te.rowid IS NULL
          AND NOT (
              (t.transcript_ensembl_id IS NULL OR t.transcript_ensembl_id = '')
              AND (t.transcript_refseq_id IS NULL OR t.transcript_refseq_id = '')
          )
    """
    return pd.read_sql_query(query, con)


def check_exon_count_mismatch(con):
    """
    Transcripts.exon_count vs the actual number of Transcript_Exon rows
    for that transcript.
    """
    query = """
        SELECT t.gene_ensembl_id, t.transcript_ensembl_id, t.transcript_refseq_id,
               t.exon_count AS declared_exon_count,
               COUNT(te.rowid) AS actual_exon_count
        FROM Transcripts t
        LEFT JOIN Transcript_Exon te
          ON (t.transcript_ensembl_id IS NOT NULL AND te.transcript_ensembl_id = t.transcript_ensembl_id)
          OR (t.transcript_ensembl_id IS NULL AND te.transcript_refseq_id = t.transcript_refseq_id)
        WHERE t.exon_count IS NOT NULL
        GROUP BY t.transcript_ensembl_id, t.transcript_refseq_id
        HAVING t.exon_count != COUNT(te.rowid)
    """
    return pd.read_sql_query(query, con)


def check_invalid_exon_coordinates(con):
    """Transcript_Exon rows where genomic start > end, or CDS start > end."""
    query = """
        SELECT *
        FROM Transcript_Exon
        WHERE genomic_start_tx > genomic_end_tx
           OR (abs_start_CDS IS NOT NULL AND abs_end_CDS IS NOT NULL
               AND abs_start_CDS > abs_end_CDS)
    """
    return pd.read_sql_query(query, con)


def check_duplicate_exon_order(con):
    """Duplicate (transcript, order_in_transcript) pairs in Transcript_Exon."""
    query = """
        SELECT transcript_ensembl_id, transcript_refseq_id, order_in_transcript,
               COUNT(*) AS n
        FROM Transcript_Exon
        GROUP BY transcript_ensembl_id, transcript_refseq_id, order_in_transcript
        HAVING COUNT(*) > 1
    """
    return pd.read_sql_query(query, con)


def check_orphan_transcript_exons(con):
    """Transcript_Exon rows whose transcript id doesn't exist in Transcripts."""
    query = """
        SELECT te.transcript_ensembl_id, te.transcript_refseq_id, COUNT(*) AS n_exons
        FROM Transcript_Exon te
        LEFT JOIN Transcripts t
          ON (te.transcript_ensembl_id IS NOT NULL AND t.transcript_ensembl_id = te.transcript_ensembl_id)
          OR (te.transcript_ensembl_id IS NULL AND t.transcript_refseq_id = te.transcript_refseq_id)
        WHERE t.rowid IS NULL
        GROUP BY te.transcript_ensembl_id, te.transcript_refseq_id
    """
    return pd.read_sql_query(query, con)


def check_canonical_without_protein(con):
    """
    Canonical transcripts that declare a protein id but have no matching row
    in Proteins.

    Canonical field: 0=not canonical, 1=RefSeq only, 2=Ensembl only, 3=both.
    canonical != 0 selects any canonical transcript.
    """
    query = """
        SELECT t.gene_ensembl_id, t.transcript_ensembl_id, t.transcript_refseq_id,
               t.protein_ensembl_id, t.protein_refseq_id
        FROM Transcripts t
        LEFT JOIN Proteins p
          ON (t.protein_ensembl_id IS NOT NULL AND p.protein_ensembl_id = t.protein_ensembl_id)
          OR (t.protein_ensembl_id IS NULL AND p.protein_refseq_id = t.protein_refseq_id)
        WHERE t.canonical != 0
          AND NOT (
              (t.protein_ensembl_id IS NULL OR t.protein_ensembl_id = '')
              AND (t.protein_refseq_id IS NULL OR t.protein_refseq_id = '')
          )
          AND p.rowid IS NULL
    """
    return pd.read_sql_query(query, con)


def check_domain_event_orphans(con):
    """DomainEvent rows whose protein id is missing from Proteins, or whose type_id is missing from DomainType."""
    query_protein = """
        SELECT de.protein_ensembl_id, de.protein_refseq_id, COUNT(*) AS n
        FROM DomainEvent de
        LEFT JOIN Proteins p
          ON (de.protein_ensembl_id IS NOT NULL AND p.protein_ensembl_id = de.protein_ensembl_id)
          OR (de.protein_ensembl_id IS NULL AND p.protein_refseq_id = de.protein_refseq_id)
        WHERE p.rowid IS NULL
        GROUP BY de.protein_ensembl_id, de.protein_refseq_id
    """
    query_type = """
        SELECT de.type_id, COUNT(*) AS n
        FROM DomainEvent de
        LEFT JOIN DomainType dt ON dt.type_id = de.type_id
        WHERE dt.type_id IS NULL
        GROUP BY de.type_id
    """
    return {
        'domain_event_unknown_protein': pd.read_sql_query(query_protein, con),
        'domain_event_unknown_type': pd.read_sql_query(query_type, con),
    }


def check_invalid_domain_ranges(con):
    """DomainEvent rows where AA_start > AA_end (or either is negative)."""
    query = """
        SELECT *
        FROM DomainEvent
        WHERE AA_start > AA_end OR AA_start < 0 OR AA_end < 0
    """
    return pd.read_sql_query(query, con)


def check_duplicate_genes(con):
    """Same gene_ensembl_id appearing more than once within the same species."""
    query = """
        SELECT gene_ensembl_id, specie, COUNT(*) AS n
        FROM Genes
        WHERE gene_ensembl_id IS NOT NULL AND gene_ensembl_id != ''
        GROUP BY gene_ensembl_id, specie
        HAVING COUNT(*) > 1
    """
    return pd.read_sql_query(query, con)


def check_orphan_transcripts(con):
    """Transcripts that don't match any Genes row on either gene_GeneID_id or gene_ensembl_id."""
    query = """
        SELECT t.gene_GeneID_id, t.gene_ensembl_id, t.transcript_ensembl_id, t.transcript_refseq_id
        FROM Transcripts t
        LEFT JOIN Genes g
          ON g.gene_GeneID_id = t.gene_GeneID_id OR g.gene_ensembl_id = t.gene_ensembl_id
        WHERE g.gene_GeneID_id IS NULL AND g.gene_ensembl_id IS NULL
    """
    return pd.read_sql_query(query, con)


def check_transcript_gene_id_mismatch(con):
    """
    Transcripts whose gene_GeneID_id doesn't match the gene_GeneID_id of the
    Genes row sharing the same gene_ensembl_id - e.g. the gene was
    re-annotated under a different NCBI GeneID after the Transcripts table
    was built.
    """
    query = """
        SELECT DISTINCT t.gene_GeneID_id AS transcript_gene_GeneID_id,
               g.gene_GeneID_id AS genes_gene_GeneID_id,
               t.gene_ensembl_id, g.gene_symbol, g.specie
        FROM Transcripts t
        JOIN Genes g ON g.gene_ensembl_id = t.gene_ensembl_id
        WHERE t.gene_GeneID_id != g.gene_GeneID_id
    """
    return pd.read_sql_query(query, con)


def check_invalid_cds_bounds(con):
    """Transcripts where cds_start > cds_end, or the CDS falls outside tx_start/tx_end."""
    query = """
        SELECT gene_ensembl_id, transcript_ensembl_id, transcript_refseq_id,
               tx_start, tx_end, cds_start, cds_end
        FROM Transcripts
        WHERE cds_start IS NOT NULL AND cds_end IS NOT NULL
          AND (cds_start > cds_end OR cds_start < tx_start OR cds_end > tx_end)
    """
    return pd.read_sql_query(query, con)


def check_protein_length_vs_cds(con):
    """
    For canonical transcripts, compare Proteins.length against the CDS
    length implied by Transcript_Exon.abs_end_CDS (max abs_end_CDS should be
    ~ protein_length * 3 + 3 for the stop codon).

    Canonical field: 0=not canonical, 1=RefSeq only, 2=Ensembl only, 3=both.
    canonical != 0 selects any canonical transcript.
    """
    query = """
        SELECT t.gene_ensembl_id, t.transcript_ensembl_id, t.transcript_refseq_id,
               MAX(p.length) AS protein_length,
               MAX(te.abs_end_CDS) AS cds_nt_length
        FROM Transcripts t
        JOIN Proteins p
          ON (t.protein_ensembl_id IS NOT NULL AND p.protein_ensembl_id = t.protein_ensembl_id)
          OR (t.protein_ensembl_id IS NULL AND p.protein_refseq_id = t.protein_refseq_id)
        JOIN Transcript_Exon te
          ON (t.transcript_ensembl_id IS NOT NULL AND te.transcript_ensembl_id = t.transcript_ensembl_id)
          OR (t.transcript_ensembl_id IS NULL AND te.transcript_refseq_id = t.transcript_refseq_id)
        WHERE t.canonical != 0 AND p.length IS NOT NULL AND te.abs_end_CDS IS NOT NULL
        GROUP BY t.transcript_ensembl_id, t.transcript_refseq_id
        HAVING ABS(MAX(te.abs_end_CDS) - (MAX(p.length) * 3 + 3)) > 3
    """
    return pd.read_sql_query(query, con)


def check_orphan_orthology_genes(con):
    """Orthology rows whose A_ensembl_id or B_ensembl_id can't be matched to a Genes row of the right species."""
    query = """
        SELECT o.A_ensembl_id, o.A_species, o.B_ensembl_id, o.B_species
        FROM Orthology o
        LEFT JOIN Genes ga ON ga.gene_ensembl_id = o.A_ensembl_id AND ga.specie = o.A_species
        LEFT JOIN Genes gb ON gb.gene_ensembl_id = o.B_ensembl_id AND gb.specie = o.B_species
        WHERE ga.gene_ensembl_id IS NULL OR gb.gene_ensembl_id IS NULL
    """
    return pd.read_sql_query(query, con)


def check_spliceindomains_orphans(con):
    """SpliceInDomains rows referencing an unknown transcript or an unknown type_id."""
    query_transcript = """
        SELECT s.transcript_ensembl_id, s.transcript_refseq_id, COUNT(*) AS n
        FROM SpliceInDomains s
        LEFT JOIN Transcripts t
          ON (s.transcript_ensembl_id IS NOT NULL AND t.transcript_ensembl_id = s.transcript_ensembl_id)
          OR (s.transcript_ensembl_id IS NULL AND t.transcript_refseq_id = s.transcript_refseq_id)
        WHERE t.rowid IS NULL
        GROUP BY s.transcript_ensembl_id, s.transcript_refseq_id
    """
    query_type = """
        SELECT s.type_id, COUNT(*) AS n
        FROM SpliceInDomains s
        LEFT JOIN DomainType dt ON dt.type_id = s.type_id
        WHERE dt.type_id IS NULL
        GROUP BY s.type_id
    """
    return {
        'spliceindomains_unknown_transcript': pd.read_sql_query(query_transcript, con),
        'spliceindomains_unknown_type': pd.read_sql_query(query_type, con),
    }


def check_proteins_without_domains(con):
    """
    Proteins (referenced by a canonical transcript) with no DomainEvent rows
    at all. Informational - a protein genuinely may have no annotated
    domains, but a large fraction of canonical proteins with zero domains
    can indicate a loading problem.

    Canonical field: 0=not canonical, 1=RefSeq only, 2=Ensembl only, 3=both.
    canonical != 0 selects any canonical transcript.
    """
    query = """
        SELECT p.protein_ensembl_id, p.protein_refseq_id,
               p.transcript_ensembl_id, p.transcript_refseq_id
        FROM Proteins p
        JOIN Transcripts t
          ON (p.transcript_ensembl_id IS NOT NULL AND t.transcript_ensembl_id = p.transcript_ensembl_id)
          OR (p.transcript_ensembl_id IS NULL AND t.transcript_refseq_id = p.transcript_refseq_id)
        LEFT JOIN DomainEvent de
          ON (p.protein_ensembl_id IS NOT NULL AND de.protein_ensembl_id = p.protein_ensembl_id)
          OR (p.protein_ensembl_id IS NULL AND de.protein_refseq_id = p.protein_refseq_id)
        WHERE t.canonical != 0
          AND de.rowid IS NULL
    """
    return pd.read_sql_query(query, con)


ANOMALY_CHECKS = {
    'multiple_canonical_transcripts': check_multiple_canonical_transcripts,
    'genes_without_canonical_transcript': check_genes_without_canonical_transcript,
    'genes_without_transcripts': check_genes_without_transcripts,
    'missing_transcript_ids': check_missing_transcript_ids,
    'transcripts_without_exons': check_transcripts_without_exons,
    'exon_count_mismatch': check_exon_count_mismatch,
    'invalid_exon_coordinates': check_invalid_exon_coordinates,
    'duplicate_exon_order': check_duplicate_exon_order,
    'orphan_transcript_exons': check_orphan_transcript_exons,
    'canonical_without_protein': check_canonical_without_protein,
    'domain_event_orphans': check_domain_event_orphans,
    'invalid_domain_ranges': check_invalid_domain_ranges,
    'duplicate_genes': check_duplicate_genes,
    'orphan_transcripts': check_orphan_transcripts,
    'transcript_gene_id_mismatch': check_transcript_gene_id_mismatch,
    'invalid_cds_bounds': check_invalid_cds_bounds,
    'protein_length_vs_cds': check_protein_length_vs_cds,
    'orphan_orthology_genes': check_orphan_orthology_genes,
    'spliceindomains_orphans': check_spliceindomains_orphans,
    'proteins_without_domains': check_proteins_without_domains,
}


def run_anomaly_checks(con, out_dir=None):
    """
    Run all local anomaly checks. Returns a dict of {check_name: DataFrame}.
    If `out_dir` is given, each non-empty result is also written to
    `<out_dir>/<check_name>.csv`.
    """
    results = {}
    for name, check in ANOMALY_CHECKS.items():
        logger.info(f"Running check: {name}")
        result = check(con)
        if isinstance(result, dict):
            results.update(result)
        else:
            results[name] = result

    for name, df in results.items():
        logger.info(f"  {name}: {len(df)} row(s)")
        if out_dir and not df.empty:
            df.to_csv(f"{out_dir.rstrip('/')}/{name}.csv", index=False)

    return results


# ---------------------------------------------------------------------------
# Random-sample comparison against Ensembl / NCBI (RefSeq)
# ---------------------------------------------------------------------------

class RateLimitedSession:
    """
    Thin wrapper around `urllib` that enforces a minimum delay between
    requests and retries on 429/5xx with the server-provided Retry-After
    (or exponential backoff) - so we stay within Ensembl's and NCBI's public
    rate limits. Uses only the standard library, no extra dependencies.
    """

    def __init__(self, min_interval, max_retries=3, timeout=15):
        self.min_interval = min_interval
        self.max_retries = max_retries
        self.timeout = timeout
        self._last_request = 0.0

    def get_json(self, url, headers=None):
        """GET `url` and parse the response body as JSON, or return None on a non-retryable error."""
        headers = headers or {}
        for attempt in range(self.max_retries + 1):
            elapsed = time.monotonic() - self._last_request
            if elapsed < self.min_interval:
                time.sleep(self.min_interval - elapsed)
            self._last_request = time.monotonic()

            request = urllib.request.Request(url, headers=headers)
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    return json.loads(response.read())
            except urllib.error.HTTPError as exc:
                if exc.code == 429 or exc.code >= 500:
                    retry_after = float(exc.headers.get('Retry-After', 2 ** attempt))
                    logger.warning(f"{url} -> HTTP {exc.code}, retrying in {retry_after}s")
                    time.sleep(retry_after)
                    continue
                if exc.code == 404:
                    return None
                logger.warning(f"{url} -> HTTP {exc.code}: {exc.reason}")
                return None
            except urllib.error.URLError as exc:
                logger.warning(f"{url} -> {exc}, retrying")
                time.sleep(2 ** attempt)
                continue
        return None


def _fetch_ensembl_gene(session, gene_ensembl_id):
    """Look up a gene (with its transcripts) on the Ensembl REST API."""
    base_id = gene_ensembl_id.split('.')[0]
    url = f"https://rest.ensembl.org/lookup/id/{base_id}?expand=1;content-type=application/json"
    return session.get_json(url)


def _fetch_ensembl_transcript(session, transcript_ensembl_id):
    """Look up a transcript (with exons/translation) on the Ensembl REST API."""
    base_id = transcript_ensembl_id.split('.')[0]
    url = f"https://rest.ensembl.org/lookup/id/{base_id}?expand=1;content-type=application/json"
    return session.get_json(url)


def _fetch_ncbi_gene_summary(session, gene_id, api_key=None, email=None):
    """Look up a gene's summary on NCBI via E-utilities (esummary)."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={gene_id}&retmode=json"
    if api_key:
        url += f"&api_key={api_key}"
    if email:
        url += f"&email={email}"
    result = session.get_json(url)
    if result is None:
        return None
    try:
        return result['result'][str(gene_id)]
    except (KeyError, ValueError):
        return None


def compare_random_genes_with_external_dbs(
    con, n, seed=None, ensembl_min_interval=0.34, ncbi_min_interval=0.35,
    ncbi_api_key=None, ncbi_email=None,
):
    """
    Pick `n` random genes from the DB and compare their locally-stored
    transcript/exon/protein/domain data against Ensembl and NCBI.

    Ensembl's public REST API allows ~15 requests/second (no key needed);
    NCBI E-utilities allows 3 requests/second without an API key, or 10/s
    with one (pass `ncbi_api_key`/`ncbi_email`). The default intervals stay
    comfortably under both limits.

    Returns a DataFrame, one row per checked gene, with the local values,
    the external values, and a `mismatches` column listing what differs.
    Genes whose species aren't in SPECIE_TO_ENSEMBL, or that have no
    canonical transcript, are skipped (and noted in the `notes` column).
    """
    rng = random.Random(seed)

    genes = pd.read_sql_query(
        "SELECT gene_GeneID_id, gene_ensembl_id, gene_symbol, chromosome, strand, specie FROM Genes",
        con,
    )
    sample = genes.sample(n=min(n, len(genes)), random_state=rng.randint(0, 2**32 - 1))

    ensembl_session = RateLimitedSession(min_interval=ensembl_min_interval)
    ncbi_session = RateLimitedSession(min_interval=ncbi_min_interval)

    rows = []
    for _, gene in sample.iterrows():
        row = {
            'gene_GeneID_id': gene.gene_GeneID_id,
            'gene_ensembl_id': gene.gene_ensembl_id,
            'gene_symbol': gene.gene_symbol,
            'specie': gene.specie,
            'notes': '',
            'mismatches': '',
        }
        mismatches = []
        notes = []

        # --- local data for this gene -------------------------------------------------
        transcripts = pd.read_sql_query(
            "SELECT * FROM Transcripts WHERE gene_GeneID_id = ?", con, params=(gene.gene_GeneID_id,)
        )
        # canonical field: 0=not canonical, 1=RefSeq only, 2=Ensembl only, 3=both
        canonical = transcripts[transcripts.canonical != 0]
        row['local_transcript_count'] = len(transcripts)
        row['local_canonical_count'] = len(canonical)

        # --- Ensembl comparison ----------------------------------------------------
        ensembl_specie = SPECIE_TO_ENSEMBL.get(gene.specie)
        if not gene.gene_ensembl_id:
            notes.append('no gene_ensembl_id - skipping Ensembl comparison')
        elif not ensembl_specie:
            notes.append(f"unknown species '{gene.specie}' - skipping Ensembl comparison")
        else:
            try:
                ensembl_gene = _fetch_ensembl_gene(ensembl_session, gene.gene_ensembl_id)
            except Exception as exc:
                ensembl_gene = None
                notes.append(f"Ensembl lookup failed: {exc}")

            if ensembl_gene is None:
                notes.append('gene not found on Ensembl')
            else:
                row['ensembl_chromosome'] = ensembl_gene.get('seq_region_name')
                row['ensembl_strand'] = '+' if ensembl_gene.get('strand', 1) > 0 else '-'
                row['ensembl_transcript_count'] = len(ensembl_gene.get('Transcript', []))

                if str(row['ensembl_chromosome']) != str(gene.chromosome):
                    mismatches.append('chromosome')
                if row['ensembl_strand'] != gene.strand:
                    mismatches.append('strand')

                # Compare the canonical transcript's exon count, if we can match it.
                for _, t in canonical.iterrows():
                    if not t.transcript_ensembl_id:
                        continue
                    base_tid = t.transcript_ensembl_id.split('.')[0]
                    ensembl_t = next(
                        (et for et in ensembl_gene.get('Transcript', [])
                         if et.get('id') == base_tid),
                        None,
                    )
                    if ensembl_t is None:
                        notes.append(f"canonical transcript {t.transcript_ensembl_id} not found on Ensembl")
                        continue
                    ensembl_exon_count = len(ensembl_t.get('Exon', []))
                    row['local_canonical_exon_count'] = t.exon_count
                    row['ensembl_canonical_exon_count'] = ensembl_exon_count
                    if t.exon_count is not None and t.exon_count != ensembl_exon_count:
                        mismatches.append('canonical_exon_count')

        # --- NCBI / RefSeq comparison -----------------------------------------------
        if gene.gene_GeneID_id and str(gene.gene_GeneID_id).isdigit():
            try:
                ncbi_summary = _fetch_ncbi_gene_summary(
                    ncbi_session, gene.gene_GeneID_id, api_key=ncbi_api_key, email=ncbi_email,
                )
            except Exception as exc:
                ncbi_summary = None
                notes.append(f"NCBI lookup failed: {exc}")

            if ncbi_summary is None:
                notes.append('gene not found on NCBI')
            else:
                ncbi_symbol = ncbi_summary.get('name')
                row['ncbi_symbol'] = ncbi_symbol
                if ncbi_symbol and ncbi_symbol != gene.gene_symbol:
                    mismatches.append('gene_symbol')
        else:
            notes.append('no numeric gene_GeneID_id - skipping NCBI comparison')

        row['mismatches'] = '; '.join(mismatches)
        row['notes'] = '; '.join(notes)
        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', required=True, help='Path to the SQLite database.')
    subparsers = parser.add_subparsers(dest='command', required=True)

    anomalies_parser = subparsers.add_parser('anomalies', help='Run local SQL-based anomaly checks.')
    anomalies_parser.add_argument('--out-dir', help='Directory to write CSV reports to (one per non-empty check).')

    compare_parser = subparsers.add_parser('compare', help='Compare a random sample of genes against Ensembl/NCBI.')
    compare_parser.add_argument('-n', type=int, default=10, help='Number of random genes to check.')
    compare_parser.add_argument('--seed', type=int, default=None, help='Random seed for reproducibility.')
    compare_parser.add_argument('--out', help='Path to write the comparison CSV to.')
    compare_parser.add_argument('--ncbi-api-key', default=None)
    compare_parser.add_argument('--ncbi-email', default=None)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')

    con = sqlite3.connect(args.db)

    if args.command == 'anomalies':
        results = run_anomaly_checks(con, out_dir=args.out_dir)
        for name, df in results.items():
            print(f"{name}: {len(df)} row(s)")
    elif args.command == 'compare':
        df = compare_random_genes_with_external_dbs(
            con, args.n, seed=args.seed,
            ncbi_api_key=args.ncbi_api_key, ncbi_email=args.ncbi_email,
        )
        if args.out:
            df.to_csv(args.out, index=False)
        print(df.to_string())


if __name__ == '__main__':
    main()
