"""
Generate gene visualization PDF similar to DoChap-web JavaScript visualization.
This module creates PDF reports showing gene transcripts with exons and protein domains.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.gridspec import GridSpec
import matplotlib.patheffects as mpe
import numpy as np
import pandas as pd
import sqlite3

import utils


# Color palette from DoChap-web
DOCHAP_COLORS = [
    "#FF4A46", "#1BE177", "#00CCFF", "#A30059", "#FFB500", "#006FA6", "#4FC601", "#FFDBE5",
    "#D16100", "#00C2A0", "#A079BF", "#C0B9B2", "#CC0744", "#549E79", "#B79762", "#B903AA",
    "#00846F", "#FF90C9", "#0AA6D8", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#BEC459",
    "#AA5199", "#0089A3", "#EEC3FF", "#8FB0FF", "#004D43", "#F4D749", "#997D87", "#3B5DFF",
    "#FF2F80", "#6B7900", "#FFAA92", "#A1C299", "#885578", "#B77B68", "#FAD09F", "#456D75",
    "#FF8A9A", "#0086ED", "#D157A0", "#00A6AA", "#B4A8BD", "#FF913F", "#636375", "#A3C8C9"
]

DOMAIN_LABEL_MAX_LEN = 14
PDF_RASTER_DPI = 110
DOMAIN_STRIP_COUNT = 80
_ANALYSED_MARK_COLOR = '#1a365d'  # analysed/drawn-only chip inside a domain ellipse

SPECIES_ALIASES = {
    'h_sapiens': 'H_sapiens',
    'human': 'H_sapiens',
    'm_musculus': 'M_musculus',
    'mouse': 'M_musculus',
    'rat': 'R_norvegicus',
    'r_norvegicus': 'R_norvegicus',   
}


def _route_representative_domain_id_to_column(domain_id):
    """Map an InterPro match-XML accession (RepresentativeDomains.domain_id) to the
    DomainType-style identifier column _format_domain_label() looks for. Only affects
    which bucket the id is displayed under, not whether a domain is drawn."""
    if domain_id.startswith('IPR'):
        return 'interpro'
    if domain_id.startswith('PF'):
        return 'pfam'
    if domain_id.startswith('SM'):
        return 'smart'
    if domain_id.startswith('TIGR'):
        return 'tigr'
    if domain_id.lower().startswith('cd'):
        # A 'cd*' accession (e.g. cd21801) is a CDD *domain* accession, so it
        # belongs in the 'cdd' bucket that _format_domain_label() displays -
        # not the numeric-PSSM 'CDD_id' bucket, which the label never reads
        # (that would drop the domain's id from the figure entirely).
        return 'cdd'
    return 'interpro'


_REPRESENTATIVE_DOMAIN_ID_COLUMNS = ('CDD_id', 'cdd', 'pfam', 'smart', 'tigr', 'interpro')


def _representative_domains_to_domain_columns(df_rep):
    """Reshape raw RepresentativeDomains rows (protein_interpro_id, domain_id, domain_name,
    start, end, score) into the AA_start/AA_end + identifier-column shape that
    _format_domain_label()/the protein view expect from DomainEvent/DomainType rows."""
    if df_rep.empty:
        return pd.DataFrame()

    df_rep = df_rep.dropna(subset=['start', 'end']).copy()
    if df_rep.empty:
        return df_rep

    domain_column = df_rep['domain_id'].map(_route_representative_domain_id_to_column)
    for col in _REPRESENTATIVE_DOMAIN_ID_COLUMNS:
        df_rep[col] = df_rep['domain_id'].where(domain_column == col)

    df_rep = df_rep.rename(columns={'start': 'AA_start', 'end': 'AA_end'})
    df_rep['AA_start'] = df_rep['AA_start'].astype(int)
    df_rep['AA_end'] = df_rep['AA_end'].astype(int)
    return df_rep


def _drop_undrawable_domains(df_domains):
    """Drop domains whose AA_start/AA_end aren't finite numbers.

    Such a domain can't be placed on the protein axis: a NaN AA_start yields a
    NaN-geometry ellipse whose empty path crashes matplotlib's tight-bbox when
    the PDF is saved (np.concatenate([])). ~1% of DoChaP's DomainEvent rows have
    a NULL coordinate, so this keeps one bad row from breaking a gene's figure.
    """
    if df_domains is None or len(df_domains) == 0:
        return df_domains
    if 'AA_start' not in df_domains.columns or 'AA_end' not in df_domains.columns:
        return df_domains
    start = pd.to_numeric(df_domains['AA_start'], errors='coerce')
    end = pd.to_numeric(df_domains['AA_end'], errors='coerce')
    keep = np.isfinite(start.to_numpy()) & np.isfinite(end.to_numpy())
    if keep.all():
        return df_domains
    return df_domains[keep].reset_index(drop=True)


def _keep_primary_entries(df_domains):
    """Keep only the InterPro Domain and Repeat entries - tier 1 of the ladder.

    The same set filter_representative_domains() compares, so the drawing and the
    events table never disagree about which domains exist. Applied here rather
    than at draw time so the protein axis, the row layout and the figure height
    are all sized to what is actually drawn.

    A frame the ladder cannot rank - DomainEvent/DomainType rows, which state no
    entry type - is returned untouched: there is no type to select on, and
    filtering on an absent column would silently blank those transcripts.
    """
    if df_domains is None or len(df_domains) == 0:
        return df_domains
    # Imported here, not at module scope: junction_analisys imports this module.
    from junction_analisys import domain_entry_tiers, TIER_PRIMARY
    tiers = domain_entry_tiers(df_domains)
    if tiers is None:
        return df_domains
    keep = (tiers == TIER_PRIMARY).to_numpy()
    if keep.all():
        return df_domains
    return df_domains[keep].reset_index(drop=True)


def _drawable_domains(df_domains):
    """The domains a figure draws: those that can be placed on the protein axis
    (_drop_undrawable_domains) and that the figure admits (_keep_primary_entries)."""
    return _keep_primary_entries(_drop_undrawable_domains(df_domains))




def prepare_gene_data_bulk(conn, gene_ensembl_ids):
    """
    Bulk-load gene/transcript/exon/protein/domain data for many genes at once.

    Domains come from RepresentativeDomains and nowhere else; a protein with no
    entry there has no domains, so the drawing matches the set the analysis
    compared. DomainEvent/DomainType are not read - the analysis does not use them.

    Returns a dict keyed by lowercased gene ensembl id, where each value is a dict
    with 'gene_data' (a Series, as `df_gene.iloc[0]`) and 'transcripts' (a list
    of {'info', 'exons', 'domains'} dicts), suitable to pass as `preloaded` to
    `GeneVisualization`. Genes not found in the database are simply absent from
    the result.
    """
    unique_ids = list(dict.fromkeys(s for s in gene_ensembl_ids if isinstance(s, str) and s))
    if not unique_ids:
        return {}

    placeholders = ','.join(['?'] * len(unique_ids))
    df_genes = pd.read_sql_query(
        f"SELECT * FROM Genes WHERE gene_ensembl_id COLLATE NOCASE IN ({placeholders})",
        conn, params=unique_ids,
    )
    if len(df_genes) == 0:
        return {}

    df_genes = df_genes.drop_duplicates(subset='gene_ensembl_id', keep='first')
    gene_ensembl_ids = df_genes['gene_ensembl_id'].tolist()
    ge_placeholders = ','.join(['?'] * len(gene_ensembl_ids))
    df_transcripts = pd.read_sql_query(
        f"""SELECT * FROM Transcripts
            WHERE gene_ensembl_id IN ({ge_placeholders})
            ORDER BY tx_start""",
        conn, params=gene_ensembl_ids,
    )

    # Combine ensembl and refseq ids so RefSeq-only transcripts (no ensembl id)
    # still get their exons/proteins matched, instead of being silently dropped.
    df_transcripts['combined_id'] = df_transcripts['transcript_ensembl_id'].fillna(
        df_transcripts['transcript_refseq_id']
    )
    transcript_ids = df_transcripts['combined_id'].dropna().tolist()
    if transcript_ids:
        tx_placeholders = ','.join(['?'] * len(transcript_ids))
        df_exons_all = pd.read_sql_query(
            f"""SELECT * FROM Transcript_exon
                WHERE transcript_ensembl_id IN ({tx_placeholders})
                   OR transcript_refseq_id IN ({tx_placeholders})""",
            conn, params=transcript_ids * 2,
        )
        df_proteins_all = pd.read_sql_query(
            f"""SELECT * FROM Proteins
                WHERE transcript_ensembl_id IN ({tx_placeholders})
                   OR transcript_refseq_id IN ({tx_placeholders})""",
            conn, params=transcript_ids * 2,
        )
        if len(df_exons_all):
            df_exons_all['combined_id'] = df_exons_all['transcript_ensembl_id'].fillna(
                df_exons_all['transcript_refseq_id']
            )
        if len(df_proteins_all):
            df_proteins_all['combined_id'] = df_proteins_all['transcript_ensembl_id'].fillna(
                df_proteins_all['transcript_refseq_id']
            )
    else:
        df_exons_all = pd.DataFrame()
        df_proteins_all = pd.DataFrame()

    protein_ids = df_proteins_all['protein_ensembl_id'].dropna().unique().tolist() if len(df_proteins_all) else []
    if protein_ids:
        pr_placeholders = ','.join(['?'] * len(protein_ids))
        df_domains_all = pd.read_sql_query(
            f"""SELECT de.*, dt.*
                FROM DomainEvent de
                JOIN DomainType dt ON de.type_id = dt.type_id
                WHERE de.protein_ensembl_id IN ({pr_placeholders})
                ORDER BY de.AA_start""",
            conn, params=protein_ids,
        )
    else:
        df_domains_all = pd.DataFrame()

    exons_by_transcript = {k: v for k, v in df_exons_all.groupby('combined_id')} if len(df_exons_all) else {}
    proteins_by_transcript = {k: v for k, v in df_proteins_all.groupby('combined_id')} if len(df_proteins_all) else {}
    domains_by_protein = {k: v for k, v in df_domains_all.groupby('protein_ensembl_id')} if len(df_domains_all) else {}

    # RepresentativeDomains and nothing else - a protein absent from that table is
    # drawn with no domains, so the PDF shows the domain set the analysis compared.
    domains_by_protein = {}

    if protein_ids and 'protein_interpro_id' in df_proteins_all.columns:
        df_protein_interpro = df_proteins_all[['protein_ensembl_id', 'protein_interpro_id']].dropna(subset=['protein_interpro_id'])
        df_protein_interpro = df_protein_interpro[df_protein_interpro.protein_interpro_id.str.strip() != '']
        interpro_ids = df_protein_interpro['protein_interpro_id'].unique().tolist()
        if interpro_ids:
            ip_placeholders = ','.join(['?'] * len(interpro_ids))
            try:
                df_rep_all = pd.read_sql_query(
                    f"SELECT * FROM RepresentativeDomains WHERE protein_interpro_id IN ({ip_placeholders})",
                    conn, params=interpro_ids,
                )
            except (sqlite3.OperationalError, pd.errors.DatabaseError):
                df_rep_all = pd.DataFrame()
            if len(df_rep_all):
                df_rep_all = df_rep_all.merge(df_protein_interpro, on='protein_interpro_id')
                df_rep_all = _representative_domains_to_domain_columns(df_rep_all)
                for protein_id, df_rep_protein in df_rep_all.groupby('protein_ensembl_id'):
                    domains_by_protein[protein_id] = df_rep_protein.reset_index(drop=True)

    def get_exons(transcript_id):
        df_exons = exons_by_transcript.get(transcript_id, pd.DataFrame(columns=df_exons_all.columns)).reset_index(drop=True)
        if len(df_exons) == 0:
            return df_exons
        if (df_exons['abs_start_CDS'] == 0).any():
            non_coding = df_exons[df_exons['abs_start_CDS'] == 0].sort_values('genomic_start_tx')
            coding = df_exons[df_exons['abs_start_CDS'] > 0].sort_values('abs_start_CDS')
            df_exons = pd.concat([non_coding, coding], ignore_index=True)
        else:
            df_exons = df_exons.sort_values('abs_start_CDS').reset_index(drop=True)
        return df_exons

    domains_columns = df_domains_all.columns if len(protein_ids) else []

    def get_domains(transcript_id):
        df_protein = proteins_by_transcript.get(transcript_id)
        if df_protein is None or len(df_protein) == 0:
            return pd.DataFrame()
        protein_id = df_protein.iloc[0]['protein_ensembl_id']
        domains = domains_by_protein.get(protein_id, pd.DataFrame(columns=domains_columns)).reset_index(drop=True)
        return _drawable_domains(domains)

    transcript_groups = {k: v for k, v in df_transcripts.groupby('gene_ensembl_id', sort=False)} if len(df_transcripts) else {}

    result = {}
    for _, gene_row in df_genes.iterrows():
        transcripts = []
        for _, transcript in transcript_groups.get(gene_row['gene_ensembl_id'], pd.DataFrame()).iterrows():
            tid = transcript['combined_id']
            transcripts.append({
                'info': transcript,
                'exons': get_exons(tid),
                'domains': get_domains(tid),
            })
        result[gene_row['gene_ensembl_id']] = {'gene_data': gene_row, 'transcripts': transcripts}

    return result


class GeneVisualization:
    """Class to create gene visualization similar to DoChap-web."""
    
    def __init__(self, conn, gene_name, preloaded=None):
        """
        Initialize gene visualization.

        Parameters:
        -----------
        conn : sqlite3.Connection
            Database connection to DoChap database
        gene_name : str
            Gene symbol to visualize
        preloaded : dict | None
            Optional pre-fetched data for this gene, as returned per-gene by
            `prepare_gene_data_bulk()`. When provided, `load_gene_data()` uses
            it instead of querying the database.
        """
        self.conn = conn
        self.gene_name = gene_name
        self.gene_data = None
        self.transcripts = []
        self.colors = {}
        self.color_index = 0
        self._preloaded = preloaded
        self.species_hint = None
        # Optional calibrated fold-change probability per transcript (vs canonical),
        # keyed by transcript_ensembl_id -> float in [0,1] (P(TM<0.5), the
        # structural companion score). Shown on the transcript label when present.

    def _normalize_species_value(self, species_value):
        """Normalize common species aliases to the values used in the DB."""
        if species_value is None or pd.isna(species_value):
            return None
        text = str(species_value).strip()
        if not text:
            return None
        if text in SPECIES_ALIASES:
            return SPECIES_ALIASES[text]
        return text

    def _infer_species_hint_from_dataframe(self, df_junction):
        """Read species hints from a junction DataFrame when available."""
        if df_junction is None or len(df_junction) == 0:
            return 'H_sapiens'  # Default to human when no data is available
        if 'specie' not in df_junction.columns:
            return 'H_sapiens'
        value = df_junction['specie'].dropna().unique()[0]
        normalized_value = self._normalize_species_value(value)
        return normalized_value   
            
    def load_gene_data(self):
        """Load gene data, from preloaded data if available, otherwise from the database."""
        if self._preloaded is not None:
            self.gene_data = self._preloaded['gene_data']
            self.transcripts = self._preloaded['transcripts']
            self._assign_exon_colors()
            return

        # Get gene info
        base_query = """
            SELECT * FROM Genes
            WHERE gene_symbol = ? COLLATE NOCASE
        """
        params = [self.gene_name]
        if self.species_hint is not None:
            query = base_query + " AND specie = ?"
            params.append(self.species_hint)
            df_gene = pd.read_sql_query(query, self.conn, params=params)
            if len(df_gene) == 0:
                df_gene = pd.read_sql_query(base_query, self.conn, params=[self.gene_name])
        else:
            df_gene = pd.read_sql_query(base_query, self.conn, params=[self.gene_name])

        if len(df_gene) == 0:
            raise ValueError(f"Gene '{self.gene_name}' not found in database")

        self.gene_data = df_gene.iloc[0]

        # Get transcripts for this gene
        trans_query = """
            SELECT * FROM Transcripts
            WHERE gene_ensembl_id = ?
            ORDER BY tx_start
        """
        df_transcripts = pd.read_sql_query(trans_query, self.conn,
                                           params=[self.gene_data['gene_ensembl_id']])

        # Load exons and domains for each transcript
        for _, transcript in df_transcripts.iterrows():
            transcript_data = {
                'info': transcript,
                'exons': self._load_exons(transcript['transcript_ensembl_id'], transcript.get('transcript_refseq_id')),
                'domains': self._load_domains(transcript['transcript_ensembl_id'], transcript.get('transcript_refseq_id'))
            }
            self.transcripts.append(transcript_data)

        # Assign colors to exons
        self._assign_exon_colors()

    def _load_exons(self, transcript_ensembl_id, transcript_refseq_id=None):
        """Load exons for a transcript, matched by either its ensembl or refseq id
        (a RefSeq-only transcript has no ensembl id, so both are checked)."""
        exon_query = """
            SELECT * FROM Transcript_exon
            WHERE transcript_ensembl_id = ? OR transcript_refseq_id = ?
            ORDER BY abs_start_CDS
        """
        df_exons = pd.read_sql_query(exon_query, self.conn, params=[transcript_ensembl_id, transcript_refseq_id])
        # If any exons have abs_start_CDS = 0 (non-coding), sort those by genomic position
        # and put them at the beginning
        if (df_exons['abs_start_CDS'] == 0).any():
            non_coding = df_exons[df_exons['abs_start_CDS'] == 0].sort_values('genomic_start_tx')
            coding = df_exons[df_exons['abs_start_CDS'] > 0].sort_values('abs_start_CDS')
            df_exons = pd.concat([non_coding, coding], ignore_index=True)
        return df_exons

    def _load_domains(self, transcript_ensembl_id, transcript_refseq_id=None):
        """Load protein domains for a transcript, matched by either its ensembl or
        refseq id (a RefSeq-only transcript has no ensembl id).

        Domains come from RepresentativeDomains only; a protein with no entry
        there has no domains. DomainEvent/DomainType are not read - the analysis
        does not use them, so drawing from them would show a set no comparison
        could refer to.
        """
        # Get protein ID
        protein_query = """
            SELECT protein_ensembl_id, protein_interpro_id FROM Proteins
            WHERE transcript_ensembl_id = ? OR transcript_refseq_id = ?
        """
        df_protein = pd.read_sql_query(protein_query, self.conn, params=[transcript_ensembl_id, transcript_refseq_id])

        if len(df_protein) == 0:
            return pd.DataFrame()

        df_rep = pd.DataFrame()
        protein_interpro_id = df_protein.iloc[0].get('protein_interpro_id')
        if pd.notna(protein_interpro_id) and str(protein_interpro_id).strip():
            rep_query = "SELECT * FROM RepresentativeDomains WHERE protein_interpro_id = ?"
            try:
                df_rep = pd.read_sql_query(rep_query, self.conn, params=[protein_interpro_id])
            except (sqlite3.OperationalError, pd.errors.DatabaseError):
                df_rep = pd.DataFrame()
            df_rep = _representative_domains_to_domain_columns(df_rep)
        if len(df_rep) == 0:
            return pd.DataFrame()
        return _drawable_domains(df_rep.reset_index(drop=True))
    
    def _assign_exon_colors(self):
        """Assign colors to unique exons across all transcripts based on genomic location."""
        seen_exons = {}
        
        for transcript in self.transcripts:
            for _, exon in transcript['exons'].iterrows():
                key = (int(exon['genomic_start_tx']), int(exon['genomic_end_tx']))
                
                if key not in seen_exons:
                    seen_exons[key] = DOCHAP_COLORS[self.color_index % len(DOCHAP_COLORS)]
                    self.color_index += 1
        
        self.colors = seen_exons
    
    def _get_exon_color(self, exon):
        """Get color for an exon, keyed on genomic location."""
        key = (int(exon['genomic_start_tx']), int(exon['genomic_end_tx']))

        color = self.colors.get(key)
        gray_like = {'#cccccc', '#808080', '#a9a9a9', 'gray', 'grey', 'darkgray', 'darkgrey'}
        if color is not None and str(color).strip().lower() not in gray_like:
            return color

        # Deterministic non-gray fallback.
        idx = abs(hash(key)) % len(DOCHAP_COLORS)
        return DOCHAP_COLORS[idx]

    def _format_domain_label(self, domain_row, compact_mode=False, max_len=DOMAIN_LABEL_MAX_LEN):
        """Return the single id shown for one domain in the figure.

        Unlike results.csv - where the analysis deliberately *collapses* a
        domain's cross-references down to one display name - the figure draws
        every domain separately, so each domain gets exactly one id here with
        no collapsing: no '*' merge/truncation marker, no dropped ids.

        The id is picked interpro-first: the columns below are scanned in the
        same order of preference the analysis uses when it collapses
        (junction_analisys.DOMAIN_NAME_PREFIX_PRIORITY = IPR, pfam, cd, smart,
        tigr), so the figure label and the results-table name agree on which id
        represents a domain. 'CDD_id' is intentionally excluded - it holds a
        numeric CDD PSSM id, not a human-facing accession (the 'cdd' column
        already carries the cd##### accession, including for RepresentativeDomains
        via _route_representative_domain_id_to_column()).
        """
        for col in ['interpro', 'pfam', 'cdd', 'smart', 'tigr']:
            if col in domain_row and pd.notna(domain_row[col]):
                value = str(domain_row[col]).strip()
                if value and value.lower() != 'nan':
                    return value

        return None
    
    def _format_number(self, num):
        """Format number with commas."""
        return f"{int(num):,}"


    def _is_negative_strand(self):
        """Return True when gene strand indicates reverse genomic orientation."""
        if self.gene_data is None or 'strand' not in self.gene_data:
            return False
        strand_value = self.gene_data['strand']
        if pd.isna(strand_value):
            return False

        strand_text = str(strand_value).strip().lower()
        return strand_text == '-' or strand_text.startswith('-') or strand_text == 'minus'

    def _transcript_produces_protein(self, transcript):
        """Return True when transcript has a protein identifier."""
        protein_id = transcript['info'].get('protein_ensembl_id')
        if protein_id is None or pd.isna(protein_id):
            return False
        protein_text = str(protein_id).strip()
        return bool(protein_text) and protein_text.lower() != 'nan'

    def _transcript_has_protein(self, transcript):
        """True when DoChaP annotates a protein for this transcript, by either id.

        Not _transcript_produces_protein(), which reads the ensembl id alone: a
        RefSeq-only transcript (XM_/XP_) has a protein under protein_refseq_id and
        would otherwise read as non-coding. Same test as the analysis's
        coding_by_transcript.
        """
        for column in ('protein_ensembl_id', 'protein_refseq_id'):
            value = transcript['info'].get(column)
            if value is None or pd.isna(value):
                continue
            text = str(value).strip()
            if text and text.lower() not in ('nan', 'none'):
                return True
        return False

    def _transcript_has_domains(self, transcript):
        """Return True when transcript protein has at least one domain row."""
        return len(transcript['domains']) > 0

    def _transcript_matches_ids(self, transcript, transcript_ids):
        """Return True when this transcript's ensembl or refseq id is in transcript_ids."""
        ensembl_id = transcript['info'].get('transcript_ensembl_id')
        refseq_id = transcript['info'].get('transcript_refseq_id')
        return (
            (pd.notna(ensembl_id) and ensembl_id in transcript_ids)
            or (pd.notna(refseq_id) and refseq_id in transcript_ids)
        )

    def _transcript_ids(self, transcript):
        """A transcript's ensembl and refseq ids, the pair every match here is made
        on - DOMAS keys a transcript by whichever of the two it has."""
        return [i for i in (transcript['info'].get('transcript_ensembl_id'),
                            transcript['info'].get('transcript_refseq_id'))
                if i is not None and pd.notna(i)]

    def _analysis_for(self, transcript, by_transcript):
        """The analysis's entry for this transcript, or None where it has none -
        a transcript the analysis never evaluated, or a standalone PDF with no
        analysis at all. Keyed by either id, as the analysis keys transcripts."""
        if not by_transcript:
            return None
        for transcript_id in self._transcript_ids(transcript):
            if transcript_id in by_transcript:
                return by_transcript[transcript_id]
        return None

    def _transcript_reading_order(self, df_results, canonical_transcript_id=None):
        """Sort key putting the figure in the order the analysis reasons about it:
        the canonical transcript first, then the event groups in order, and within
        a group the transcript the selection rule picked - most-like-canonical
        ahead of longest-CDS - ahead of the ones it passed over. Transcripts
        belonging to no group (no unique feature, or not carrying the event at all)
        come last.

        `canonical_transcript_id` is the transcript DOMAS compared against, which
        is not always the one the DB flags: a gene can carry two canonical
        transcripts - one flagged by RefSeq and one by Ensembl - and _resolve_
        canonical() picks between them. Sorting on the DB flag alone drew the
        other one first and left the actual reference in the middle of the figure,
        unmarked. Falls back to the flag where no id is given.
        """
        rows_by_id = {}
        if df_results is not None and 'transcript_id' in df_results.columns:
            for transcript_id, rows in df_results.groupby('transcript_id'):
                rows_by_id[transcript_id] = rows

        def key(transcript):
            ids = self._transcript_ids(transcript)
            if canonical_transcript_id is not None:
                is_canonical = canonical_transcript_id in ids
            else:
                is_canonical = bool(transcript['info'].get('canonical'))
            if is_canonical:
                return (0, 0, 0, '')

            group, rank = None, 2
            for transcript_id in ids:
                rows = rows_by_id.get(transcript_id)
                if rows is None or len(rows) == 0:
                    continue
                if 'group' in rows.columns and rows['group'].notna().any():
                    value = rows['group'].dropna().min()
                    group = value if group is None else min(group, value)
                if rows.get('is_most_like_canonical') is not None and rows['is_most_like_canonical'].any():
                    rank = min(rank, 0)
                elif rows.get('is_longest_cds') is not None and rows['is_longest_cds'].any():
                    rank = min(rank, 1)
            # No group: drawn after every group, in the same relative order as before.
            return (1, float(group) if group is not None else float('inf'), rank,
                    ids[0] if ids else '')

        return key

    def _compute_domain_label_positions(self, label_items, axis_max, base_y, lane_step=0.08, lanes=4):
        """Assign below-domain label positions across lanes to reduce overlap."""
        if not label_items:
            return []

        # Keep deterministic ordering for stable output.
        ordered = sorted(label_items, key=lambda item: item['center'])
        lane_right_edges = [-1e12] * lanes
        placed = []

        for item in ordered:
            # Approximate label footprint in axis units from text length.
            est_half_width = max(item['width'] * 0.25, len(item['text']) * 2.2)
            left_edge = item['center'] - est_half_width
            right_edge = item['center'] + est_half_width

            lane_idx = None
            for candidate in range(lanes):
                if left_edge >= lane_right_edges[candidate]:
                    lane_idx = candidate
                    break

            if lane_idx is None:
                lane_idx = min(range(lanes), key=lambda idx: lane_right_edges[idx])

            lane_right_edges[lane_idx] = right_edge + max(3.0, axis_max * 0.004)
            placed.append({
                'center': item['center'],
                'text': item['text'],
                'label_y': base_y - lane_idx * lane_step,
            })

        return placed

    def _domain_row_layout(self, domains_df, max_protein_length):
        """
        Assign each domain in domains_df to a display row (0-based) so that
        domains sharing a row never overlap - classic greedy interval-
        scheduling row placement (sort by start, place each interval in the
        first row whose last-placed interval already ends before this one
        starts).

        Uses each domain's *effective* (post min-width-padding) span rather
        than its raw AA_start/AA_end, matching what _draw_protein_view()
        actually renders - otherwise a narrow domain padded up to
        min_domain_width could still visually collide with a "different
        row" neighbor that was only clear of its raw (unpadded) span.

        Both create_pdf() (to size each transcript's GridSpec row before
        drawing anything) and _draw_protein_view() (to actually place
        ellipses) call this, so planning and drawing always agree on row
        count/assignment.

        Returns (row_of: {domains_df index -> row number}, num_rows: int).
        Returns ({}, 1) for an empty/None domains_df, so callers can use
        num_rows directly as a height multiplier without a special case.
        """
        if domains_df is None or len(domains_df) == 0:
            return {}, 1

        min_domain_width = max(8.0, max_protein_length * 0.02)
        effective = {}
        for idx, domain in domains_df.iterrows():
            start = float(domain['AA_start'])
            end = float(domain['AA_end'])
            width = max(end - start, min_domain_width)
            center = (start + end) / 2
            effective[idx] = (center - width / 2, center + width / 2)

        order = sorted(effective, key=lambda i: effective[i])
        row_last_end = []
        row_of = {}
        for idx in order:
            start, end = effective[idx]
            row = next((r for r, last_end in enumerate(row_last_end) if start >= last_end), None)
            if row is None:
                row_last_end.append(end)
                row_of[idx] = len(row_last_end) - 1
            else:
                row_last_end[row] = end
                row_of[idx] = row

        return row_of, len(row_last_end)

    def _get_coding_exon_segments(self, transcript):
        """Return coding exon segments in amino-acid coordinate space."""
        segments = []
        for exon_idx, (_, exon) in enumerate(transcript['exons'].iterrows()):
            if exon['abs_start_CDS'] > 0 and exon['abs_end_CDS'] >= exon['abs_start_CDS']:
                start_nuc = int(exon['abs_start_CDS'])
                end_nuc = int(exon['abs_end_CDS'])
                # Convert nucleotide CDS coordinates to amino-acid axis.
                start_aa = (start_nuc - 1) / 3.0
                end_aa = end_nuc / 3.0
                color = self._get_exon_color(exon)
                if color == '#CCCCCC':
                    color = DOCHAP_COLORS[exon_idx % len(DOCHAP_COLORS)]
                segments.append({
                    'start_aa': start_aa,
                    'end_aa': end_aa,
                    'color': color,
                })
        return sorted(segments, key=lambda segment: segment['start_aa'])

    def _prepare_junction_display_df(self, df_junction):
        """Return a copy of junction metadata with stable idx/color columns."""
        if df_junction is None or len(df_junction) == 0:
            return None

        display_df = df_junction.copy().reset_index(drop=True)
        display_df['idx'] = np.arange(1, len(display_df) + 1)
        display_df['junction_color'] = [
            DOCHAP_COLORS[i % len(DOCHAP_COLORS)] for i in range(len(display_df))
        ]
        return display_df

    def _normalize_junction_df(self, df_junction):
        """Return a junction DataFrame with normalized start/end columns."""
        if df_junction is None or len(df_junction) == 0:
            return None

        normalized = df_junction.copy()
        if 'start' not in normalized.columns and 'start_position' in normalized.columns:
            normalized['start'] = normalized['start_position']
        if 'end' not in normalized.columns and 'end_position' in normalized.columns:
            normalized['end'] = normalized['end_position']
        if 'start' not in normalized.columns or 'end' not in normalized.columns:
            return None

        normalized['start'] = pd.to_numeric(normalized['start'], errors='coerce')
        normalized['end'] = pd.to_numeric(normalized['end'], errors='coerce')
        return normalized

    def _is_transcript_canonical(self, transcript):
        """Check if the transcript is marked canonical."""
        canonical_val = transcript['info'].get('canonical')
        if canonical_val is None or pd.isna(canonical_val):
            return False
        return int(canonical_val) != 0

    def _filter_junctions_for_transcript(self, transcript, junction_items):
        """Keep only event features that are truly supported by this transcript.

        Delegates to utils.find_matching_junction_indices, the same predicate the
        analysis uses, so the visualization and the analysis cannot disagree about
        which features are mapped.

        A feature carrying utils.FEATURE_RETAINED_INTRON is matched by containment
        (one exon spanning the intron) rather than by adjacency, so a retained-intron
        isoform is drawn against the feature it actually supports.
        """
        if not junction_items:
            return []

        pairs, feature_types = [], []
        keep_positions = []
        for position, junction in enumerate(junction_items):
            if pd.isna(junction.get('start')) or pd.isna(junction.get('end')):
                continue
            pairs.append((int(junction['start']), int(junction['end'])))
            feature_types.append(junction.get('feature_type') or utils.FEATURE_JUNCTION)
            keep_positions.append(position)

        if not pairs:
            return []

        strand = '-' if self._is_negative_strand() else '+'
        matched = utils.find_matching_junction_indices(
            transcript['exons'], pairs, strand=strand, feature_types=feature_types)
        return [junction_items[keep_positions[i]] for i in sorted(matched)]

    def _junction_items(self, df_junction):
        """The cluster's features as drawable items, before any per-transcript
        selection. Split out of _get_matching_junctions() so a caller holding the
        analysis's own answer can select from the same items instead of matching
        them again."""
        normalized = self._normalize_junction_df(df_junction)
        if normalized is None or len(normalized) == 0:
            return []

        items = []
        for _, junction in normalized.iterrows():
            if pd.isna(junction['start']) or pd.isna(junction['end']):
                continue
            items.append({
                'start': int(junction['start']),
                'end': int(junction['end']),
                'idx': int(junction['idx']) if 'idx' in junction and pd.notna(junction['idx']) else None,
                'color': junction['junction_color'] if 'junction_color' in junction else 'red',
                # Carried through so _filter_junctions_for_transcript() can apply the
                # containment predicate to a retained intron. Without it every feature
                # defaults to a junction and the retaining transcript - the one the
                # feature exists to identify - is drawn as not supporting it.
                'feature_type': (str(junction[utils.FEATURE_TYPE_COLUMN])
                                 if utils.FEATURE_TYPE_COLUMN in junction
                                 and pd.notna(junction[utils.FEATURE_TYPE_COLUMN])
                                 else utils.FEATURE_JUNCTION),
            })

        return items

    def _get_matching_junctions(self, transcript, df_junction):
        """The features this transcript carries, matched here. Used where no
        analysis result is available - a standalone PDF, or a transcript the
        analysis never evaluated."""
        return self._filter_junctions_for_transcript(
            transcript, self._junction_items(df_junction))

    def _draw_junction_table(self, ax, df_junction):
        """Draw a compact first-page table for provided junction metadata."""
        ax.axis('off')
        if df_junction is None or len(df_junction) == 0:
            return

        table_df = df_junction.copy()
        if 'junction_color' in table_df.columns:
            junction_colors = list(table_df['junction_color'])
            table_df = table_df.drop(columns=['junction_color'])
        else:
            junction_colors = ['black'] * len(table_df)
        if 'idx' in table_df.columns:
            idx_column = list(table_df.columns).index('idx')
            front_cols = ['idx'] + [col for col in table_df.columns if col != 'idx']
            table_df = table_df[front_cols]
            idx_column = 0
        else:
            idx_column = None
        if 'cluster_name' in table_df.columns:
            table_df = table_df.drop(columns=['cluster_name'])
        if table_df.empty:
            return

        display_df = table_df.fillna('')
        column_names = list(display_df.columns)
        table = ax.table(
            cellText=display_df.astype(str).values,
            colLabels=column_names,
            cellLoc='center',
            bbox=[0, 0, 1, 1],
        )
        table.auto_set_font_size(False)
        table.set_fontsize(5.5)

        width_map = {
            'idx': 0.035,
            'start': 0.095,
            'end': 0.095,
        }
        flexible_columns = [name for name in column_names if name not in width_map]
        remaining_width = max(0.10, 1.0 - sum(width_map.get(name, 0.0) for name in column_names))
        default_width = remaining_width / max(1, len(flexible_columns))

        for (row, col), cell in table.get_celld().items():
            column_name = column_names[col]
            cell.set_width(width_map.get(column_name, default_width))
            cell.set_linewidth(0.6)
            cell.set_edgecolor('#888888')
            if row == 0:
                cell.set_facecolor('#F1F1F1')
                cell.set_text_props(weight='bold')
            elif idx_column is not None and col == idx_column:
                cell.set_text_props(color=junction_colors[row - 1], weight='bold')

    def _draw_results_table(self, ax, df_results):
        """Draw a compact first-page table summarizing domain comparison results."""
        ax.axis('off')
        if df_results is None or len(df_results) == 0:
            return

        # astype(object) first: the results frame carries a nullable-Int64 'group'
        # column, and filling pd.NA with '' in place would raise on its dtype.
        table_df = df_results.astype(object).fillna('')
        column_names = list(table_df.columns)
        table = ax.table(
            cellText=table_df.astype(str).values,
            colLabels=column_names,
            cellLoc='center',
            bbox=[0, 0, 1, 1],
        )
        table.auto_set_font_size(False)
        table.set_fontsize(5.5)

        # Relative widths, normalised over the columns actually present (the
        # optional ones come and go with the run's flags). The identifying
        # columns are what a reader scans; the counts and lengths are at most a
        # few characters and were getting the same share as a transcript id.
        weights = {
            'event': 3.4, 'transcript_id': 3.4, 'domain_id': 1.9, 'rank': 1.5,
            'group': 0.6,
            'c_domain_length': 0.85, 't_domain_length': 0.85,
            'c_domains_number': 0.85, 't_domains_number': 0.85,
            'c_junction_in_cds': 1.0, 't_junction_in_cds': 1.0,
        }
        total = sum(weights.get(name, 1.0) for name in column_names)

        for (row, col), cell in table.get_celld().items():
            cell.set_width(weights.get(column_names[col], 1.0) / total)
            cell.set_linewidth(0.6)
            cell.set_edgecolor('#888888')
            if row == 0:
                cell.set_facecolor('#F1F1F1')
                cell.set_text_props(weight='bold')






    def _draw_genomic_junctions(self, ax, junctions, exon_y, exon_height):
        """Draw colored bracket-like junction markers above the genomic exon track."""
        if not junctions:
            return

        baseline_top = exon_y + exon_height / 2
        baseline_bottom = exon_y - exon_height / 2
        for index, junction in enumerate(junctions):
            left = min(junction['start'], junction['end'])
            right = max(junction['start'], junction['end'])
            color = junction.get('color', 'red')
            linestyle = junction.get('linestyle', 'solid')

            if junction.get('feature_type') == utils.FEATURE_RETAINED_INTRON:
                # A retained intron is CONTAINED in one exon rather than spliced
                # across: drawn as a raised bracket it would read as the excision
                # it is the opposite of. Marked instead as a dotted span inside the
                # exon body, bounded by dotted ticks at the intron's own edges, so
                # it stays within the exon's height and disturbs no layout.
                ax.plot([left, right], [exon_y, exon_y], color=color, linewidth=1.6,
                        linestyle='dotted', zorder=6)
                for edge in (left, right):
                    ax.plot([edge, edge], [baseline_bottom, baseline_top], color=color,
                            linewidth=1.0, linestyle='dotted', zorder=6)
                continue

            junction_top = min(0.95, baseline_top + 0.10 + index * 0.08)
            ax.plot([left, left], [baseline_top, junction_top], color=color, linewidth=1.4,
                    linestyle=linestyle, zorder=5)
            ax.plot([right, right], [baseline_top, junction_top], color=color, linewidth=1.4,
                    linestyle=linestyle, zorder=5)
            ax.plot([left, right], [junction_top, junction_top], color=color, linewidth=1.4,
                    linestyle=linestyle, zorder=5)
    
    def create_pdf(self, output_file='gene_visualization.pdf', transcripts_per_page=4,
                   protein_only=False, domains_only=False, df_junction=None,
                   df_results=None, show_canonical_non_relevant_junctions=True,
                   transcript_ids=None, no_comparison_note=None,
                   canonical_transcript_id=None, analysis_domains=None,
                   analysis_features=None):
        """
        Create PDF visualization of the gene, one page per transcripts_per_page transcripts.
        Each page has its own axis scales at the top so they never overlap with transcript rows.

        Parameters:
        -----------
        output_file : str
            Path to save the PDF.
        transcripts_per_page : int
            Number of transcript rows to draw per page.
        protein_only : bool
            If True, include only transcripts with a protein identifier.
        domains_only : bool
            If True, include only transcripts whose protein has at least one domain.
        df_junction : pandas.DataFrame | None
            Optional junction metadata table with at least `start` and `end` columns.
        df_results : pandas.DataFrame | None
            Optional domain comparison results table, displayed on the first page.
        show_canonical_non_relevant_junctions : bool
            If True, the canonical transcript shows all junctions with dashed lines
            for those that are not relevant to that transcript.
        transcript_ids : set | list | None
            If provided, only draw transcripts whose ensembl or refseq id is in
            this collection - typically the canonical transcript plus the
            transcripts that were actually compared to it. All other transcripts
            of the gene are omitted entirely. If None (default), every transcript
            of the gene is eligible (subject to protein_only/domains_only).
        no_comparison_note : str | None
            If given, this text is drawn as an extra line below the last
            transcript on the last page - typically used together with
            transcript_ids to explain that no other transcript qualified for
            comparison against the canonical one.
        """
        from matplotlib.backends.backend_pdf import PdfPages

        if df_junction is not None and len(df_junction) > 0:
            inferred_species = self._infer_species_hint_from_dataframe(df_junction)
            if inferred_species is not None:
                self.species_hint = inferred_species

        if self.gene_data is None:
            self.load_gene_data()

        junction_display_df = self._prepare_junction_display_df(df_junction)

        # Skip empty/invalid transcript entries, then order them the way the
        # analysis reasons about them (see _transcript_reading_order).
        valid_transcripts = [t for t in self.transcripts if len(t['exons']) > 0]
        valid_transcripts.sort(
            key=self._transcript_reading_order(df_results, canonical_transcript_id))
        if transcript_ids is not None:
            valid_transcripts = [t for t in valid_transcripts if self._transcript_matches_ids(t, transcript_ids)]
        if protein_only:
            valid_transcripts = [t for t in valid_transcripts if self._transcript_produces_protein(t)]
        if domains_only:
            valid_transcripts = [t for t in valid_transcripts if self._transcript_has_domains(t)]
        num_transcripts = len(valid_transcripts)
        if num_transcripts == 0:
            if protein_only and domains_only:
                print(
                    f"No transcripts with protein and domains found for gene {self.gene_name}"
                )
            elif domains_only:
                print(f"No domain-containing transcripts found for gene {self.gene_name}")
            elif protein_only:
                print(f"No protein-producing transcripts found for gene {self.gene_name}")
            elif transcript_ids is not None:
                print(f"None of the specified transcripts were found for gene {self.gene_name}")
            else:
                print(f"No transcripts found for gene {self.gene_name}")
            return

        canonical_transcript = next(
            (t for t in valid_transcripts if t['info'].get('canonical')), None)

        genomic_start = min(t['info']['tx_start'] for t in valid_transcripts)
        genomic_end   = max(t['info']['tx_end']   for t in valid_transcripts)

        # Calculate max AA axis from all transcripts once.
        max_protein_length = 0
        for transcript in valid_transcripts:
            protein_length = transcript['info'].get('protein_length')
            if protein_length is not None and not pd.isna(protein_length):
                max_protein_length = max(max_protein_length, float(protein_length))
            coding_segments = self._get_coding_exon_segments(transcript)
            if coding_segments:
                max_protein_length = max(max_protein_length,
                                         max(s['end_aa'] for s in coding_segments))
            if len(transcript['domains']) > 0:
                max_protein_length = max(max_protein_length,
                                         float(transcript['domains']['AA_end'].max()))
        if max_protein_length <= 0:
            max_protein_length = 1.0

        # Every domain the ladder can rank is marked, whether or not this figure
        # demotes one: the tier is worth reading in its own right, and marking only
        # demoted figures would leave an unmarked domain meaning either "kept" or
        # "unrankable", with no way to tell which.
        ladder_marks_by_transcript = {
            id(t): self.domain_ladder_marks(t, self._analysis_for(t, analysis_domains))
            for t in valid_transcripts
        }
        marks_active = any(v for v in ladder_marks_by_transcript.values())
        # Under DomainEvent/DomainType the rows state no entry type, so the ladder
        # cannot rank them and their domains carry no mark. Next to a ranked
        # transcript that needs explaining, or it reads as an inconsistency rather
        # than a different annotation source.
        marks_partial = marks_active and any(
            ladder_marks_by_transcript[id(t)] is None and len(t['domains']) > 0
            for t in valid_transcripts)

        gene_title = (
            f"{self.gene_name} - {self.gene_data['specie']}  |  "
            f"chr{self.gene_data['chromosome']}:"
            f"{self._format_number(genomic_start)}-{self._format_number(genomic_end)}  |  "
            f"{num_transcripts} transcript(s)"
        )
        if junction_display_df is not None and 'cluster_name' in junction_display_df.columns:
            cluster_values = junction_display_df['cluster_name'].dropna()
            if len(cluster_values) > 0:
                gene_title += f"  |  {cluster_values.iloc[0]}"

        # Split into pages
        pages = [valid_transcripts[i:i + transcripts_per_page]
                 for i in range(0, num_transcripts, transcripts_per_page)]
        num_pages = len(pages)

        # Cluster-level events (not tied to a specific transcript) shown once above canonical
        _CLUSTER_EVENTS = {'feature_not_mapped', 'no_canonical_features',
                           'no_canonical_transcript', 'only_one_transcript',
                           'no_gene_specified'}
        cluster_events_df = None
        if df_results is not None and 'event' in df_results.columns:
            _mask = df_results['event'].isin(_CLUSTER_EVENTS)
            if _mask.any():
                _ce = df_results[_mask][['event', 'transcript_id']].copy()
                _ce['transcript_id'] = _ce.apply(
                    lambda r: f"junction #{int(r['transcript_id'])}"
                              if r['event'] == 'feature_not_mapped' and pd.notna(r['transcript_id'])
                              else '',
                    axis=1,
                )
                cluster_events_df = _ce.rename(columns={'transcript_id': 'details'})

        with PdfPages(output_file) as pdf:
            for page_idx, page_transcripts in enumerate(pages):
                n = len(page_transcripts)
                show_junction_table = page_idx == 0 and junction_display_df is not None and len(junction_display_df) > 0

                # Per-transcript results tables, sized to fit their row count
                transcript_results = []
                for transcript in page_transcripts:
                    transcript_id = transcript['info']['transcript_ensembl_id']
                    transcript_refseq_id = transcript['info'].get('transcript_refseq_id')
                    if df_results is not None and 'transcript_id' in df_results.columns:
                        ids_to_match = [i for i in [transcript_id, transcript_refseq_id]
                                        if i is not None and not pd.isna(i)]
                        rows = df_results[df_results['transcript_id'].isin(ids_to_match)]
                    else:
                        rows = None
                    transcript_results.append(rows)

                # Height layout per page:
                #   row 0  : page title
                #   row 1  : optional junction table on first page
                #   next   : axis scales
                #   next   : spacer
                #   next   : optional cluster-level events table (page 1 only)
                #   next   : spacer after cluster table (page 1 only)
                #   rows.. : n * (results table row + exon/genomic row + protein/domain row + label row)
                show_cluster_events = page_idx == 0 and cluster_events_df is not None
                height_ratios = [0.55]
                if show_junction_table:
                    height_ratios.append(1.50)
                scale_row = len(height_ratios)
                height_ratios += [0.75, 0.28]
                if show_cluster_events:
                    _num_ce = len(cluster_events_df)
                    _ce_height = 0.42 + 0.30 * _num_ce
                    cluster_events_row = scale_row + 2
                    height_ratios += [_ce_height, 0.20]
                    transcript_start_row = cluster_events_row + 2
                else:
                    transcript_start_row = scale_row + 2
                for transcript, rows in zip(page_transcripts, transcript_results):
                    num_rows = len(rows) if rows is not None else 0
                    results_height = 0.18 if num_rows == 0 else 0.42 + 0.30 * num_rows
                    # A transcript with several overlapping domains needs
                    # more than one domain row (see _domain_row_layout()) -
                    # grow its genomic/protein row here so the extra rows
                    # get real page space instead of being squeezed into a
                    # fixed box and spilling into the next transcript.
                    _, num_domain_rows = self._domain_row_layout(transcript['domains'], max_protein_length)
                    protein_row_height = 1.2 + 0.55 * (num_domain_rows - 1)
                    height_ratios += [results_height, 0.7, protein_row_height, 0.18]
                show_no_comparison_note = (
                    page_idx == num_pages - 1 and no_comparison_note
                )
                if show_no_comparison_note:
                    note_row = len(height_ratios)
                    height_ratios.append(0.30)
                total_rows = len(height_ratios)

                fig = plt.figure(figsize=(11, 8.5))
                gs = GridSpec(
                    total_rows, 2,
                    figure=fig,
                    height_ratios=height_ratios,
                    width_ratios=[1.2, 1],
                    hspace=0.30,
                    wspace=0.35,
                    left=0.08,
                    right=0.95,
                    top=0.97,
                    bottom=0.04,
                )

                # ── title row ──────────────────────────────────────────────
                ax_title = fig.add_subplot(gs[0, :])
                ax_title.axis('off')
                page_label = f"  (page {page_idx + 1}/{num_pages})" if num_pages > 1 else ""
                ax_title.text(0.5, 0.5, gene_title + page_label,
                              ha='center', va='center',
                              fontsize=11, fontweight='bold',
                              transform=ax_title.transAxes)

                if show_junction_table:
                    ax_table = fig.add_subplot(gs[1, :])
                    self._draw_junction_table(ax_table, junction_display_df)

                if show_cluster_events:
                    ax_cluster = fig.add_subplot(gs[cluster_events_row, :])
                    self._draw_results_table(ax_cluster, cluster_events_df)

                # ── axis scale row (repeated on every page) ────────────────
                self._draw_genomic_scale(fig, gs[scale_row, 0], genomic_start, genomic_end)
                self._draw_protein_scale(fig, gs[scale_row, 1], max_protein_length)

                # ── transcript rows ────────────────────────────────────────
                stride = 4
                for i, transcript in enumerate(page_transcripts):
                    row = transcript_start_row + i * stride

                    ax_results = fig.add_subplot(gs[row, :])
                    # is_longest_cds/is_most_like_canonical show as a tag on the
                    # label below, not as True/False table columns.
                    table_rows = transcript_results[i]
                    if table_rows is not None:
                        table_rows = table_rows.drop(columns=['is_longest_cds', 'is_most_like_canonical'], errors='ignore')
                    self._draw_results_table(ax_results, table_rows)

                    grow = row + 1

                    ax_genomic = fig.add_subplot(gs[grow:grow + 2, 0])
                    self._draw_genomic_view(
                        ax_genomic,
                        transcript,
                        genomic_start,
                        genomic_end,
                        df_junction=junction_display_df,
                        show_canonical_non_relevant_junctions=show_canonical_non_relevant_junctions,
                        analysis_features=analysis_features,
                    )

                    ax_protein = fig.add_subplot(gs[grow:grow + 2, 1])
                    self._draw_protein_view(
                        ax_protein, transcript, max_protein_length,
                        ladder_marks=(ladder_marks_by_transcript.get(id(transcript))
                                      if marks_active else None),
                    )

                    ax_label = fig.add_subplot(gs[grow + 2, :])
                    ax_label.axis('off')
                    transcript_name = transcript['info']['transcript_ensembl_id']
                    if transcript_name is None or pd.isna(transcript_name):
                        transcript_name = transcript['info'].get('transcript_refseq_id')
                    protein_name = transcript['info'].get('protein_ensembl_id')
                    if protein_name is None or pd.isna(protein_name) or str(protein_name).strip() == '':
                        protein_name = transcript['info'].get('protein_refseq_id')
                    if protein_name is None or pd.isna(protein_name) or str(protein_name).strip() == '':
                        protein_name = 'N/A'
                    rows_for_transcript = transcript_results[i]
                    tie_break_tags = []
                    # Named first, so the reference the whole figure is measured
                    # against is identifiable rather than merely first.
                    if canonical_transcript_id is not None and \
                            canonical_transcript_id in self._transcript_ids(transcript):
                        tie_break_tags.append('canonical')
                    if not self._transcript_has_protein(transcript):
                        tie_break_tags.append('no protein')
                    if rows_for_transcript is not None and len(rows_for_transcript) > 0:
                        if 'is_longest_cds' in rows_for_transcript.columns and rows_for_transcript['is_longest_cds'].any():
                            tie_break_tags.append('longest CDS')
                        if 'is_most_like_canonical' in rows_for_transcript.columns and rows_for_transcript['is_most_like_canonical'].any():
                            tie_break_tags.append('most like canonical')
                    tag_suffix = f"  [{', '.join(tie_break_tags)}]" if tie_break_tags else ""
                    ax_label.text(
                        0.02, 0.98,
                        f"Transcript: {transcript_name}  |  Protein: {protein_name}{tag_suffix}",
                        fontsize=8, va='top', transform=ax_label.transAxes,
                    )

                if show_no_comparison_note:
                    ax_note = fig.add_subplot(gs[note_row, :])
                    ax_note.axis('off')
                    ax_note.text(
                        0.02, 0.5, no_comparison_note,
                        fontsize=9, fontstyle='italic', va='center', transform=ax_note.transAxes,
                    )

                if marks_active:
                    legend = (
                        "Only InterPro Domain and Repeat entries are drawn, and only those are "
                        "compared — so every domain the events table names appears here. Neither "
                        "drawn nor compared: InterPro Family and Homologous_superfamily entries, the "
                        "residue features (active/binding/conserved site, PTM), and member-database "
                        "signatures (G3DSA/PTHR/SSF/cd/PF), which are not curated structural units. "
                        "FILLED = compared. HOLLOW = removed before the comparison, so no row of the "
                        "events table can refer to it — an entry that duplicates a longer kept entry "
                        "of the same accession overlapping it."
                    )
                    if marks_partial:
                        legend += (
                            " A transcript whose domains carry NO mark has no RepresentativeDomains "
                            "entry: its annotation comes from DomainEvent/DomainType, which states no "
                            "entry type, so the ladder can neither rank nor select on it - every one "
                            "of its domains is drawn, and compared as-is."
                        )
                    fig.text(
                        0.02, 0.004, legend,
                        fontsize=5.6, style='italic', color='#555555', va='bottom', wrap=True,
                    )
                pdf.savefig(fig, bbox_inches='tight', dpi=PDF_RASTER_DPI)
                plt.close(fig)

        print(f"PDF saved to: {output_file}  ({num_pages} page(s))")
    
    def _draw_genomic_scale(self, fig, gs_position, start, end):
        """Draw genomic coordinate scale with alternating above/below labels."""
        ax = fig.add_subplot(gs_position)
        left = min(start, end)
        right = max(start, end)
        if self._is_negative_strand():
            ax.set_xlim(right, left)
        else:
            ax.set_xlim(left, right)
        ax.set_ylim(-0.5, 1.2)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.set_xlabel('Genomic Position (bp)', fontsize=9, labelpad=18)
        ax.grid(axis='x', alpha=0.3, linestyle='--')

        # Compute tick positions then draw labels alternating above/below axis.
        locator = plt.MaxNLocator(nbins=8, prune='both')
        locator.set_axis(ax.xaxis)
        ticks = locator.tick_values(left, right)
        ticks = [t for t in ticks if left <= t <= right]
        ax.set_xticks(ticks)
        ax.set_xticklabels([])  # hide default labels
        ax.tick_params(axis='x', direction='out', length=4)

        for i, t in enumerate(ticks):
            label = f"{int(round(t)):,}"
            if i % 2 == 0:
                # above the axis line
                ax.text(t, 0.15, label, ha='center', va='bottom', fontsize=6,
                        transform=ax.get_xaxis_transform())
            else:
                # below the axis line
                ax.text(t, -0.35, label, ha='center', va='top', fontsize=6,
                        transform=ax.get_xaxis_transform())
        
    def _draw_protein_scale(self, fig, gs_position, max_length):
        """Draw protein coordinate scale in amino acids."""
        ax = fig.add_subplot(gs_position)
        ax.set_xlim(0, max_length)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.xaxis.tick_bottom()
        ax.xaxis.set_label_position('bottom')
        ax.set_xlabel('Amino Acid Position', fontsize=9, labelpad=3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, prune='both'))
        ax.tick_params(axis='x', labelsize=6.5, rotation=40)
    
    def _draw_genomic_view(self, ax, transcript, genomic_start, genomic_end,
                          df_junction=None, show_canonical_non_relevant_junctions=False,
                          analysis_features=None):
        """Draw genomic view with exons."""
        left = min(genomic_start, genomic_end)
        right = max(genomic_start, genomic_end)
        if self._is_negative_strand():
            ax.set_xlim(right, left)
        else:
            ax.set_xlim(left, right)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        
        exon_y = 0.5
        exon_height = 0.3
        non_cds_height = exon_height * 0.5

        cds_start_tx = transcript['info'].get('cds_start') if 'cds_start' in transcript['info'] else None
        cds_end_tx = transcript['info'].get('cds_end') if 'cds_end' in transcript['info'] else None
        has_transcript_cds = (
            cds_start_tx is not None
            and cds_end_tx is not None
            and pd.notna(cds_start_tx)
            and pd.notna(cds_end_tx)
            and int(cds_start_tx) > 0
            and int(cds_end_tx) > 0
        )
        if has_transcript_cds:
            cds_left = min(int(cds_start_tx), int(cds_end_tx))
            cds_right = max(int(cds_start_tx), int(cds_end_tx))
        else:
            cds_left = None
            cds_right = None
        
        # Draw baseline
        ax.plot([left, right], [exon_y, exon_y], 'k-', linewidth=1.5, zorder=1)
        
        # Draw exons
        for exon_idx, (_, exon) in enumerate(transcript['exons'].iterrows()):
            start = min(exon['genomic_start_tx'], exon['genomic_end_tx'])
            end = max(exon['genomic_start_tx'], exon['genomic_end_tx'])
            width = end - start
            exon_has_cds = (
                pd.notna(exon.get('abs_start_CDS'))
                and pd.notna(exon.get('abs_end_CDS'))
                and int(exon['abs_start_CDS']) > 0
                and int(exon['abs_end_CDS']) >= int(exon['abs_start_CDS'])
            )
            
            # Get exon color
            color = self._get_exon_color(exon)
            if color == '#CCCCCC':
                color = DOCHAP_COLORS[exon_idx % len(DOCHAP_COLORS)]
            
            if not has_transcript_cds:
                # No coding region for this transcript: draw exon as reduced-height non-CDS.
                if exon_has_cds:
                    coding_rect = Rectangle(
                        (start, exon_y - exon_height / 2),
                        width,
                        exon_height,
                        facecolor=color,
                        edgecolor='none',
                        linewidth=0,
                        zorder=4,
                    )
                    ax.add_patch(coding_rect)
                    continue
                gray_rect = Rectangle(
                    (start, exon_y - non_cds_height / 2),
                    width,
                    non_cds_height,
                    facecolor=color,
                    edgecolor='none',
                    alpha=0.35,
                    zorder=3,
                )
                ax.add_patch(gray_rect)
                continue

            coding_start = max(start, cds_left)
            coding_end = min(end, cds_right)

            if coding_end <= coding_start:
                # Exon is outside CDS: draw reduced-height non-CDS.
                if exon_has_cds:
                    coding_rect = Rectangle(
                        (start, exon_y - exon_height / 2),
                        width,
                        exon_height,
                        facecolor=color,
                        edgecolor='none',
                        linewidth=0,
                        zorder=4,
                    )
                    ax.add_patch(coding_rect)
                    continue
                gray_rect = Rectangle(
                    (start, exon_y - non_cds_height / 2),
                    width,
                    non_cds_height,
                    facecolor=color,
                    edgecolor='none',
                    alpha=0.35,
                    zorder=3,
                )
                ax.add_patch(gray_rect)
                continue

            # Draw non-CDS subsegments (5' and/or 3' UTR) at half height.
            if start < coding_start:
                left_non_cds = Rectangle(
                    (start, exon_y - non_cds_height / 2),
                    coding_start - start,
                    non_cds_height,
                    facecolor=color,
                    edgecolor='none',
                    alpha=0.35,
                    zorder=3,
                )
                ax.add_patch(left_non_cds)

            if coding_end < end:
                right_non_cds = Rectangle(
                    (coding_end, exon_y - non_cds_height / 2),
                    end - coding_end,
                    non_cds_height,
                    facecolor=color,
                    edgecolor='none',
                    alpha=0.35,
                    zorder=3,
                )
                ax.add_patch(right_non_cds)

            # Draw CDS part at full height with exon color.
            coding_rect = Rectangle(
                (coding_start, exon_y - exon_height / 2),
                coding_end - coding_start,
                exon_height,
                facecolor=color,
                edgecolor='none',
                linewidth=0,
                zorder=4,
            )
            ax.add_patch(coding_rect)

        normalized = self._normalize_junction_df(df_junction)
        if normalized is None or len(normalized) == 0:
            return

        if self._is_transcript_canonical(transcript) and show_canonical_non_relevant_junctions:
            all_junctions = []
            for _, junction in normalized.iterrows():
                if pd.isna(junction['start']) or pd.isna(junction['end']):
                    continue
                all_junctions.append({
                    'start': int(junction['start']),
                    'end': int(junction['end']),
                    'idx': int(junction['idx']) if 'idx' in junction and pd.notna(junction['idx']) else None,
                    'color': junction['junction_color'] if 'junction_color' in junction else 'red',
                    'linestyle': 'dashed',
                    # Needed both to match a retained intron by containment and to
                    # draw it as a span; without it this branch silently treats
                    # every feature on the canonical transcript as a junction.
                    'feature_type': (str(junction[utils.FEATURE_TYPE_COLUMN])
                                     if utils.FEATURE_TYPE_COLUMN in junction
                                     and pd.notna(junction[utils.FEATURE_TYPE_COLUMN])
                                     else utils.FEATURE_JUNCTION),
                })
            relevant_junctions = self._filter_junctions_for_transcript(transcript, all_junctions)
            relevant_ids = {(item['start'], item['end']) for item in relevant_junctions}
            for item in all_junctions:
                if (item['start'], item['end']) in relevant_ids:
                    item['linestyle'] = 'solid'
            self._draw_genomic_junctions(ax, all_junctions, exon_y, exon_height)
        else:
            analysis_matched = self._analysis_for(transcript, analysis_features)
            if analysis_matched is not None:
                # Which features a transcript carries is a decision the analysis
                # already made, with the same predicate; re-running it here risked
                # the drawing and the events table disagreeing. Selected by
                # coordinate rather than by position: the analysis and the drawing
                # build their feature lists separately, and one skips rows the
                # other keeps.
                wanted = {(int(start), int(end)) for start, end in analysis_matched}
                matched_junctions = [item for item in self._junction_items(normalized)
                                     if (item['start'], item['end']) in wanted]
            else:
                matched_junctions = self._get_matching_junctions(transcript, normalized)
            self._draw_genomic_junctions(ax, matched_junctions, exon_y, exon_height)

    
    @staticmethod
    def _domain_key(df, index):
        """Identity of a domain row across two independent loads of the same DB
        rows: the analysis and the drawing each read them separately, so the frame
        index cannot be compared - the accession and its residue span can."""
        return (str(df.at[index, 'domain_id']),
                df.at[index, 'AA_start'], df.at[index, 'AA_end'])

    def domain_ladder_marks(self, transcript, analysis_kept=None):
        """`{row key: (tier, was_kept)}` for `transcript`'s domains, or None when the
        frame carries no InterPro entry types to rank by.

        `analysis_kept` is the frame the analysis kept for this transcript. Given
        one, `was_kept` is read straight off it - the drawing states what the
        comparison did rather than re-deciding it. Without one (a standalone PDF,
        or a transcript the analysis never evaluated) the filter is run here, which
        is the only way to answer at all.

        The frame reaching here has already been reduced to tier 1 by
        _keep_primary_entries(), so every mark reads '1'; what it still carries is
        `was_kept`. compare_domains() runs filter_representative_domains() over the
        same frame before comparing anything, and that filter drops a Domain/Repeat
        entry duplicating a longer one of the same accession - such an entry is
        drawn hollow, and no row of the events table refers to it. `was_kept` is
        read back from the filter's own output rather than recomputed, so the mark
        states what the analysis did.
        """
        domains = transcript.get('domains')
        if domains is None or len(domains) == 0:
            return None
        # Imported here, not at module scope: junction_analisys imports this module.
        from junction_analisys import domain_entry_tiers, filter_representative_domains
        tiers = domain_entry_tiers(domains)
        if tiers is None:
            return None
        if analysis_kept is not None and len(analysis_kept) >= 0 and 'domain_id' in getattr(analysis_kept, 'columns', []):
            kept_keys = {self._domain_key(analysis_kept, i) for i in analysis_kept.index}
            return {key: (tiers[key], self._domain_key(domains, key) in kept_keys)
                    for key in domains.index}
        kept = set(filter_representative_domains(domains).index)
        return {key: (tiers[key], key in kept) for key in domains.index}

    def _draw_protein_view(self, ax, transcript, max_protein_length, ladder_marks=None):
        """Draw protein/domain view on the right (protein above domains).

        `ladder_marks`, when given, labels each domain with its ladder tier and
        whether it was kept - see domain_ladder_marks().
        """
        ax.set_xlim(0, max_protein_length)
        # Keep labels/borders as vector, but rasterize dense protein/domain fills.
        ax.set_rasterization_zorder(2.5)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])

        protein_y = 0.72
        protein_height = 0.18
        domain_y_top = 0.32  # y of the first (topmost) domain row's ellipse center
        domain_height = 0.25
        domain_row_step = 0.42  # vertical distance between successive domain rows
        min_domain_width = max(8.0, max_protein_length * 0.02)

        row_of, num_domain_rows = self._domain_row_layout(transcript['domains'], max_protein_length)


        # Extra rows push the bottom of the axis down so a transcript with
        # several overlapping domains gets a taller box (see create_pdf(),
        # which sizes the GridSpec row the same way) instead of its extra
        # rows spilling into/being masked by the next transcript.
        domain_bottom = -0.35 - (num_domain_rows - 1) * domain_row_step
        ax.set_ylim(domain_bottom, 1)

        protein_length_aa = transcript['info'].get('protein_length')
        if protein_length_aa is None or pd.isna(protein_length_aa):
            protein_length_aa = max_protein_length
        protein_length_aa = min(max_protein_length, float(protein_length_aa))

        # Draw protein backbone rectangle and color it by contributing coding exons.
        coding_segments = self._get_coding_exon_segments(transcript)

        # A transcript with no annotated protein is drawn grey and labelled. It
        # still HAS a bar: DoChaP fills cds_start/cds_end and the exons' CDS
        # offsets for a non-coding transcript as though it were coding, so the
        # geometry exists and can even be longer than a real neighbouring protein
        # - mouse RGS19's ENSMUST00000143510.8 draws 277 aa against the 195 aa of
        # the transcript actually selected. Colouring it as a protein invites the
        # reader to compare lengths that are not comparable; the analysis does not
        # (protein-coding candidates come first in the selection).
        has_protein = self._transcript_has_protein(transcript)

        protein_bg = Rectangle(
            (0, protein_y - protein_height / 2),
            protein_length_aa,
            protein_height,
            facecolor='#F2F2F2' if not has_protein else 'white',
            edgecolor='#AAAAAA' if not has_protein else 'black',
            linewidth=1.2,
            zorder=1,
            rasterized=True,
        )
        ax.add_patch(protein_bg)

        for segment in coding_segments:
            seg_start = max(0, float(segment['start_aa']))
            seg_end = min(protein_length_aa, float(segment['end_aa']))
            if seg_end <= seg_start:
                continue
            seg_rect = Rectangle(
                (seg_start, protein_y - protein_height / 2),
                seg_end - seg_start,
                protein_height,
                facecolor='#DDDDDD' if not has_protein else segment['color'],
                edgecolor='none',
                alpha=1.0,
                zorder=2,
                rasterized=True,
            )
            ax.add_patch(seg_rect)

        if not has_protein:
            ax.text(protein_length_aa / 2, protein_y, 'no annotated protein',
                    ha='center', va='center', fontsize=6, style='italic',
                    color='#777777', zorder=3)

        if len(transcript['domains']) == 0:
            ax.text(0.5, 0.62, 'No domains', transform=ax.transAxes,
                   ha='center', va='center', fontsize=10, style='italic', color='black')
            return

        domains_sorted = transcript['domains'].sort_values('AA_end', ascending=False)

        # Draw each domain as ellipse with gradient, one row per overlap
        # group (row_of, from _domain_row_layout() above) so overlapping
        # domains no longer draw on top of one another - every domain gets
        # its own ellipse and label instead of one being picked to
        # represent the group.
        label_items_by_row = {}
        for domain_idx, (row_key, domain) in enumerate(domains_sorted.iterrows()):
            domain_start_aa = float(domain['AA_start'])
            domain_end_aa = float(domain['AA_end'])
            domain_width_raw = domain_end_aa - domain_start_aa
            domain_width = max(domain_width_raw, min_domain_width)
            domain_center = (domain_start_aa + domain_end_aa) / 2
            domain_y = domain_y_top - row_of[row_key] * domain_row_step

            overlapping_segments = []
            for segment in coding_segments:
                overlap_start = max(domain_start_aa, segment['start_aa'])
                overlap_end = min(domain_end_aa, segment['end_aa'])
                if overlap_start < overlap_end:
                    overlapping_segments.append({
                        'start_aa': overlap_start,
                        'end_aa': overlap_end,
                        'color': segment['color'],
                    })
            
            if len(overlapping_segments) == 0:
                # No exon overlap: use a stable non-gray palette fallback.
                fallback_color = DOCHAP_COLORS[domain_idx % len(DOCHAP_COLORS)]
                ellipse = Ellipse(
                    (domain_center, domain_y),
                    domain_width,
                    domain_height,
                    facecolor=fallback_color,
                    edgecolor='black',
                    linewidth=1.2,
                    alpha=0.8,
                    zorder=2,
                )
                ax.add_patch(ellipse)
            elif len(overlapping_segments) == 1:
                # Single exon, solid color
                color = overlapping_segments[0]['color']
                ellipse = Ellipse(
                    (domain_center, domain_y),
                    domain_width,
                    domain_height,
                    facecolor=color,
                    edgecolor='black',
                    linewidth=1.2,
                    alpha=0.8,
                    zorder=2,
                )
                ax.add_patch(ellipse)
            else:
                # Multiple exons: color the ellipse in amino-acid strips.
                num_strips = DOMAIN_STRIP_COUNT
                strip_width = domain_width / num_strips
                a = domain_width / 2  # Semi-major axis
                b = domain_height / 2  # Semi-minor axis

                strip_patches = []
                for i in range(num_strips):
                    strip_aa_pos = domain_start_aa + (i + 0.5) * strip_width
                    strip_color = overlapping_segments[0]['color']

                    for segment in overlapping_segments:
                        if segment['start_aa'] <= strip_aa_pos <= segment['end_aa']:
                            strip_color = segment['color']
                            break

                    # Calculate ellipse segment height (varies across width)
                    # Using ellipse equation: y = b * sqrt(1 - (x/a)^2)
                    x_from_center = strip_aa_pos - domain_center

                    if abs(x_from_center) < a:
                        # Height at this x position
                        relative_height = np.sqrt(1 - (x_from_center / a) ** 2)
                        strip_height = 2 * b * relative_height

                        # Draw vertical strip
                        rect = Rectangle((strip_aa_pos - strip_width/2, domain_y - strip_height/2),
                                       strip_width * 1.1, strip_height,  # Slight overlap to avoid gaps
                                       facecolor=strip_color)
                        strip_patches.append(rect)

                if strip_patches:
                    strips = PatchCollection(strip_patches, match_original=True, edgecolor='none',
                                             alpha=1.0, zorder=2, rasterized=True)
                    ax.add_collection(strips)

                # Draw border ellipse on top
                ellipse_border = Ellipse((domain_center, domain_y), domain_width, domain_height,
                                       facecolor='none', edgecolor='black', linewidth=1.2, zorder=3)
                ax.add_patch(ellipse_border)
            
            # Ladder mark: a small rectangle inside the ellipse carrying the
            # domain's tier - 1 for an InterPro Domain/Repeat, 2 for a member-DB
            # hit, - for an entry no tier admits - filled when the domain reaches
            # compare_domains() and hollow when the ladder removes it. At 18% of the
            # width it clears the ellipse edge and the centred length label, and the
            # ellipse is over 3/4 of full height there, so the chip fits inside.
            if ladder_marks is not None and row_key in ladder_marks:
                tier, was_kept = ladder_marks[row_key]
                # A fixed fraction of the axis, so the digit stays legible at any
                # domain length; scaling to the domain gives a slab on a 300-aa one
                # and something unreadable on a 30-aa one.
                chip_w = max_protein_length * 0.034
                chip_h = domain_height * 0.55
                chip_x = domain_center - domain_width / 2 + domain_width * 0.18
                ax.add_patch(Rectangle(
                    (chip_x - chip_w / 2, domain_y - chip_h / 2), chip_w, chip_h,
                    facecolor=_ANALYSED_MARK_COLOR if was_kept else 'white',
                    edgecolor=_ANALYSED_MARK_COLOR,
                    linewidth=0.7, zorder=3.6,
                ))
                ax.text(chip_x, domain_y, tier,
                        ha='center', va='center', fontsize=4.2, fontweight='bold',
                        color='white' if was_kept else _ANALYSED_MARK_COLOR,
                        zorder=3.7, clip_on=True)

            # Domain length (aa) printed in the centre of the ellipse; a white
            # stroke keeps it legible over any fill colour.
            ax.text(domain_center, domain_y, str(int(domain_end_aa - domain_start_aa + 1)),
                    ha='center', va='center', fontsize=5.5, fontweight='bold',
                    color='black', zorder=4, clip_on=True,
                    path_effects=[mpe.withStroke(linewidth=1.3, foreground='white')])

            # One label per domain: overlapping domains sit on separate rows, so
            # none is hidden and none has to stand for its group.
            domain_name = self._format_domain_label(domain, compact_mode=False, max_len=DOMAIN_LABEL_MAX_LEN)
            if domain_name:
                label_items_by_row.setdefault(row_of[row_key], []).append({
                    'center': domain_center,
                    'text': domain_name,
                    'width': max(1.0, domain_width_raw),
                })

        for row, label_items in label_items_by_row.items():
            row_domain_y = domain_y_top - row * domain_row_step
            label_base_y = row_domain_y - (domain_height / 2) - 0.05
            placed_labels = self._compute_domain_label_positions(
                label_items,
                max_protein_length,
                label_base_y,
                lane_step=0.08,
                lanes=4,
            )

            for label in placed_labels:
                ax.plot(
                    [label['center'], label['center']],
                    [row_domain_y - domain_height / 2, label['label_y'] + 0.01],
                    color='gray',
                    linewidth=0.5,
                    zorder=4,
                    alpha=0.8,
                )
                ax.text(
                    label['center'],
                    label['label_y'],
                    label['text'],
                    ha='center',
                    va='top',
                    fontsize=6.2,
                    fontweight='bold',
                    style='italic',
                    zorder=5,
                    rotation=32,
                    color='black',
                    clip_on=True,
                )


def generate_gene_pdf(gene_name, conn, output_file=None,
                      protein_only=False, domains_only=False):
    """
    Generate a PDF visualization for a gene similar to DoChap-web.

    Parameters:
    -----------
    gene_name : str
        Gene symbol (e.g., 'PUF60', 'BRCA1')
    conn : sqlite3.Connection
        Open SQLite connection to the DoChap database
    output_file : str, optional
        Path to save the PDF. If None, uses gene_name.pdf
    protein_only : bool, optional
        If True, include only transcripts that produce protein.
    domains_only : bool, optional
        If True, include only transcripts whose protein has at least one domain.

    Returns:
    --------
    str
        Path to the generated PDF file

    Example:
    --------
    >>> conn = sqlite3.connect('../DoChaP-web/DB_merged.sqlite')
    >>> generate_gene_pdf('PUF60', conn)
    >>> conn.close()
    """
    if output_file is None:
        output_file = f"{gene_name}_visualization.pdf"

    # Create visualization using provided DB connection
    viz = GeneVisualization(conn, gene_name)
    viz.create_pdf(
        output_file,
        protein_only=protein_only,
        domains_only=domains_only,
    )
    return output_file


if __name__ == "__main__":
    # Example usage
    import argparse

    parser = argparse.ArgumentParser(description="Generate gene visualization PDF")
    parser.add_argument("gene_name", type=str, help="Gene symbol (e.g., PUF60)")
    parser.add_argument("-db", "--database", type=str, required=True,
                       help="Path to DoChap database")
    parser.add_argument("-o", "--output", type=str, default=None,
                       help="Output PDF file path")
    parser.add_argument("--protein-only", action="store_true",
                       help="Include only transcripts that produce protein")
    parser.add_argument("--domains-only", action="store_true",
                       help="Include only transcripts whose protein has domains")

    args = parser.parse_args()

    conn = sqlite3.connect(args.database)
    try:
        output = generate_gene_pdf(
            args.gene_name,
            conn,
            args.output,
            protein_only=args.protein_only,
            domains_only=args.domains_only,
        )
    finally:
        conn.close()
    print(f"Generated: {output}")
