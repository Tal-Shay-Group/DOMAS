import logging
import os

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Milestones worth putting on the console. Sits between INFO and WARNING so a
# console handler set to this level shows progress, warnings and errors, while
# the log file set to INFO keeps the per-chunk detail as well.
PROGRESS = 25
logging.addLevelName(PROGRESS, 'PROGRESS')


# rMATS-turbo per-event files, junction-count variant. JCEC carries identical
# coordinates, so JC is enough for DOMAS's coordinate mapping.
#
# An RI record names one junction, the spliced form, because the retaining
# isoform is defined by that junction's absence. It is emitted as two features at
# the same coordinates - one FEATURE_JUNCTION, one FEATURE_RETAINED_INTRON - so a
# transcript can hold a feature the canonical one lacks (_rmats_event_feature_types).
_RMATS_EVENT_FILES = {
    'SE': 'SE.MATS.JC.txt',
    'A5SS': 'A5SS.MATS.JC.txt',
    'A3SS': 'A3SS.MATS.JC.txt',
    'MXE': 'MXE.MATS.JC.txt',
    'RI': 'RI.MATS.JC.txt',
}


# ---------------------------------------------------------------------------
# Gene identity.
#
# DoChaP identifies a gene by gene_ensembl_id OR gene_GeneID_id (the NCBI side).
# 16,672 of 128,454 genes - 13% - carry only a GeneID, so a lookup keyed on the
# Ensembl id alone drops them and all of their transcripts. Absent ids are NULL,
# never an empty string.
#
# The junctions frame keeps one 'gene_ensembl_id' column holding whichever id the
# input supplied; these helpers accept either kind, mirroring the
# transcript_ensembl_id/transcript_refseq_id fallback used for transcripts.
# ---------------------------------------------------------------------------

GENE_ID_COLUMNS = ('gene_ensembl_id', 'gene_GeneID_id')


def combined_gene_ids(df):
    """The gene id to key on per row: gene_ensembl_id where present, else
    gene_GeneID_id - matching what the readers put in the junctions frame."""
    if 'gene_GeneID_id' not in df.columns:
        return df['gene_ensembl_id']
    return df['gene_ensembl_id'].fillna(df['gene_GeneID_id'])


def gene_id_clause(gene_ids):
    """(sql_fragment, params) matching `gene_ids` against either gene id column,
    e.g. "(gene_ensembl_id IN (?,?) OR gene_GeneID_id IN (?,?))" with the ids
    repeated once per column."""
    ids = [str(g) for g in gene_ids]
    placeholders = ','.join(['?'] * len(ids))
    fragment = ' OR '.join(f'{column} IN ({placeholders})' for column in GENE_ID_COLUMNS)
    return f'({fragment})', ids * len(GENE_ID_COLUMNS)


# ---------------------------------------------------------------------------
# The junctions frame: one contract for every reader.
#
# Each reader parses a different file and hands back the same frame, whose shape
# is declared here. Optional columns are filled with their default at the
# boundary (normalize_junctions_frame), so code past that point reads them
# directly.
# ---------------------------------------------------------------------------

REQUIRED_JUNCTION_COLUMNS = ('gene_ensembl_id', 'start_position', 'end_position', 'cluster_name')

def _optional_junction_defaults():
    """Optional columns and the value to fill them with.

    A function rather than a module constant so it can name the feature-type
    constants, which are defined further down with the matcher they belong to."""
    return {
        'chromosome': None,
        'gene_symbol': None,
        'event_type': None,      # the AS type the tool assigned (SE / A3SS / ...)
        'specie': None,          # 'human', 'mouse', ... - see _SPECIE_DB_NAME
        'junction_name': None,   # provenance only; nothing reads it
        FEATURE_TYPE_COLUMN: FEATURE_JUNCTION,
    }

# The species DOMAS supports: the label carried on the junctions frame mapped to
# the DoChaP Genes.specie value. Lives here so junction_analisys.py can check a
# stated species against the database without importing alternative_splicing.py.
SPECIE_DB_NAME = {
    'human': 'H_sapiens',
    'mouse': 'M_musculus',
    'rat': 'R_norvegicus',
    'zebrafish': 'D_rerio',
    'frog': 'X_tropicalis',
}

SPECIE_FROM_DB_NAME = {db_name: label for label, db_name in SPECIE_DB_NAME.items()}

# NCBI publishes no representative-transcript tag for zebrafish or frog, so
# DoChaP marks a canonical for them only where Ensembl supplies one. A gene with
# none falls back to its longest-CDS transcript (ClusterAnalysisResult._resolve_canonical),
# which is why both are supported: it affects 17.7% of comparable zebrafish genes
# and 6.1% of frog ones, against 0.4% for human.


# Ensembl stamps the species into its gene ids. The order is explicit so a
# species added later cannot shadow another by prefix.
_ENSEMBL_GENE_PREFIX_SPECIE = (
    ('ENSMUSG', 'mouse'),
    ('ENSRNOG', 'rat'),
    ('ENSDARG', 'zebrafish'),
    ('ENSXETG', 'frog'),
    ('ENSG', 'human'),
)


def specie_from_gene_id(gene_id):
    """The species an Ensembl gene id belongs to, or None if it is not one.

    rMATS, MAJIQ and SUPPA embed an Ensembl gene id per event but no species, so
    their frames carried no specie column. That is not cosmetic: clusters are
    grouped by (specie, cluster_name) where the column exists and by cluster_name
    alone where it does not, so a multi-species run of those formats would merge
    same-named clusters from different species into one.
    """
    if gene_id is None or (not isinstance(gene_id, str) and pd.isna(gene_id)):
        return None
    text = str(gene_id).strip().upper()
    for prefix, specie in _ENSEMBL_GENE_PREFIX_SPECIE:
        if text.startswith(prefix):
            return specie
    return None


def normalize_junctions_frame(df, specie=None):
    """Validate the required columns and fill every optional one with its default.

    `specie` is the species the caller states the input belongs to - required at
    the CLI, since three of the five formats carry no species field. It is
    authoritative: it fills the column, and any gene id that demonstrably
    contradicts it aborts the run.

    Contradiction is not the same as silence. An Ensembl id names its species in
    its prefix, so it can disagree; a GeneID or other non-Ensembl id names nothing
    and simply takes the stated value. Treating the second as an error would
    reject exactly the GeneID-only genes DOMAS goes out of its way to support.

    Where no species is stated the prefix is used, so library callers that predate
    the flag keep working.
    """
    missing = [c for c in REQUIRED_JUNCTION_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Columns {missing} are required in df_junctions but not found.")

    for column, default in _optional_junction_defaults().items():
        if column not in df.columns:
            df[column] = default

    derived = df['gene_ensembl_id'].map(specie_from_gene_id)
    if specie is None:
        df['specie'] = df['specie'].fillna(derived) if df['specie'].notna().any() else derived
        return df

    if specie not in SPECIE_DB_NAME:
        raise ValueError(f"Unknown specie {specie!r}. Expected one of: "
                         f"{', '.join(sorted(SPECIE_DB_NAME))}.")

    conflicting = df.loc[derived.notna() & (derived != specie), 'gene_ensembl_id']
    if not conflicting.empty:
        found = sorted(derived[conflicting.index].unique())
        raise ValueError(
            f"Input does not match -species {specie}: {len(conflicting)} of {len(df)} "
            f"rows carry gene ids from {', '.join(found)} "
            f"(e.g. {', '.join(str(g) for g in conflicting.unique()[:3])}). "
            f"Re-run with the species the data actually came from."
        )

    # Rows whose species could not be derived keep the stated one; rows that agree
    # are unaffected by writing it.
    df['specie'] = specie
    return df


def _clean_ensembl_id(gene_id):
    """rMATS GeneID like '"ENSG00000156256.15"' -> 'ENSG00000156256'."""
    gid = str(gene_id).strip().strip('"').strip("'")
    return gid.split('.')[0]


def _ordered_pair(a, b):
    """(a, b) as a genomically ordered (low, high) int pair."""
    a, b = int(a), int(b)
    return (a, b) if a <= b else (b, a)


def _rmats_event_junctions(event_type, r):
    """The (start, end) intron junctions implied by one rMATS event row `r`
    (a pandas Series indexable by the MATS column names). A junction is a
    (genomic_low, genomic_high) pair; JunctionsAnalysis's matcher resolves strand.

    No base adjustment is applied to the coordinates. Despite the `_0base` column
    names, these rMATS start values were verified to align *exactly* with DoChaP's
    stored exon coordinates (626/626 exact across all event types, zero
    off-by-one), so adding +1 would introduce a 1bp error rather than fix one. The
    matcher's 1bp tolerance is for exon- vs intron-boundary conventions, not to
    compensate for any base-offset here.

    SE and MXE name their exons by *genomic* position - rMATS reports
    upstreamEE <= exonStart_0base and exonEnd <= downstreamES on both strands
    (verified 100% over the fixture) - so one set of pairs is correct either way.

    A5SS and A3SS instead name the long/short forms of one exon plus its flanking
    exon, so which genomic boundary distinguishes the two forms flips with strand:
    the alternative donor of an A5SS event is the exon's genomic END on the plus
    strand but its genomic START on the minus strand, and the flanking exon sits
    on the other side. Reading the plus-strand boundary regardless of strand made
    both forms yield the *same* junction, collapsing the event to a single
    junction - which can never produce a comparable transcript, so every
    minus-strand A5SS/A3SS event was silently unanalysable."""
    strand = str(r['strand']) if 'strand' in r else '+'
    if event_type == 'SE':
        return [
            _ordered_pair(r['upstreamEE'], r['exonStart_0base']),  # inclusion: upstream -> exon
            _ordered_pair(r['exonEnd'], r['downstreamES']),        # inclusion: exon -> downstream
            _ordered_pair(r['upstreamEE'], r['downstreamES']),     # skipping
        ]
    if event_type in ('A5SS', 'A3SS'):
        # The two alternative-splice-site types exchange roles under strand
        # reversal: a minus-strand A5SS has the same geometry as a plus-strand
        # A3SS. Built only in this branch - MXE has no longExon*/short*/flanking*
        # columns, and a KeyError here is swallowed by rmats2junctions().
        varying_boundary_is_upper = (event_type == 'A5SS') == (strand != '-')
        if varying_boundary_is_upper:
            return [
                _ordered_pair(r['longExonEnd'], r['flankingES']),          # long form
                _ordered_pair(r['shortEE'], r['flankingES']),              # short form
            ]
        return [
            _ordered_pair(r['flankingEE'], r['longExonStart_0base']),      # long form
            _ordered_pair(r['flankingEE'], r['shortES']),                  # short form
        ]
    if event_type == 'MXE':
        return [
            _ordered_pair(r['upstreamEE'], r['1stExonStart_0base']),
            _ordered_pair(r['1stExonEnd'], r['downstreamES']),
            _ordered_pair(r['upstreamEE'], r['2ndExonStart_0base']),
            _ordered_pair(r['2ndExonEnd'], r['downstreamES']),
        ]
    if event_type == 'RI':
        # One interval, emitted twice: once as the junction the spliced isoform
        # carries, once as the intron the retained isoform contains. Same
        # coordinates - _rmats_event_feature_types() supplies the distinction.
        intron = _ordered_pair(r['upstreamEE'], r['downstreamES'])
        return [intron, intron]
    return []


# Columns _rmats_event_junctions() reads per event type, on top of the common
# GeneID/geneSymbol/chr/strand. Checked once per file so a renamed column in a
# future rMATS release fails on the header rather than row by row, where the
# per-row handler would swallow it into a silently truncated analysis.
_RMATS_REQUIRED_COLUMNS = {
    'SE': ('upstreamEE', 'exonStart_0base', 'exonEnd', 'downstreamES'),
    'A5SS': ('longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE',
             'flankingES', 'flankingEE'),
    'A3SS': ('longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE',
             'flankingES', 'flankingEE'),
    'MXE': ('upstreamEE', 'downstreamES', '1stExonStart_0base', '1stExonEnd',
            '2ndExonStart_0base', '2ndExonEnd'),
    'RI': ('upstreamEE', 'downstreamES'),
}

_RMATS_COMMON_COLUMNS = ('GeneID', 'geneSymbol', 'chr', 'strand')


def _check_rmats_columns(event_type, columns, path):
    """Raise if `path` lacks a column the event type needs."""
    required = _RMATS_COMMON_COLUMNS + _RMATS_REQUIRED_COLUMNS.get(event_type, ())
    missing = [c for c in required if c not in columns]
    if missing:
        raise ValueError(
            f"{os.path.basename(path)} is missing the column(s) {missing} that DOMAS "
            f"needs to build {event_type} junctions. Found: {sorted(columns)}. "
            f"This usually means the file is not rMATS-turbo {event_type} output, or "
            f"the format changed."
        )


def _rmats_event_feature_types(event_type):
    """Feature type per junction returned by _rmats_event_junctions(), in order.

    Only RI mixes types: its two entries share coordinates and are distinguished
    solely by type, which is what lets one event match both isoforms."""
    if event_type == 'RI':
        return [FEATURE_JUNCTION, FEATURE_RETAINED_INTRON]
    return None  # all plain junctions


def rmats2junctions(rmats_dir):
    """Parse an rMATS-turbo output directory (the five [Event].MATS.JC.txt files)
    into a junctions DataFrame ready for JunctionsAnalysis.analyze_junctions().

    Carries a FEATURE_TYPE_COLUMN: RI rows come in junction/retained-intron pairs,
    every other event type is plain junctions.

    rMATS provides the Ensembl GeneID and gene symbol directly, so no
    symbol->ensembl DB lookup is needed. ALL events are taken (no
    FDR/IncLevelDifference filtering) - pre-filter the input files if you only
    want significant events.
    """
    junctions = []
    idx = 0
    for event_type, filename in _RMATS_EVENT_FILES.items():
        path = os.path.join(rmats_dir, filename)
        if not os.path.exists(path):
            logger.warning(f"rMATS {event_type} file not found at {path} - no {event_type} "
                           f"events will be analysed.")
            continue
        df = pd.read_csv(path, sep='\t')
        _check_rmats_columns(event_type, df.columns, path)
        malformed = 0
        for _, r in df.iterrows():
            gene_ensembl_id = _clean_ensembl_id(r['GeneID'])
            gene_symbol = str(r['geneSymbol']).strip().strip('"')
            chromosome = str(r['chr'])
            if chromosome.startswith('chr'):
                chromosome = chromosome[3:]  # DoChaP stores bare '21', not 'chr21'
            idx += 1
            cluster_name = f'{event_type}_{gene_ensembl_id}_{chromosome}_{idx}'
            try:
                event_junctions = _rmats_event_junctions(event_type, r)
            except (ValueError, TypeError):
                # One row with an unparseable coordinate. The columns are known to
                # exist by now, so this is bad data, not a schema problem; it is
                # counted and reported below rather than passed over.
                malformed += 1
                continue
            feature_types = _rmats_event_feature_types(event_type)
            for position, (start, end) in enumerate(event_junctions):
                feature_type = FEATURE_JUNCTION if feature_types is None else feature_types[position]
                junctions.append([chromosome, gene_ensembl_id, gene_symbol,
                                  event_type, start, end, cluster_name, feature_type])

        if malformed:
            logger.warning(f"{os.path.basename(path)}: skipped {malformed} of {len(df)} "
                           f"{event_type} row(s) with unparseable coordinates.")

    if not junctions:
        raise ValueError(f"No rMATS MATS.JC.txt events found under {rmats_dir}")

    return pd.DataFrame(junctions, columns=[
        'chromosome', 'gene_ensembl_id', 'gene_symbol', 'event_type',
        'start_position', 'end_position', 'cluster_name', FEATURE_TYPE_COLUMN])


def _parse_coord_pairs(text):
    """'a-b;c-d' -> [(low, high), ...] as int pairs; blanks/malformed skipped."""
    pairs = []
    for token in str(text).split(';'):
        token = token.strip()
        if not token or '-' not in token:
            continue
        try:
            a, b = token.split('-')[:2]
            a, b = int(a), int(b)
        except ValueError:
            continue
        pairs.append((a, b) if a <= b else (b, a))
    return pairs


def voila2junctions(tsv_path):
    """Parse a MAJIQ voila TSV (`voila tsv` output, deltapsi or psi) into a
    junctions DataFrame ready for JunctionsAnalysis.analyze_junctions().

    One LSV -> one cluster; each coordinate in 'Junctions coords' (plus 'IR
    coords', when present) -> one junction row. voila embeds the Ensembl Gene ID,
    so - like the rMATS path - no symbol->ensembl lookup is needed. ALL LSVs are
    taken (no probability/dPSI filtering); pre-filter the TSV for significant LSVs
    if desired. event_type is built from MAJIQ's A5SS/A3SS/ES/IR classification.
    """
    df = pd.read_csv(tsv_path, sep='\t', dtype=str)  # header line starts with '#'

    gene_name_col = '#Gene Name' if '#Gene Name' in df.columns else 'Gene Name'

    junctions = []
    for _, r in df.iterrows():
        gene_ensembl_id = str(r['Gene ID']).split('.')[0].strip().strip('"')
        gene_symbol = str(r[gene_name_col]).strip()
        chromosome = str(r['chr']).strip()
        if chromosome.startswith('chr'):
            chromosome = chromosome[3:]  # DoChaP stores bare '5', not 'chr5'
        cluster_name = str(r['LSV ID']).strip()

        ir_coords = str(r['IR coords']).strip() if 'IR coords' in df.columns and pd.notna(r['IR coords']) else ''
        types = [name for name in ('A5SS', 'A3SS', 'ES')
                 if name in df.columns and str(r[name]).strip() == 'True']
        if ir_coords:
            types.append('IR')
        event_type = '+'.join(types) if types else 'LSV'

        # MAJIQ quantifies retention as one of the LSV's edges, so the retained
        # intron also appears in 'Junctions coords'; 'IR coords' names which edge
        # it is. That edge is emitted once, as a containment feature, and its
        # 'Junctions coords' copy dropped - the spliced junction for the same
        # intron is listed separately anyway, either at exactly (a-1, b+1) or
        # inside a larger junction spanning it (860 and 686 of 1,546 over the
        # fixture, none without either).
        retained_intron_pairs = _parse_coord_pairs(ir_coords) if ir_coords else []
        retained_intron_set = set(retained_intron_pairs)
        event_features = [(pair, FEATURE_JUNCTION)
                          for pair in _parse_coord_pairs(r['Junctions coords'])
                          if pair not in retained_intron_set]
        event_features += [(pair, FEATURE_RETAINED_INTRON) for pair in retained_intron_pairs]

        for (start, end), feature_type in event_features:
            junctions.append([chromosome, gene_ensembl_id, gene_symbol,
                              event_type, start, end, cluster_name, feature_type])

    if not junctions:
        raise ValueError(f"No LSV junctions found in {tsv_path}")

    return pd.DataFrame(junctions, columns=[
        'chromosome', 'gene_ensembl_id', 'gene_symbol', 'event_type',
        'start_position', 'end_position', 'cluster_name', FEATURE_TYPE_COLUMN])


def ioe2junctions(file_path):
    """
    Parses an IOE file and returns a DataFrame with the relevant data.

    SUPPA writes a retained intron as
    `<gene>;RI:<chr>:<s1>:<e1>-<s2>:<e2>:<strand>`, in which only one token holds
    a junction - the spliced form - because the retained isoform is defined by that
    junction's absence. Such an event is emitted as two features at the same
    coordinates, one FEATURE_JUNCTION and one FEATURE_RETAINED_INTRON, so each
    isoform matches the feature it actually carries. Emitted as a lone junction it
    was structurally unanalysable: 9,092 of 9,092 RI clusters in
    events_RI_strict.ioe held a single junction, while every other SUPPA event type
    held none.

    Args:
        file_path (str): The path to the IOE file.
    """
    junctions = []
    df_ioe = pd.read_csv(file_path, sep='\t')
    count = 0
    retained_intron_events = 0
    for row in df_ioe.itertuples():
        count += 1
        chromosome = row.seqname
        event_parts = row.event_id.split(';')
        gene_ensembl_id = event_parts[0]
        event_parts2 = event_parts[1].split(':')
        event_type = event_parts2[0]
        cluster_name = f'{event_type}_{gene_ensembl_id}_{chromosome}_{count}'
        
        current_junctions = []   
        for junction in event_parts2[1:]:
            if junction == '-':
                continue
            junction_parts = junction.split('-')
            if len(junction_parts) != 2:
                continue
            
            start = int(junction_parts[0])
            end = int(junction_parts[1])
            current_junctions.append((start, end))
        if event_type == 'RI':
            # One interval, two features: the junction the spliced isoform carries
            # and the intron the retained isoform contains.
            for start, end in current_junctions:
                junctions.append([chromosome, gene_ensembl_id, event_type, start, end,
                                  cluster_name, FEATURE_JUNCTION])
                junctions.append([chromosome, gene_ensembl_id, event_type, start, end,
                                  cluster_name, FEATURE_RETAINED_INTRON])
                retained_intron_events += 1
        elif event_type == 'SE':
            # SE:<chr>:<e1>-<s2>:<e2>-<s3> names the two inclusion junctions -
            # the skipped exon to each flanking exon - and leaves the skipping
            # junction (e1-s3) to be constructed. All three are emitted, matching
            # what rMATS gives for the same event; tests/test_cross_format.py
            # holds the two readers to that.
            if len(current_junctions) != 2:
                continue
            upstream_inclusion = current_junctions[0]
            downstream_inclusion = current_junctions[1]
            skipping = (current_junctions[0][0], current_junctions[1][1])
            for start, end in (upstream_inclusion, downstream_inclusion, skipping):
                junctions.append([chromosome, gene_ensembl_id, event_type, start, end,
                                  cluster_name, FEATURE_JUNCTION])
        else:
            for start, end in current_junctions:
                junctions.append([chromosome, gene_ensembl_id, event_type, start, end,
                                  cluster_name, FEATURE_JUNCTION])


    if retained_intron_events:
        logger.info(f"{os.path.basename(file_path)}: {retained_intron_events} RI event(s) "
                    f"emitted as junction + retained-intron feature pairs.")
    if not junctions:
        logger.warning(f"No analysable events found in {file_path}.")

    df_junctions = pd.DataFrame(junctions, columns=['chromosome', 'gene_ensembl_id', 'event_type',
                                                    'start_position', 'end_position', 'cluster_name',
                                                    FEATURE_TYPE_COLUMN])
    return df_junctions

def get_gene_symbols(con, gene_ensembl_ids):
    """
    Retrieves gene symbols for a list of gene IDs, each either an Ensembl gene id
    or a GeneID (see GENE_ID_COLUMNS).

    Args:
        gene_ensembl_ids (list): A list of gene IDs of either kind.
    Returns:
        dict: A dictionary mapping gene IDs to gene symbols.
    """
    clause, params = gene_id_clause(gene_ensembl_ids)
    query = f"SELECT gene_ensembl_id, gene_GeneID_id, gene_symbol FROM genes WHERE {clause}"
    df = pd.read_sql_query(query, con, params=params)
    gene_symbol_dict = dict(zip(combined_gene_ids(df), df['gene_symbol']))
    return gene_symbol_dict


def get_genes_number_of_transcripts(con, gene_ensembl_ids):
    """
    Retrieves the number of transcripts for a list of Ensembl gene IDs.

    Args:
        gene_ensembl_ids (list): A list of Ensembl gene IDs.
    Returns:
        dict: A dictionary mapping Ensembl gene IDs to the number of transcripts.
    """
    clause, params = gene_id_clause(gene_ensembl_ids)
    query = (f"SELECT gene_ensembl_id, gene_GeneID_id, "
             f"COUNT(COALESCE(transcript_ensembl_id, transcript_refseq_id)) AS num_transcripts "
             f"FROM transcripts WHERE {clause} "
             f"GROUP BY COALESCE(gene_ensembl_id, gene_GeneID_id)")
    df = pd.read_sql_query(query, con, params=params)
    num_transcripts_dict = dict(zip(combined_gene_ids(df), df['num_transcripts']))
    return num_transcripts_dict



def get_canonical_exon_counts(con, gene_ensembl_ids):
    """
    Given a list of gene ensembl ids, return a dict mapping each id to the
    number of exons (exon_count) in its canonical transcript.

    @param gene_ensembl_ids: list/iterable of gene_ensembl_id strings
    @param con: SQLite connection object
    @return: dict {gene_ensembl_id: exon_count}. Genes with no canonical
             transcript found in the db are mapped to None.
    """
    result = {gid: None for gid in gene_ensembl_ids}

    cur = con.cursor()
    clause, params = gene_id_clause(gene_ensembl_ids)
    query = f'''
        SELECT COALESCE(gene_ensembl_id, gene_GeneID_id), exon_count
            FROM Transcripts
            WHERE {clause}
              AND canonical != 0
        '''
    cur.execute(query, params)
    for gene_id, exon_count in cur.fetchall():
        result[gene_id] = exon_count

    return result


# ---------------------------------------------------------------------------
# Event feature <-> exon matching.
#
# Lives here so junction_analisys.py and generate_gene_pdf.py share one
# implementation: junction_analisys imports generate_gene_pdf, so the PDF layer
# cannot import back from it.
# ---------------------------------------------------------------------------

# An event is a list of features, each a (low, high) genomic coordinate pair
# carrying the type below. Both describe the same interval - an intron - but ask
# opposite questions of a transcript, which is what lets one event distinguish
# two isoforms:
#
#   FEATURE_JUNCTION       the intron is spliced OUT: two adjacent exons abut it
#   FEATURE_RETAINED_INTRON the intron is retained: one exon contains it
#
# A retained-intron event emits both, at identical coordinates. The type is what
# keeps them apart: matched type-blind, a spliced transcript would satisfy both
# and no transcript could hold a feature the canonical one lacks - exactly the
# collapse that makes intron retention unanalysable today.
FEATURE_JUNCTION = 'junction'
FEATURE_RETAINED_INTRON = 'retained_intron'

# Optional column on the junctions DataFrame carrying the type per row. Absent
# means every row is a plain junction, so a frame without the column is valid.
FEATURE_TYPE_COLUMN = 'feature_type'


def find_matching_junction_indices(df_transcript_exons, junctions, strand='+', feature_types=None):
    """
    Return the set of indices (into `junctions`) of event features that match
    this transcript's exon structure.

    A junction (start_position, end_position) matches the transcript if there
    are two exons, adjacent in transcript order, such that one exon's
    genomic_end_tx is within 1bp of the junction's intron-left boundary and the
    other exon's genomic_start_tx is within 1bp of the junction's intron-right
    boundary.

    On the positive strand the intron-left boundary is start_position and the
    intron-right boundary is end_position.  On the negative strand the
    transcript runs right-to-left in genomic coordinates, so the roles are
    reversed: the intron-left boundary (in genomic terms) is end_position and
    the intron-right boundary is start_position.

    `feature_types` is an optional sequence parallel to `junctions`, giving each
    feature's type (see FEATURE_JUNCTION / FEATURE_RETAINED_INTRON). None - the
    default - treats every feature as a junction, so existing callers and any
    junctions frame without the column behave exactly as before.

    A FEATURE_RETAINED_INTRON feature matches when a SINGLE exon contains the
    interval, within the same 1bp tolerance. Strand is irrelevant to it:
    containment is a genomic-coordinate test with no orientation.
    """
    if df_transcript_exons.empty or not junctions:
        return set()

    exon_starts = df_transcript_exons['genomic_start_tx'].to_numpy()
    exon_ends = df_transcript_exons['genomic_end_tx'].to_numpy()
    exon_orders = df_transcript_exons['order_in_transcript'].to_numpy()

    matched = set()
    for idx, (start_position, end_position) in enumerate(junctions):
        if feature_types is not None and feature_types[idx] == FEATURE_RETAINED_INTRON:
            # Retained: one exon spans the whole intron. Exons are disjoint, so
            # at most one can, and `any` is enough.
            low, high = min(start_position, end_position), max(start_position, end_position)
            if ((exon_starts <= low + 1) & (exon_ends >= high - 1)).any():
                matched.add(idx)
            continue
        if strand == '-':
            # Negative strand: transcript runs from high to low genomic coords.
            # DB always stores genomic_start_tx < genomic_end_tx, so:
            #   - the upstream exon's intron boundary is its genomic_start_tx (lower bound of the higher exon)
            #   - the downstream exon's intron boundary is its genomic_end_tx (upper bound of the lower exon)
            # Junction end_position is near the upstream exon's genomic_start_tx.
            # Junction start_position is near the downstream exon's genomic_end_tx.
            intron_left, intron_right = end_position, start_position
            upstream_orders = exon_orders[np.abs(exon_starts - intron_left) <= 1]
            downstream_orders = exon_orders[np.abs(exon_ends - intron_right) <= 1]
        else:
            intron_left, intron_right = start_position, end_position
            upstream_orders = exon_orders[np.abs(exon_ends - intron_left) <= 1]
            downstream_orders = exon_orders[np.abs(exon_starts - intron_right) <= 1]
        if len(upstream_orders) == 1 and len(downstream_orders) == 1:
            if abs(int(downstream_orders[0]) - int(upstream_orders[0])) == 1:
                matched.add(idx)
    return matched


# ---------------------------------------------------------------------------
# DoChaP database readers.
#
# Here rather than in alternative_splicing.py so junction_analisys.py can import
# them directly: alternative_splicing.py imports junction_analisys, so the
# dependency only runs one way.
# ---------------------------------------------------------------------------

def get_exons_for_transcripts(con, transcript_ids):
    # 450 ids per chunk: each is bound twice, and SQLite allows 999 parameters.
    chunk_size = 450
    df_list = []
    
    t_ids = list(transcript_ids)  # Ensure it's a list if it's a different iterable
    for i in range(0, len(t_ids), chunk_size):
        chunk_ids = t_ids[i:i + chunk_size]
        
        placeholders = ','.join(['?'] * len(chunk_ids))
        
        query = f'''
            SELECT * FROM Transcript_exon 
            WHERE transcript_ensembl_id IN ({placeholders}) 
            OR transcript_refseq_id IN ({placeholders})
        '''
        
        # Once per IN clause.
        params = chunk_ids * 2
        
        df_chunk = pd.read_sql_query(query, con, params=params)
        df_list.append(df_chunk)

    if not df_list:
        # No transcript ids, so the chunk loop never ran and pd.concat([]) would
        # raise. Happens whenever no event resolves a gene - every LeafCutter
        # cluster unannotated, or a gene set absent from the database. Return the
        # table's empty shape so the run finishes and reports each event's reason.
        return pd.read_sql_query('SELECT * FROM Transcript_exon LIMIT 0', con)

    df_exons = pd.concat(df_list, ignore_index=True)
    return df_exons


def get_genes_df_transcripts(con, gene_ids):
    """
    Load all transcripts (with their canonical flag) for a list of genes.

    Args:
        con: Database connection
        gene_ids: List of gene IDs, each either an Ensembl gene id or a GeneID

    Returns:
        DataFrame with one row per transcript, including transcript_ensembl_id,
        transcript_refseq_id, gene_ensembl_id, gene_GeneID_id and canonical columns.
    """
    dfs = []
    # 450 per batch: gene_id_clause() binds each id once per gene id column, so
    # 900 parameters, under SQLite's 999 limit.
    batch_size = 450
    for i in range(0, len(gene_ids), batch_size):
        batch = gene_ids[i:i + batch_size]
        clause, params = gene_id_clause(batch)
        query = f'SELECT * FROM Transcripts WHERE {clause}'
        dfs.append(pd.read_sql_query(query, con, params=params))
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def _link_proteins_by_refseq(df_transcript, df_protein, df_all_proteins):
    """Recover the protein of a transcript whose protein_ensembl_id names no
    Proteins row, by falling back to its protein_refseq_id.

    Proteins holds one row per protein, with UNIQUE protein_refseq_id and
    protein_ensembl_id, so it can only back-link to one transcript. Where two
    transcripts encode the same protein - the pseudoautosomal genes, whose X and Y
    copies are identical - the second transcript's protein_ensembl_id matches
    nothing and its domains are unreachable, even though the protein and its
    domains are in the database under the sibling. 35 coding transcripts across 12
    PAR genes (SHOX, PLCXD1, CRLF2, CSF2RA, ...) are in that state; each one's
    protein IS found by its protein_refseq_id.

    The transcript's own protein_ensembl_id is the broken half, so it is rewritten
    to the one the Proteins row carries - that is the id DomainEvent and
    RepresentativeDomains are keyed on, and the point of the exercise is to reach
    those domains.
    """
    linked = set(df_protein['transcript_ensembl_id'])
    unlinked = df_transcript[~df_transcript['transcript_ensembl_id'].isin(linked)]
    if 'protein_refseq_id' not in unlinked.columns or unlinked.empty:
        return df_transcript, df_protein
    refseq = unlinked['protein_refseq_id'].astype(str).str.strip()
    unlinked = unlinked[refseq.notna() & ~refseq.isin(('', 'nan', 'None'))]
    if unlinked.empty:
        return df_transcript, df_protein

    by_refseq = df_all_proteins.dropna(subset=['protein_refseq_id'])
    by_refseq = by_refseq[by_refseq['protein_ensembl_id'].notna()]
    by_refseq = by_refseq.drop_duplicates('protein_refseq_id').set_index('protein_refseq_id')

    recovered, rewritten = [], {}
    for _, transcript in unlinked.iterrows():
        key = str(transcript['protein_refseq_id']).strip()
        if key not in by_refseq.index:
            continue
        protein = by_refseq.loc[key].copy()
        protein['protein_refseq_id'] = key
        # The record is real; only its back-link names the sibling transcript.
        protein['transcript_ensembl_id'] = transcript['transcript_ensembl_id']
        protein['transcript_refseq_id'] = transcript['transcript_refseq_id']
        recovered.append(protein)
        rewritten[transcript['transcript_ensembl_id']] = protein['protein_ensembl_id']

    if not recovered:
        return df_transcript, df_protein

    df_protein = pd.concat([df_protein, pd.DataFrame(recovered)], ignore_index=True)
    df_transcript = df_transcript.copy()
    mask = df_transcript['transcript_ensembl_id'].isin(rewritten)
    df_transcript.loc[mask, 'protein_ensembl_id'] = (
        df_transcript.loc[mask, 'transcript_ensembl_id'].map(rewritten))
    logger.log(PROGRESS, 'Linked %d transcript(s) to their protein by protein_refseq_id', len(rewritten))
    return df_transcript, df_protein


def _read_transcripts_and_proteins(con, transcript_ids):
    """Transcripts/Proteins rows for transcript_ids, filtered to a non-empty
    protein_ensembl_id. Shared by get_transcript_domains_db() and
    get_representative_domains_db() so get_domains_db() reads these tables once."""
    df_transcript = pd.read_sql_query('select * from Transcripts', con)
    df_transcript = df_transcript[df_transcript.transcript_ensembl_id.isin(transcript_ids)]
    df_all_proteins = pd.read_sql_query('select * from Proteins', con)
    df_protein = df_all_proteins[df_all_proteins.transcript_ensembl_id.isin(transcript_ids)]
    df_protein = df_protein.dropna(subset=['protein_ensembl_id'])
    df_protein = df_protein[df_protein.protein_ensembl_id.str.strip() != '']
    return _link_proteins_by_refseq(df_transcript, df_protein, df_all_proteins)


def get_transcript_domains_db(con, transcript_ids, df_transcript=None, df_protein=None):
    logger.log(PROGRESS, 'Reading domains from DomainEvent/DomainType')
    if df_transcript is None or df_protein is None:
        df_transcript, df_protein = _read_transcripts_and_proteins(con, transcript_ids)
    proteins_ids = np.unique(df_protein.protein_ensembl_id.values).tolist()
    df_domain_event = pd.read_sql_query('select * from DomainEvent', con)

    df_domain_event = df_domain_event[df_domain_event.protein_ensembl_id.isin(proteins_ids)]
    df_domain_event = df_domain_event.dropna(subset=['protein_ensembl_id'])
    df_domain_event = df_domain_event[df_domain_event.protein_ensembl_id.str.strip() != '']
    df_domain_type = pd.read_sql_query('select * from DomainType', con)
    type_ids = np.unique(df_domain_event.type_id.values).tolist()
    df_domain_type = df_domain_type[df_domain_type.type_id.isin(type_ids)]

    merged_df = pd.merge(df_protein, df_transcript, on=['protein_ensembl_id', 'transcript_ensembl_id'])
    merged_df = merged_df.drop(columns=['gene_GeneID_id', 'synonyms'])
    merged_df = pd.merge(merged_df, df_domain_event, on='protein_ensembl_id')
    merged_df = merged_df.drop(columns=['protein_refseq_id_x', 'length',
                               'protein_refseq_id_y', 'nuc_start','nuc_end',
                               'total_length','splice_junction', 'complete_exon'])
    merged_df = pd.merge(merged_df, df_domain_type, on='type_id')
    merged_df = merged_df.dropna(subset=['AA_start', 'AA_end'])
    merged_df = merged_df.astype(str)
    merged_df = merged_df.fillna('nan')
    merged_df['AA_start'] = merged_df['AA_start'].astype(float).astype(int)
    merged_df['AA_end'] = merged_df['AA_end'].astype(float).astype(int)

    merged_df = merged_df.rename(columns={'protein_ensembl_id': 'protein_ensembl_id_version',
                                          'transcript_ensembl_id': 'transcript_ensembl_id_version',
                                          'description_y': 'short_description',
                                          'gene_ensembl_id' : 'gene_ensembl_id'})
    # DomainType.description is what this path already carries as short_description
    # (description_y above; description_x is the protein's). Expose it under
    # `description` too, so a domain frame from either source answers to the same
    # column name and the results CSV does not care which source it came from.
    merged_df['description'] = merged_df['short_description']
    merged_df = merged_df.drop(columns=['type_id', 'ext_id', 'name', 'other_name', 'description_x'])
    merged_df = merged_df.drop(columns=['transcript_refseq_id_x', 'tx_start', 'tx_end', 'cds_start', 'cds_end', 'exon_count'])
    merged_df = merged_df.drop(columns=['transcript_refseq_id_y','protein_refseq_id'])
    logger.log(PROGRESS, 'Read %d domain rows from DomainEvent/DomainType', len(merged_df))
    return merged_df


REPRESENTATIVE_DOMAINS_COLUMNS = [
    'protein_ensembl_id_version', 'transcript_ensembl_id_version', 'protein_interpro_id',
    'gene_ensembl_id', 'canonical', 'AA_start', 'AA_end', 'short_description',
    # The entry's prose description, reported per domain in the results CSV.
    # RepresentativeDomains carries its own; the DomainEvent/DomainType path
    # supplies DomainType.description under the same name, so a frame from
    # either source answers to 'description'.
    'description',
    'CDD_id', 'cdd', 'pfam', 'smart', 'tigr', 'interpro',
    # domain_id and the InterPro entry `type` are carried through so
    # junction_analisys.filter_representative_domains() can rank the domain set
    # by curated type (Domain/Repeat above Family/Homologous_superfamily). `type`
    # is NULL in DBs whose RepresentativeDomains lacks the column; the filter
    # then passes the frame through unchanged.
    'domain_id', 'type',
]


def _route_domain_id_to_column(domain_id):
    """Map an InterPro match-XML accession (RepresentativeDomains.domain_id) to the
    DomainType-style identifier column it belongs to. Only affects which bucket the id is
    displayed under - domain matching in junction_analisys.py checks all buckets together."""
    if domain_id.startswith('IPR'):
        return 'interpro'
    if domain_id.startswith('PF'):
        return 'pfam'
    if domain_id.startswith('SM'):
        return 'smart'
    if domain_id.startswith('TIGR'):
        return 'tigr'
    if domain_id.lower().startswith('cd'):
        return 'CDD_id'
    return 'interpro'


def get_representative_domains_db(con, transcript_ids, df_transcript=None, df_protein=None):
    """
    Domains sourced from the RepresentativeDomains table (populated by
    DoChaP-db/InterProRepresentativeDomains.py) instead of DomainEvent/DomainType.

    Returns a DataFrame with the same columns as get_transcript_domains_db() so it's a
    drop-in replacement for build_domain_lookup(), but only contains rows for proteins
    that actually have a RepresentativeDomains entry - see get_domains_db() for combining
    it with the DomainEvent/DomainType fallback for proteins that don't.
    """
    logger.log(PROGRESS, 'Reading domains from RepresentativeDomains')
    if df_transcript is None or df_protein is None:
        df_transcript, df_protein = _read_transcripts_and_proteins(con, transcript_ids)

    if 'protein_interpro_id' not in df_protein.columns:
        logger.warning('Proteins table has no protein_interpro_id column '
                       '(InterProRepresentativeDomains.py has not been run against this DB).')
        return pd.DataFrame(columns=REPRESENTATIVE_DOMAINS_COLUMNS)

    df_protein = df_protein.dropna(subset=['protein_interpro_id'])
    df_protein = df_protein[df_protein.protein_interpro_id.str.strip() != '']
    if df_protein.empty:
        return pd.DataFrame(columns=REPRESENTATIVE_DOMAINS_COLUMNS)

    interpro_ids = np.unique(df_protein.protein_interpro_id.values).tolist()
    try:
        df_rep = pd.read_sql_query('select * from RepresentativeDomains', con)
    except (sqlite3.OperationalError, pd.errors.DatabaseError):
        logger.warning('RepresentativeDomains table not found in this DB.')
        return pd.DataFrame(columns=REPRESENTATIVE_DOMAINS_COLUMNS)

    df_rep = df_rep[df_rep.protein_interpro_id.isin(interpro_ids)]
    df_rep = df_rep.dropna(subset=['start', 'end'])
    if df_rep.empty:
        return pd.DataFrame(columns=REPRESENTATIVE_DOMAINS_COLUMNS)

    # Some DBs have no `type` column; add it as NULL to keep the frame uniform.
    if 'type' not in df_rep.columns:
        df_rep['type'] = None

    merged_df = pd.merge(df_protein, df_transcript, on=['protein_ensembl_id', 'transcript_ensembl_id'])
    merged_df = pd.merge(merged_df, df_rep, on='protein_interpro_id')

    domain_column = merged_df['domain_id'].map(_route_domain_id_to_column)
    for col in ('CDD_id', 'cdd', 'pfam', 'smart', 'tigr', 'interpro'):
        merged_df[col] = merged_df['domain_id'].where(domain_column == col)

    merged_df = merged_df.rename(columns={
        'protein_ensembl_id': 'protein_ensembl_id_version',
        'transcript_ensembl_id': 'transcript_ensembl_id_version',
        'start': 'AA_start',
        'end': 'AA_end',
        'domain_name': 'short_description',
        # Proteins and RepresentativeDomains both have a `description`, so the
        # merge suffixes them: _x is the protein's, _y the domain entry's. It is
        # the domain's that belongs on a domain row.
        'description_y': 'description',
    })
    merged_df['AA_start'] = merged_df['AA_start'].astype(int)
    merged_df['AA_end'] = merged_df['AA_end'].astype(int)

    logger.log(PROGRESS, 'Read %d domain rows from RepresentativeDomains', len(merged_df))
    return merged_df[REPRESENTATIVE_DOMAINS_COLUMNS]


def get_domains_db(con, transcript_ids, use_representative_domains=False):
    """
    Domain source used by JunctionsAnalysis.analyze_junctions().

    use_representative_domains=False (default): unchanged - domains come from
    DomainEvent/DomainType exactly as get_transcript_domains_db() always has, so the
    existing algorithm runs without any change.

    use_representative_domains=True: domains come from RepresentativeDomains and
    nowhere else. A protein with no entry there has no domains, rather than falling
    back to DomainEvent/DomainType: the two sources carry different coordinates and
    different notions of what a domain is, so mixing them within one run makes a
    comparison's result depend on which source each side happened to be drawn from.
    A protein that drops out is visible as no_domains_in_region, not as a silent
    substitution.
    """
    df_transcript, df_protein = _read_transcripts_and_proteins(con, transcript_ids)
    if use_representative_domains:
        return get_representative_domains_db(con, transcript_ids, df_transcript=df_transcript,
                                             df_protein=df_protein)
    return get_transcript_domains_db(con, transcript_ids, df_transcript=df_transcript,
                                     df_protein=df_protein)
