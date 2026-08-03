import logging
import math
import os

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# rMATS-turbo per-event result files (junction-count variant). The JCEC files
# carry identical coordinates, so JC is enough for DOMAS's coordinate mapping.
#
# RI is read again. A retained intron has only one junction - the spliced form -
# because the retained isoform is defined by that junction's ABSENCE, so the event
# used to collapse to a single junction and no transcript could hold a feature the
# canonical one lacked. It is now emitted as two features at the same coordinates,
# one FEATURE_JUNCTION and one FEATURE_RETAINED_INTRON, which the matcher tells
# apart (see _rmats_event_feature_types()).
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
# DoChaP identifies a gene by gene_ensembl_id OR gene_GeneID_id (the NCBI/RefSeq
# side). 16,672 of 128,454 genes - 13% - carry ONLY a GeneID, so a lookup keyed on
# gene_ensembl_id alone silently drops them along with all of their transcripts.
# Absent ids are stored as NULL, never as an empty string. (The share was highest
# in D_rerio and X_tropicalis, but it is still 4,515 genes across human, mouse
# and rat.)
#
# The junctions frame keeps one 'gene_ensembl_id' column holding whichever id the
# input supplied; these helpers let a lookup accept either kind, mirroring the
# transcript_ensembl_id/transcript_refseq_id fallback already used for transcripts.
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
# Each reader parses a different file, but they all hand back the same frame, and
# this is where its shape is declared. Previously each reader carried its own
# column list and the differences were absorbed downstream by a dozen scattered
# "if 'x' in df.columns" defaults - which is how three of the five formats ended
# up with no specie column at all (see specie_from_gene_id below).
#
# Optional columns are filled with their default at the boundary
# (normalize_junctions_frame), so code past that point can read them directly.
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

# The species DOMAS supports, mapping the label carried on the junctions frame to
# the DoChaP Genes.specie value. Lives here rather than in alternative_splicing.py
# so that junction_analisys.py can cross-check the stated species against the
# database without importing it (which would be circular).
SPECIE_DB_NAME = {
    'human': 'H_sapiens',
    'mouse': 'M_musculus',
    'rat': 'R_norvegicus',
    'zebrafish': 'D_rerio',
    'frog': 'X_tropicalis',
}

SPECIE_FROM_DB_NAME = {db_name: label for label, db_name in SPECIE_DB_NAME.items()}

# zebrafish and frog were dropped for a while, because NCBI publishes no
# representative-transcript tag (MANE Select or RefSeq Select) for them, so
# DoChaP can mark a canonical only where Ensembl supplies one - never for a gene
# held on the RefSeq side alone. That is no longer disqualifying: a gene with no
# flagged canonical now falls back to its longest-CDS transcript as a stand-in
# (see ClusterAnalysisResult.analyze_cluster), so such a gene is analysed rather
# than abandoned.


# Ensembl stamps the species into its gene ids. Longest-prefix-first is not
# needed - 'ENSMUSG...' does not start with 'ENSG' - but the order is kept
# explicit so adding a species cannot silently shadow another.
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
            f"Input does not match -specie {specie}: {len(conflicting)} of {len(df)} "
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


# Columns _rmats_event_junctions() reads, per event type, on top of the common
# GeneID/geneSymbol/chr/strand. Checked once per file so that a schema change -
# a renamed column in a future rMATS release - fails loudly on the header instead
# of raising KeyError on every row and being swallowed row by row, which would
# silently produce an empty or truncated analysis.
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
                # A single row with an unparseable coordinate. Counted and reported
                # below rather than passed over in silence: the columns are known to
                # exist by now, so this is bad data in one row, not a schema problem.
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

        # MAJIQ quantifies intron retention as one of the LSV's edges, so the
        # retained intron appears in 'Junctions coords' too - and 'IR coords'
        # names which edge it is. That edge is NOT a splice junction: the actual
        # junction for the same intron is listed separately, either at exactly
        # (a-1, b+1) - the same intron in the flanking-exon-boundary convention -
        # or inside a larger junction spanning it (verified over the whole
        # fixture: 860 and 686 of 1546, none without one).
        #
        # So the IR interval is emitted ONCE, as a containment feature. Emitting
        # it as a junction as well asked the opposite question of it: matched by
        # adjacency it finds the transcript that splices the intron out, never the
        # one that retains it, and it merely duplicated the real junction that the
        # 1bp tolerance already matches.
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
            # SE:<chr>:<e1>-<s2>:<e2>-<s3> names the two INCLUSION junctions - the
            # skipped exon to each flanking exon - and leaves the skipping junction
            # (e1-s3) to be constructed. All three are emitted.
            #
            # The second inclusion junction used to be dropped: it supplied only its
            # right coordinate to the constructed skipping junction and was never
            # added in its own right, so SUPPA produced 2 junctions where rMATS
            # produces 3 for the same event (verified across all 51,450 SE clusters
            # of events_SE_strict.ioe). A transcript carrying only that downstream
            # junction - an alternative first exon starting inside the skipped exon -
            # was therefore comparable when the event came from rMATS but not from
            # SUPPA. See tests/test_cross_format.py.
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


# ============================================================================
# Domain-change scoring (analysis-flow functions; kept out of the PDF code).
#
# These are pure functions so the analysis can compute a per-comparison
# `impact` and a per-isoform `spade_score`, write them as columns in the result
# table, and let the GUI and the PDF read them instead of recomputing.
# ============================================================================

def hmm_change_impact(canonical_cov, alt_cov, canonical_plddt,
                      fold_state=None, hits_functional_site=False, region_am=None,
                      region_rsa=None):
    """Deterministic likely-impact label for a change to one HMM (Pfam) element.

    Base signal: how much of the Pfam model the alternative isoform loses
    (coverage), weighted by whether the changed region is structured.

        coverage lost (pts)      base level
          lost entirely / >= 30    high
          10 - 30                  moderate
          < 10                     low
        structure weighting
          structured -> up one level (max high)
          disordered -> down one level (min low)

    "Structured" is read best-signal-first: UniProt `fold_state` (curated) >
    AlphaFold **burial** of the changed region (`region_rsa`, mean relative
    solvent accessibility) > AlphaFold `canonical_plddt` (confidence). Burial
    beats pLDDT because it is region-specific and directional: a buried-core loss
    destabilises the fold and (empirically) is where pathogenic variants sit, so
    it raises impact; a solvent-exposed surface loss is usually tolerated, so it
    lowers it. Thresholds: region_rsa < 0.30 -> buried/structured; > 0.50 ->
    exposed. (fold_state='folded' or pLDDT >= 70 -> structured;
    fold_state='disordered' or pLDDT < 50 -> disordered.)

    hits_functional_site: the changed region overlaps a UniProt functional
    residue (ACT_SITE/BINDING/MOD_RES/SITE/MOTIF/...). Losing such a residue is
    significant regardless of how much coverage/structure changed, so it forces
    at least 'high' when coverage is lost, and 'moderate' even when the domain
    coverage is unchanged (the disordered-functional-motif case pLDDT alone
    would dismiss).

    Returns: 'none', 'gain', 'low', 'moderate', 'high', or 'n/a'.
    """
    if canonical_cov is None:
        return 'gain' if alt_cov is not None else 'n/a'
    a = 0 if alt_cov is None else alt_cov
    loss = canonical_cov - a
    if loss <= 0:
        # no domain coverage lost; still flag if a functional residue is affected
        # or the changed region is functionally constrained (high AlphaMissense).
        if hits_functional_site or (region_am is not None and region_am >= 0.564):
            return 'moderate'
        return 'none'
    if hits_functional_site:
        return 'high'
    if a == 0 or loss >= 30:
        level = 3
    elif loss >= 10:
        level = 2
    else:
        level = 1
    # constraint/structure weighting: take the FIRST signal that has an opinion,
    # in priority order AlphaMissense constraint > curated fold_state > AlphaFold
    # burial > pLDDT confidence. A signal that is present but *inconclusive*
    # (e.g. mid-range AlphaMissense 0.34-0.564) yields to the next, so burial can
    # still nudge when AM is ambiguous. (region_am >= 0.564 = AM "likely
    # pathogenic"; <= 0.34 = "likely benign"; region_rsa < 0.30 = buried core;
    # > 0.50 = solvent-exposed; pLDDT >= 70 structured, < 50 disordered.)
    ladder = []
    if region_am is not None:
        ladder.append((region_am >= 0.564, region_am <= 0.34))
    if fold_state == 'folded':
        ladder.append((True, False))
    elif fold_state == 'disordered':
        ladder.append((False, True))
    if region_rsa is not None:
        ladder.append((region_rsa < 0.30, region_rsa > 0.50))
    if canonical_plddt is not None:
        ladder.append((canonical_plddt >= 70, canonical_plddt < 50))
    raise_lvl, lower_lvl = False, False
    for r, l in ladder:
        if r or l:
            raise_lvl, lower_lvl = r, l
            break
    if raise_lvl and level < 3:
        level += 1
    elif lower_lvl and level > 1:
        level -= 1
    return {1: 'low', 2: 'moderate', 3: 'high'}[level]


def insertion_impact(inserted_len, inside_domain, junction_fold_state=None):
    """Deterministic impact label for residues present ONLY in the alternative
    isoform (an insertion, i.e. the added part of a longer isoform).

    hmm_change_impact is loss-based: an insertion removes no canonical Pfam
    coverage, so it scores 'none' there even for a large added segment. An
    insertion is scored instead on where the added residues land - because
    insertions are usually tolerated: empirically only ~1 in 5 isoform insertions
    changes the fold (most sit in exposed termini/linkers). So SIZE ALONE does not
    imply disruption; the escalating factor is whether the residues land inside a
    folded/Pfam domain.

        placement                              base level
          outside any domain (terminus/linker)
            >= 30 aa                             low
            < 30 aa                              none
          inside a Pfam domain
            >= 30 aa                             moderate
            < 30 aa                              low
        junction fold-state
          'folded'      -> up one level (max high)   [inserted into ordered structure]
          'disordered'  -> down one level (min none) [inserted into a tolerant region]

    inside_domain: the inserted residues overlap a Pfam domain on the
    alternative sequence. An insertion within a structured domain disrupts the
    fold; one at a terminus or in a linker usually does not.

    junction_fold_state: UniProt fold-state ('folded'/'disordered'/None) of the
    canonical residues flanking the insertion point.

    Returns: 'none', 'low', 'moderate', or 'high'.
    """
    if inserted_len <= 0:
        return 'none'
    if inside_domain:
        level = 2 if inserted_len >= 30 else 1
    else:
        level = 1 if inserted_len >= 30 else 0
    if junction_fold_state == 'folded' and level < 3:
        level += 1
    elif junction_fold_state == 'disordered' and level > 0:
        level -= 1
    return {0: 'none', 1: 'low', 2: 'moderate', 3: 'high'}[level]


# Calibrated logistic model for the continuous impact probability. Trained on
# 3,282 isoform changed-regions labelled by UniProt humsavar variant overlap
# (pathogenic LP/P vs benign LB/B), features standardised, missing values imputed
# to the training median. Gene-grouped 5-fold CV AUC 0.766 (every isoform of a gene
# confined to one fold, so no protein is in both train and test; the random 80/20
# hold-out gives 0.778 - the gap is the leakage check, not the headline) and
# well-calibrated (predicted prob tracks observed pathogenic rate). Read it as an
# enrichment/prioritisation signal, not a per-event classifier.
#
# HOW MUCH IS THIS MODEL WORTH OVER ITS DOMINANT FEATURE (E45, measured): ranking the
# same regions by raw region_am with no model at all gives AUC 0.753, against the
# model's 0.766 - delta +0.013, gene-clustered bootstrap 95% CI [+0.001, +0.024]. Real
# but marginal: impact_prob is region-averaged AlphaMissense with a small correction.
# Nothing leaks (AlphaMissense is not trained on clinical labels; it uses population
# frequency as weak supervision) but its thresholds were calibrated on ClinVar, which
# humsavar aggregates, so the 0.766 is largely AlphaMissense's own performance.
# The AlphaMissense-FREE combination [loeuf, buried_frac, max_cov_loss] reaches 0.715 -
# that is the number to quote for signal DOMAS contributes independently. Note also
# max_cov_loss alone = 0.523, i.e. chance: Pfam coverage change, the one thing DOMAS
# measures natively, carries no functional signal against this label.
# Reproduce with alphafold_benchmark/am_alone.py.
# Coefficients are on standardised features; region_am dominates (+0.80), gnomAD
# LOEUF is the non-circular add (-0.42, lower LOEUF = more constrained = more
# pathogenic), burial +0.24, coverage-loss ~0.
# Regenerate with alphafold_benchmark/fit_calibrated.py.
_IMPACT_PROB_MODEL = {
    'features': ['region_am', 'loeuf', 'max_cov_loss', 'buried_frac'],
    'median':   [0.482, 0.983, 41.0, 0.3402],
    'mean':     [0.482, 0.9788, 46.071, 0.3091],
    'std':      [0.1554, 0.4429, 41.0077, 0.2193],
    'coef':     [0.7983, -0.423, -0.069, 0.2363],
    'intercept': -1.19,
}


def impact_probability(region_am=None, loeuf=None, max_cov_loss=None, buried_frac=None):
    """Calibrated probability (0-1) that the region the alternative isoform changes
    contains a PATHOGENIC missense variant rather than only benign ones - estimated
    from the constraint on those residues (mean AlphaMissense, gnomAD LOEUF) together
    with their burial in the canonical fold and the Pfam coverage lost. The functional
    companion to fold_change_probability (which is structural), and the *continuous*
    companion to the categorical `hmm_change_impact`: it does not quantise, and it
    folds in the gene-level constraint.

    Label - the operational definition (alphafold_benchmark/fit_calibrated.py). For a
    canonical changed span [lo, hi], over UniProt humsavar single-residue variants:
      pat  - an LP/P (pathogenic) position falls inside the span
      ben  - an LB/B (benign) position falls inside the span
      kept - only rows with (pat or ben): a span carrying no classified variant is
             not a training example at all
      y    - int(pat), so a span carrying BOTH classes counts as positive
    Defined for substitutions and deletions; pure insertions are left blank (they
    have no canonical changed region, and were dropped from training too).

    Conditional in training, unconditional at runtime. Every training row already
    contained a classified variant - the base rate is P(pathogenic | overlaps a
    variant) = 0.277. DOMAS then scores EVERY changed region, including ones with no
    variant evidence at all, so read the output as "how much does this region resemble
    the pathogenic-variant-carrying ones", not as a literal probability for an
    arbitrary region.

    What it is NOT: a verdict on the splice event. The label is variant-overlap of a
    SPAN - whatever survives prefix/suffix trimming in `changed_region`, so its length
    varies a lot - and one feature (loeuf) is gene-level, shared by every comparison of
    that gene. The score says "this region is the kind of place pathogenic variants
    sit", not "this isoform is pathogenic": a prioritisation signal (gene-grouped AUC
    0.766), not a classifier. Nothing in it involves the isoform's fold, expression or
    consequence; the fold is fold_change_probability, and the two can disagree in both
    directions.

    Features (see `_IMPACT_PROB_MODEL`) - note only two are region-specific:
      region_am    - mean AlphaMissense pathogenicity over the changed residues.
                     Dominant (+0.80), but partly circular: AlphaMissense is itself
                     calibrated against clinical variant sets overlapping humsavar.
      buried_frac  - fraction of the changed residues buried in the canonical
                     AlphaFold model (RSA < 0.25).
      loeuf        - gnomAD LoF constraint of the GENE, not the region, so every
                     comparison of a gene shares this term (-0.42, lower = more
                     constrained). The non-circular signal.
      max_cov_loss - max Pfam coverage (%) lost. Coefficient ~0 (-0.07): domain loss
                     carries essentially no functional signal here, unlike in
                     fold_change_probability where it is a real term.
    A missing feature is imputed to the training median (a neutral value), so the
    score degrades gracefully. See `_IMPACT_PROB_MODEL` for provenance and
    calibration.

    Returns a float in [0, 1], or None if every feature is missing.
    """
    m = _IMPACT_PROB_MODEL
    vals = [region_am, loeuf, max_cov_loss, buried_frac]
    if all(v is None for v in vals):
        return None
    z = m['intercept']
    for v, med, mean, std, coef in zip(vals, m['median'], m['mean'], m['std'], m['coef']):
        x = med if v is None else v
        z += coef * ((x - mean) / std)
    return round(1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, z)))), 4)


# Calibrated logistic model for the fold-change probability (structural axis):
# P(the alternative isoform has a different fold, TM-score < 0.5). Labels are the
# Genome Biology 2025 (Song et al.) AlphaFold TM-scores; features standardised,
# missing imputed to the training median.
#
# Fit on the 6,973 pairs whose ISOFORM structure the source paper classed `high` or
# `confident` - the same quality filter the paper applies to its own metric analyses
# (E44). The earlier fit used all 10,227 pairs, but 32% of those are `low` /
# `unstructured` rows whose TM<0.5 rate is 0.74 vs 0.15 in the confident ones: isoform
# pLDDT correlates with TM at rho +0.755, so those rows largely label AlphaFold's
# failure to model the isoform rather than a splicing-induced fold change. Including
# them inflated the intercept (-1.26 -> -2.51 once removed) and broke calibration on
# exactly the well-predicted proteins users care about: the old model's top quintile
# there predicted 0.68 against an observed 0.46, and its 'change' calls were only 55%
# correct. After the refit the top quintile predicts 0.48 vs observed 0.45 and 'change'
# calls are 64% correct. Gene-grouped 5-fold CV on the confident set: AUC 0.865 /
# accuracy 0.871 / R2 0.578 (the all-pairs fit's headline 0.894 was partly the easy
# low-quality rows). The superseded coefficients are kept in
# alphafold_benchmark/foldchange_model_allpairs.json.
#
# `pae_global` (whole canonical structure mean predicted-aligned-error) still dominates
# (std coef +1.32, down from +1.96 - some of that weight was the quality confound);
# `identity` (trimmed-region identity %) is #2; `max_cov_loss` and `protL` add small
# real signal. Burial / region_am / loeuf are ~0 marginal here once PAE is in (they
# drive impact_probability, the functional axis - the fold-vs-function duality, and the
# two scores are empirically independent at rho -0.05).
# Regenerate with alphafold_benchmark/fit_foldchange_hq.py.
_FOLD_CHANGE_MODEL = {
    'features': ['pae_global', 'identity', 'max_cov_loss', 'protL'],
    'median':   [13.8059, 78.1609, 10.0, 388.0],
    'mean':     [13.9957, 69.3368, 31.0288, 381.6161],
    'std':      [6.218, 26.4663, 37.817, 128.6137],
    'coef':     [1.3167, -0.5837, 0.5167, -0.5125],
    'intercept': -2.5145,
}


def fold_change_probability(pae_global=None, identity=None, max_cov_loss=None, protL=None):
    """Calibrated probability (0-1) that the canonical and alternative isoforms have
    DIFFERENT folds - operationally, AlphaFold TM-score < 0.5 - from the rigidity and
    length of the canonical structure (mean PAE, protein length) together with the
    extent of the change (sequence identity, maximum Pfam coverage lost). The
    structural companion to impact_probability (which is functional/pathogenicity).

    What TM<0.5 actually means. It is not a clean 'the domain refolds' label: it fires
    on three different things, and the model cannot tell them apart.
      1. genuine refolding of the retained sequence;
      2. rigid-body RE-ORIENTATION of domains that individually keep their folds -
         the paper's own explanation for most of its high-identity/low-TM outliers
         ("variation in the linker regions ... different orientation of domains");
      3. AlphaFold simply failing on the isoform (long IDRs, isolated helices). The
         refit excludes the worst of these, but the label still carries some.
    TM is also averaged over BOTH length normalisations, so a large deletion depresses
    it mechanically: length ratio alone predicts TM<0.5 at AUC 0.676. Read the score as
    "structurally divergent", not "refolded".

    Calibrated for AlphaFold-CONFIDENT canonical structures (see `_FOLD_CHANGE_MODEL`):
    training kept only pairs whose isoform model was `high`/`confident` quality. On a
    floppy, poorly-predicted protein the score is extrapolating - canonical mean pLDDT
    (afdb_plddt) separates that regime at AUC 0.873 if a guard is ever wanted.

    Defined for substitutions and deletions; pure insertions are left blank (they
    have no canonical changed region, and were dropped from training too).

    Features (see `_FOLD_CHANGE_MODEL`) - two describe the canonical alone, two
    compare it to the alternative, so the score cannot be computed from the
    canonical by itself:
      pae_global   - mean predicted-aligned-error over the canonical AlphaFold
                     structure (higher = floppier/multi-domain -> more fold change).
                     The dominant signal; provisioned via the afdb_pae table.
      protL        - canonical protein length (aa).
      identity     - sequence identity (%) between canonical and alternative protein,
                     from the untouched flanks: how much of the protein survives.
      max_cov_loss - max Pfam coverage (%) lost in the changed region.
    Note these carry the EXTENT of the change, not its LOCATION: position along the
    protein (dist_term / rel_center / frac_before / is_terminal) was evaluated and
    dropped - pae_global absorbs it (+0.017 AUC before PAE, +0.002 after). Placement
    enters only indirectly, when the change overlaps a Pfam domain and so shows up in
    max_cov_loss. A missing feature is imputed to the training median.

    Returns a float in [0, 1], or None if every feature is missing.
    """
    m = _FOLD_CHANGE_MODEL
    vals = [pae_global, identity, max_cov_loss, protL]
    if all(v is None for v in vals):
        return None
    z = m['intercept']
    for i, v in enumerate(vals):
        x = m['median'][i] if v is None else v          # impute missing to the median
        z += m['coef'][i] * ((x - m['mean'][i]) / m['std'][i])
    return round(1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, z)))), 4)


_PLDDT_FLOOR = 70.0   # AlphaFold's own "confident" cutoff (below 70 = low confidence)

# The feature each score is mostly made of. When it is absent the model would fall back
# to the training median - a HUMAN BENCHMARK CONSTANT - and still emit a number, which
# is indistinguishable from a measurement. That is not sporadic missingness: pae_global
# is absent for ~81% of human genes and 100% of non-human ones (afdb_pae was provisioned
# for the 5.5k benchmark accessions only), and region_am for every non-human gene
# (AlphaMissense is human-only). With region_am imputed, impact_prob collapses to a
# constant 0.22-0.25 whatever else changes; with pae_global imputed, fold_change_prob is
# systematically pulled towards 'preserved'. So the score is BLANKED instead, and the
# note says which feature was unavailable.
_DOMINANT = {'fold_change_prob': 'pae_global', 'impact_prob': 'region_am'}
_FEATURE_SOURCE = {           # table the value would have come from
    'pae_global': 'afdb_pae',
    'region_am': 'am_pathogenicity',
    'loeuf': 'gene_constraint',
    'buried_frac': 'afdb_rsa',
    'identity': 'ensembl_sequence',
    'protL': 'ensembl_sequence',
    'max_cov_loss': 'pfam',
}


NO_ACCESSION_NOTE = ('not scored: canonical transcript has no UniProt accession in '
                     'DoChaP, so no structural or constraint data can be looked up')


def missing_features(**values):
    """Names of the features whose value is absent (None or ''), in the order given."""
    return tuple(k for k, v in values.items() if v is None or v == '')


def _note_parts(missing, score):
    """Shared note text: whether the score was withheld, and what was imputed."""
    dominant = _DOMINANT[score]
    parts = []
    withheld = dominant in missing
    if withheld:
        parts.append(f'not scored: {dominant} unavailable '
                     f'- no {_FEATURE_SOURCE[dominant]} entry for this protein')
    others = [m for m in missing if m != dominant]
    if others:
        # only claim imputation when a score was actually produced from it
        parts.append(('also missing: ' if withheld else 'imputed to training median: ')
                     + ', '.join(others))
    return parts


def impact_note(missing=()):
    """Reason string for impact_prob: blank when every feature was measured, else why
    the score was withheld and/or which features fell back to the training median.

    region_am is the dominant feature (+0.80 of a model whose next term is -0.42); with
    it imputed the score is a constant 0.22-0.25 regardless of the other inputs, so it
    is withheld rather than emitted. loeuf and buried_frac missing only weaken it, so
    those are imputed and declared.
    """
    return '; '.join(_note_parts(tuple(missing), 'impact_prob'))


def fold_change_note(canonical_plddt=None, missing=(), plddt_floor=_PLDDT_FLOOR):
    """Reason string when fold_change_prob is outside the regime it was calibrated
    for, or '' when it is inside.

    The model is fit on pairs whose ISOFORM structure AlphaFold resolved well - a
    quality DOMAS can never check, because it does not fold isoforms (AlphaFold DB
    holds one model per canonical accession; there is no P12345-2 entry). Canonical
    mean pLDDT is the available proxy: it predicts the isoform's quality class at AUC
    0.873, because a protein AlphaFold models poorly is modelled poorly in any isoform.

    Below the floor the score is extrapolating, and its errors are not symmetric: on
    those proteins a 'preserved' call is wrong 37.6% of the time versus 6.1% inside the
    regime, while 'change' calls stay 0.946 precise. So fold_change_call downgrades
    'preserved' (only) to 'uncertain', and this note explains why, rather than emitting
    a confident-looking "nothing happened" about a structure nobody can trust.

    The note is emitted on EVERY row below the floor, including ones still called
    'change' or 'uncertain' - it flags an unreliable structure, which is worth knowing
    regardless of which way the call went.

    Caveat: the proxy is imperfect (AUC 0.873), so it is a reduction in exposure, not a
    fix. It removes 46% of the unreliable-regime 'preserved' calls; the survivors are no
    more accurate than before (39.7% wrong). Rows carrying this note deserve a real
    isoform structure (ESMFold) before any structural claim is made about them.

    Returns '' when canonical_plddt is None (unknown != unreliable - the guard only
    fires on positive evidence that the canonical model is poor).
    """
    missing = tuple(missing)
    parts = _note_parts(missing, 'fold_change_prob')
    if _DOMINANT['fold_change_prob'] not in missing:
        # the reliability guard only matters when a call is actually being made
        try:
            v = float(canonical_plddt) if canonical_plddt not in (None, '') else None
        except (TypeError, ValueError):
            v = None
        if v is not None and v < plddt_floor:
            parts.insert(0, f'canonical structure poorly resolved (mean pLDDT {v:.0f})')
    return '; '.join(parts)


def fold_change_call(fold_change_prob, low=0.2, high=0.5, canonical_plddt=None,
                     plddt_floor=_PLDDT_FLOOR):
    """Triage a fold_change_probability into an actionable call. The mid-probability
    band is genuinely ambiguous - the 'does the remainder refold' class that only
    actual folding resolves - so route it out rather than force a call:

        prob >= high  -> 'change'     (structurally divergent)
        prob <= low   -> 'preserved'  (same fold)
        otherwise     -> 'uncertain'  (fold/inspect to resolve)

    The band is ASYMMETRIC because the calibrated base rate is 0.150, not 0.5: on the
    confident set a 0.5 reading is already well above average risk. Re-derived with the
    refit model (E44); out-of-fold on 6,973 confident pairs it gives called-set accuracy
    0.911 at 83% coverage, 'preserved' NPV 0.938, and 'change' precision 0.642 at recall
    0.316. The old 0.4/0.6 band on the old model made 894 'change' calls on the same
    rows at 0.548 precision - i.e. a coin flip labelled 'confident'; 380 of those now
    correctly land in 'uncertain'. Widen towards 0.15/0.45 for a stricter called set
    (same accuracy, 78% coverage).

    canonical_plddt: mean pLDDT of the canonical AlphaFold model. Below `plddt_floor`
    a 'preserved' call is downgraded to 'uncertain' - see fold_change_note for why.
    Emit that note alongside the call so the reason travels with the verdict. Pass None
    (the default) to disable the guard.

    The guard is deliberately ONE-SIDED. 'preserved' asserts that nothing happened
    structurally, which is only as good as the structure it rests on: below the floor
    it is wrong 39.7% of the time. 'change' is driven by sequence-level facts (a large
    deletion, lost Pfam coverage) that hold whether or not the model resolved well, and
    it stays 0.946 precise there - so guarding it would discard good calls for nothing.
    Measured on the benchmark, one-sided guarding cuts 'preserved' error 0.100 -> 0.084
    while keeping all 1,664 'change' calls at 0.851; guarding both sides gives the same
    'preserved' gain but leaves only 492 'change' calls at 0.744.

    Returns None if prob is None/blank.
    """
    if fold_change_prob is None or fold_change_prob == '':
        return None
    p = float(fold_change_prob)
    call = 'change' if p >= high else 'preserved' if p <= low else 'uncertain'
    if call == 'preserved' and fold_change_note(canonical_plddt, plddt_floor=plddt_floor):
        return 'uncertain'
    return call


def calc_spade_score(domains_by_isoform):
    """SPADE-style per-isoform Pfam domain-integrity score.

    SPADE (APPRIS) scores each isoform by how intact its Pfam domains are,
    using the HMMER per-domain **bitscore** as the integrity measure and taking
    the most-intact isoform as the per-domain reference. This reproduces that
    mechanic from our `pfam_match` data (which stores pfam_acc + bitscore per
    protein); it is SPADE-*style*, not byte-identical to the APPRIS code (whose
    exact aggregation formula isn't published).

    Algorithm:
      1. Per isoform, sum the bitscores of each Pfam family's hits
         (so an isoform missing one of three zinc fingers scores lower on that
         family, and a weakened/partial domain scores lower too).
      2. For each family, the reference = the highest summed bitscore across the
         gene's isoforms (the most-intact copy).
      3. An isoform's penalty for a family = reference - its own summed bitscore
         (0 if it holds the best copy; = reference if the family is absent = lost).
      4. spade_score = -sum(penalties): 0 means every domain is held at its best
         (most intact, domain-wise); more negative = more damaged/lost domains.

    Input : {isoform_id: [(pfam_acc, bitscore), ...]}
    Output: {isoform_id: {'spade_score': float, 'lost': [acc,...],
                          'damaged': [acc,...], 'intact': bool}}
    """
    # family -> {isoform: summed bitscore}
    fam_bits = {}
    for iso, doms in domains_by_isoform.items():
        per_fam = {}
        for pfam_acc, bitscore in doms:
            if bitscore is None:
                continue
            per_fam[pfam_acc] = per_fam.get(pfam_acc, 0.0) + float(bitscore)
        for fam, b in per_fam.items():
            fam_bits.setdefault(fam, {})[iso] = b

    best = {fam: max(iso_b.values()) for fam, iso_b in fam_bits.items()}

    result = {}
    for iso in domains_by_isoform:
        penalty = 0.0
        lost, damaged = [], []
        for fam, ref in best.items():
            b = fam_bits[fam].get(iso, 0.0)
            deficit = ref - b
            if deficit > 1e-9:
                penalty += deficit
                (lost if b == 0 else damaged).append(fam)
        result[iso] = {
            'spade_score': round(-penalty, 1),
            'lost': lost,
            'damaged': damaged,
            'intact': penalty <= 1e-9,
        }
    return result

# ---------------------------------------------------------------------------
# Event feature <-> exon matching.
#
# Lives here so junction_analisys.py and generate_gene_pdf.py can share ONE
# implementation: junction_analisys imports generate_gene_pdf, so the PDF layer
# cannot import back from it. The two used to carry separate copies of this
# predicate, with a comment asking that they be kept in step by hand.
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
# means every row is a plain junction, so every existing input file and every
# reader that has not been taught to emit spans keeps working unchanged.
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
# These live here rather than in alternative_splicing.py so that
# junction_analisys.py can import them directly: alternative_splicing.py
# imports junction_analisys, so reaching them from there required a
# function-local import to dodge the circular dependency.
# ---------------------------------------------------------------------------

def get_exons_for_transcripts(con, transcript_ids):
    # Define the maximum chunk size (SQLite limit is 999, so 450 is safe since 450 * 2 = 900)
    chunk_size = 450
    df_list = []
    
    t_ids = list(transcript_ids)  # Ensure it's a list if it's a different iterable
    # Loop through the IDs in chunks
    for i in range(0, len(t_ids), chunk_size):
        chunk_ids = t_ids[i:i + chunk_size]
        
        # Create the correct number of placeholders for this specific chunk
        placeholders = ','.join(['?'] * len(chunk_ids))
        
        # Construct the query dynamically for the chunk
        query = f'''
            SELECT * FROM Transcript_exon 
            WHERE transcript_ensembl_id IN ({placeholders}) 
            OR transcript_refseq_id IN ({placeholders})
        '''
        
        # Duplicate the chunk IDs to match the two IN clauses
        params = chunk_ids * 2
        
        # Execute and append the DataFrame chunk
        df_chunk = pd.read_sql_query(query, con, params=params)
        df_list.append(df_chunk)

    if not df_list:
        # No transcript ids at all, so the chunk loop never ran and there is nothing
        # to concatenate - pd.concat([]) raises "No objects to concatenate", which
        # surfaced as an opaque pandas error rather than an empty result. It happens
        # whenever no event resolves a gene: every LeafCutter cluster unannotated, or
        # a gene set absent from the database. Returning the table's own empty shape
        # lets the run finish and report each event's reason (no_gene_specified /
        # gene_not_in_db) instead of dying.
        return pd.read_sql_query('SELECT * FROM Transcript_exon LIMIT 0', con)

    # Combine all chunks into one final DataFrame
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
    # 450, not 500: gene_id_clause() binds each id once per gene id column, so a
    # batch costs 2 parameters per gene. The same arithmetic as
    # get_exons_for_transcripts() - 450 * 2 = 900, under SQLite's 999 limit.
    batch_size = 450
    for i in range(0, len(gene_ids), batch_size):
        batch = gene_ids[i:i + batch_size]
        clause, params = gene_id_clause(batch)
        query = f'SELECT * FROM Transcripts WHERE {clause}'
        dfs.append(pd.read_sql_query(query, con, params=params))
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def _read_transcripts_and_proteins(con, transcript_ids):
    """Transcripts/Proteins rows for transcript_ids, filtered to a non-empty
    protein_ensembl_id. Shared by get_transcript_domains_db() and
    get_representative_domains_db() so get_domains_db() reads these tables once."""
    df_transcript = pd.read_sql_query('select * from Transcripts', con)
    df_transcript = df_transcript[df_transcript.transcript_ensembl_id.isin(transcript_ids)]
    df_protein = pd.read_sql_query('select * from Proteins', con)
    df_protein = df_protein[df_protein.transcript_ensembl_id.isin(transcript_ids)]
    # Drop rows with empty or NaN protein_ensembl_id
    df_protein = df_protein.dropna(subset=['protein_ensembl_id'])
    df_protein = df_protein[df_protein.protein_ensembl_id.str.strip() != '']
    return df_transcript, df_protein


def get_transcript_domains_db(con, transcript_ids, df_transcript=None, df_protein=None):
    print('Starting getting domains from dochap')
    if df_transcript is None or df_protein is None:
        df_transcript, df_protein = _read_transcripts_and_proteins(con, transcript_ids)
    proteins_ids = np.unique(df_protein.protein_ensembl_id.values).tolist()
    df_domain_event = pd.read_sql_query('select * from DomainEvent', con)

    df_domain_event = df_domain_event[df_domain_event.protein_ensembl_id.isin(proteins_ids)]
    # Drop rows with NaN or empty protein_ensembl_id
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
    print('Done getting domains from dochap')
    print(f'df columns: {merged_df.columns}')
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
    # domain_id + InterPro entry `type` are carried through so
    # junction_analisys.filter_representative_domains() can reduce the domain
    # set by curated type (Domain/Repeat vs Family/Homologous_superfamily)
    # instead of the geometric collapse heuristic. `type` is NULL for DBs built
    # before RepresentativeDomains gained the column (filter then no-ops).
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
    print('Starting getting domains from RepresentativeDomains')
    if df_transcript is None or df_protein is None:
        df_transcript, df_protein = _read_transcripts_and_proteins(con, transcript_ids)

    if 'protein_interpro_id' not in df_protein.columns:
        print('Proteins table has no protein_interpro_id column '
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
        print('RepresentativeDomains table not found in this DB.')
        return pd.DataFrame(columns=REPRESENTATIVE_DOMAINS_COLUMNS)

    df_rep = df_rep[df_rep.protein_interpro_id.isin(interpro_ids)]
    df_rep = df_rep.dropna(subset=['start', 'end'])
    if df_rep.empty:
        return pd.DataFrame(columns=REPRESENTATIVE_DOMAINS_COLUMNS)

    # `type` is absent in DBs built before RepresentativeDomains gained the
    # column; add it as NULL so downstream columns/filtering stay uniform.
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

    print('Done getting domains from RepresentativeDomains')
    return merged_df[REPRESENTATIVE_DOMAINS_COLUMNS]


def get_domains_db(con, transcript_ids, use_representative_domains=False):
    """
    Domain source used by JunctionsAnalysis.analyze_junctions().

    use_representative_domains=False (default): unchanged - domains come from
    DomainEvent/DomainType exactly as get_transcript_domains_db() always has, so the
    existing algorithm runs without any change.

    use_representative_domains=True: domains come from RepresentativeDomains where a
    protein has an entry there; a protein with no RepresentativeDomains entry falls back
    to its DomainEvent/DomainType domains, so no protein silently loses domain coverage.
    """
    df_transcript, df_protein = _read_transcripts_and_proteins(con, transcript_ids)
    df_event_domains = get_transcript_domains_db(con, transcript_ids, df_transcript=df_transcript, df_protein=df_protein)
    if not use_representative_domains:
        return df_event_domains

    df_rep_domains = get_representative_domains_db(con, transcript_ids, df_transcript=df_transcript, df_protein=df_protein)
    proteins_with_rep_domains = set(df_rep_domains['protein_ensembl_id_version'].unique())
    df_fallback = df_event_domains[
        ~df_event_domains['protein_ensembl_id_version'].isin(proteins_with_rep_domains)
    ]
    return pd.concat([df_rep_domains, df_fallback], ignore_index=True)
