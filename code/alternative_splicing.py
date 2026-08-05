import logging
import os
import re
import pandas as pd

from junction_analisys import JunctionsAnalysis
import utils

logger = logging.getLogger(__name__)

# read only relevant columns from input file and
# parse h_junction into chromosome, start_position, end_position
human_relevant_columns_names = [
    'h_junction',
    'symbol_h',
    'ensembl_h',
    'cluster',
]
mouse_relevant_columns_names = [
    'm_junction',
    'genes',
    'cluster'
]


def hadas_read_input_file(con, input_path):
    try:
        df = pd.read_excel(input_path)
        logger.info("Completed reading input file.")
        # create human df with relevant columns and parsed junction coordinates
        df_h = df[human_relevant_columns_names].copy()
        df_h[['chromosome', 'start_position', 'end_position']] = df_h['h_junction'].str.split(':', expand=True)
        df_h['start_position'] = df_h['start_position'].astype(int)
        df_h['end_position'] = df_h['end_position'].astype(int)
        df_h = df_h[(df_h.start_position >= 0) & (df_h.end_position >= 0)]
        df_h['specie'] = 'human'
        df_h.rename(columns={'symbol_h': 'gene_symbol', 'ensembl_h': 'gene_ensembl_id',
                             'h_junction': 'junction_name', 'cluster': 'cluster_name'}, inplace=True)
        df_h = df_h[['gene_symbol', 'gene_ensembl_id', 'junction_name', 'chromosome',
                      'start_position', 'end_position', 'specie', 'cluster_name']]

        # create mouse df with relevant columns and parsed junction coordinates
        df_m = df[mouse_relevant_columns_names].copy()
        df_m[['chromosome', 'start_position', 'end_position']] = df_m['m_junction'].str.split(':', expand=True)
        df_m['start_position'] = df_m['start_position'].astype(int)
        df_m['end_position'] = df_m['end_position'].astype(int)
        df_m = df_m[(df_m.start_position >= 0) & (df_m.end_position >= 0)]
        df_m['specie'] = 'mouse'
        # get gene ensembl id from genes column, which is in the format "gene_symbol (gene_ensembl_id)"
        # gene_GeneID_id is the fallback for the 13% of genes with no
        # gene_ensembl_id, which otherwise resolve to NaN and are dropped.
        query = "SELECT gene_ensembl_id, gene_GeneID_id, gene_symbol FROM Genes WHERE specie = 'M_musculus' AND UPPER(gene_symbol) IN ({gene_symbols})"
        gene_symbols = df_m['genes'].unique().tolist()
        df_mouse_genes = pd.read_sql_query(query.format(gene_symbols=','.join(['?']*len(gene_symbols))), con, params=gene_symbols)
        df_mouse_genes['gene_ensembl_id'] = utils.combined_gene_ids(df_mouse_genes)
        df_mouse_genes = df_mouse_genes.drop(columns=['gene_GeneID_id'])
        df_mouse_genes['gene_symbol'] = df_mouse_genes['gene_symbol'].str.upper()
        df_m = pd.merge(df_m, df_mouse_genes, left_on='genes', right_on='gene_symbol', how='left')
        df_m.drop(columns=['genes'], inplace=True)
        df_m.rename(columns={ 'm_junction': 'junction_name', 'cluster': 'cluster_name'}, inplace=True)
        df_m = df_m[['gene_symbol', 'gene_ensembl_id', 'junction_name', 'chromosome',
                     'start_position', 'end_position', 'specie', 'cluster_name']]
        merged_df = pd.concat([df_h, df_m], ignore_index=True)

        logger.info("Completed parsing junctions.")
        return merged_df
    except Exception as e:
        logger.error("Error reading input file: %s", e)
        raise(e)


# --- LeafCutter differential-splicing output readers -------------------------
# The two standard leafcutter_ds outputs, reshaped into the per-junction schema
# the other readers produce:
#   - leafcutter_ds_effect_sizes.txt: one row per intron, keyed by an `intron`
#     id of the form chr:start:end:clu_<n>_<strand>.
#   - leafcutter_ds_cluster_significance.txt: one row per cluster, keyed by
#     chr:clu_<n>_<strand>, carrying the `genes` column.

# Maps the specie label carried on the junctions (matching the hadas convention:
# 'human'/'mouse') to the DoChaP Genes.specie value used for symbol lookups.
# Moved to utils.SPECIE_DB_NAME so junction_analisys.py can cross-check the
# stated species against the database without a circular import.
_SPECIE_DB_NAME = utils.SPECIE_DB_NAME


def _split_leafcutter_intron(intron):
    """chr:start:end:clu_<n>_<strand> -> (chromosome, start, end, cluster_name),
    where cluster_name is chr:clu_<n> (strand dropped, to match the hadas
    `cluster` format)."""
    parts = str(intron).split(':')
    chromosome = parts[0] if len(parts) > 0 else ''
    start = parts[1] if len(parts) > 1 else ''
    end = parts[2] if len(parts) > 2 else ''
    clu = parts[3] if len(parts) > 3 else ''
    # clu is clu_<n>_<strand>; _leafcutter_cluster_key() drops the trailing
    # _<strand>, so both files' ids are normalised by one function.
    return chromosome, start, end, _leafcutter_cluster_key(f'{chromosome}:{clu}')


def _leafcutter_gene_symbols(genes):
    """The gene symbols a leafcutter_ds significance row names, as a list.

    An unannotated cluster gives []. The missing-value check matters: str() on a
    missing value yields the string 'nan', which would reach the output CSV as a
    gene symbol.
    """
    if genes is None or (not isinstance(genes, str) and pd.isna(genes)):
        return []
    symbols = []
    for symbol in str(genes).split(','):
        symbol = symbol.strip()
        if symbol and symbol.lower() not in ('nan', 'na', '.', 'none') and symbol not in symbols:
            symbols.append(symbol)
    return symbols


def _leafcutter_cluster_key(cluster):
    """chr:clu_<n>_<strand> (significance file) -> chr:clu_<n>, the same key
    _split_leafcutter_intron() builds from the effect-sizes intron id.

    Older LeafCutter versions write the cluster without the strand suffix
    (clu_<n>), so only a trailing _+ / _- is dropped. Splitting on the last
    underscore unconditionally turned every clu_<n> into 'clu', which merged
    every cluster on a chromosome into a single one - junctions from unrelated
    genes analysed together, silently.
    """
    parts = str(cluster).split(':')
    if len(parts) < 2:
        return str(cluster)
    chromosome, clu = parts[0], parts[1]
    clu_no_strand = clu
    if clu.startswith('clu_'):
        head, _sep, tail = clu.rpartition('_')
        if tail in ('+', '-'):
            clu_no_strand = head
    return f'{chromosome}:{clu_no_strand}'


def _leafcutter_attach_genes(con, df, cluster_to_symbols, specie):
    """Attach the gene annotation to a per-junction frame keyed by cluster_name,
    given the gene symbols each cluster names.

    LeafCutter clusters introns annotation-free and only then names the genes each
    overlaps, so a cluster can name several - 6.6% of clusters in a real run do, up
    to nine. Each is analysed once per named gene, so a change in any of them is
    reported and a name absent from DoChaP does not cost the others. Only a
    multi-gene cluster carries the symbol in its identifier.

    Shared by the two LeafCutter-shaped readers: the pair of leafcutter_ds files
    and the internal2 Excel export of the same results.
    """
    df = df.copy()
    df['gene_symbols'] = df['cluster_name'].map(cluster_to_symbols)
    # A cluster naming no gene keeps one row with a missing symbol - it is a real
    # LeafCutter outcome (a novel cluster overlapping nothing annotated), not bad
    # input, and is labelled no_gene_specified downstream rather than dropped.
    df['gene_symbols'] = df['gene_symbols'].apply(lambda s: s if isinstance(s, list) and s else [None])
    df = df.explode('gene_symbols', ignore_index=True)
    df = df.rename(columns={'gene_symbols': 'gene_symbol'})

    multi_gene_clusters = {key for key, symbols in cluster_to_symbols.items() if len(symbols) > 1}
    is_multi = df['cluster_name'].isin(multi_gene_clusters)
    df.loc[is_multi, 'cluster_name'] = (df.loc[is_multi, 'cluster_name']
                                        + ':' + df.loc[is_multi, 'gene_symbol'].astype(str))

    # resolve gene symbol -> gene_ensembl_id via the DoChaP Genes table.
    # The same symbol exists across species in the DB (e.g. USP16 -> human,
    # mouse, rat, ...), so the lookup must be restricted to this input's species.
    db_specie = _SPECIE_DB_NAME.get(specie, specie)
    symbols = [s for s in df['gene_symbol'].dropna().unique().tolist()
               if s and str(s).lower() not in ('nan', 'na', '.')]
    symbol_to_ensembl = {}
    if symbols:
        placeholders = ','.join(['?'] * len(symbols))
        # gene_GeneID_id is selected as a fallback: 13% of DoChaP genes carry no
        # gene_ensembl_id, and resolving a symbol to that column alone yields NaN,
        # dropping the gene before analysis (see utils.combined_gene_ids).
        query = (f"SELECT gene_ensembl_id, gene_GeneID_id, gene_symbol FROM Genes "
                 f"WHERE specie = ? AND UPPER(gene_symbol) IN ({placeholders})")
        df_genes = pd.read_sql_query(query, con, params=[db_specie] + [s.upper() for s in symbols])
        symbol_to_ensembl = {sym.upper(): gid
                             for gid, sym in zip(utils.combined_gene_ids(df_genes),
                                                 df_genes['gene_symbol'])}
    df['gene_ensembl_id'] = df['gene_symbol'].apply(
        lambda s: symbol_to_ensembl.get(str(s).upper()) if pd.notna(s) else None)

    df['specie'] = specie
    return df


def leafcutter_read_input_files(con, significance_file, effect_sizes_file, specie='human'):
    """Read a pair of LeafCutter differential-splicing outputs
    (leafcutter_ds_cluster_significance.txt + leafcutter_ds_effect_sizes.txt) and
    return a per-junction DataFrame with the same columns hadas_read_input_file()
    produces, ready for analyze_junctions().

    The effect-sizes file supplies the junction coordinates and cluster; the
    significance file supplies the gene symbol per cluster, which is resolved to a
    gene_ensembl_id via the DoChaP `Genes` table.
    """
    # leafcutter_ds writes tab-separated files; the significance `status` column
    # contains spaces (e.g. "<2 introns used in >=min_samples_per_intron
    # samples"), so splitting on whitespace would corrupt it.
    df_effect = pd.read_csv(effect_sizes_file, sep='\t')
    df_sig = pd.read_csv(significance_file, sep='\t')

    if 'intron' not in df_effect.columns:
        raise ValueError(f"'intron' column not found in effect-sizes file {effect_sizes_file}")
    if 'cluster' not in df_sig.columns or 'genes' not in df_sig.columns:
        raise ValueError(f"'cluster'/'genes' columns not found in significance file {significance_file}")

    # per-junction coordinates + cluster from the effect-sizes file
    df = df_effect.copy()
    # Unpack the parsed intron id into named columns in one pass, rather than
    # re-applying a positional lambda over the tuple Series four times.
    parsed = pd.DataFrame(
        df['intron'].map(_split_leafcutter_intron).tolist(),
        columns=['chromosome', 'start_position', 'end_position', 'cluster_name'],
        index=df.index,
    )
    df['chromosome'] = parsed['chromosome']
    df['start_position'] = parsed['start_position'].astype(int)
    df['end_position'] = parsed['end_position'].astype(int)
    df['cluster_name'] = parsed['cluster_name']
    df['junction_name'] = (df['chromosome'] + ':' + df['start_position'].astype(str)
                           + ':' + df['end_position'].astype(str))
    df = df[(df.start_position >= 0) & (df.end_position >= 0)].copy()

    # cluster -> the gene symbols the significance file names for it.
    df_sig['cluster_key'] = df_sig['cluster'].apply(_leafcutter_cluster_key)
    cluster_to_symbols = dict(zip(df_sig['cluster_key'],
                                  df_sig['genes'].apply(_leafcutter_gene_symbols)))
    df = _leafcutter_attach_genes(con, df, cluster_to_symbols, specie)

    columns = ['gene_symbol', 'gene_ensembl_id', 'junction_name', 'chromosome',
               'start_position', 'end_position', 'specie', 'cluster_name']
    return df[columns].reset_index(drop=True)


def analyze_leafcutter_input(con, significance_file, effect_sizes_file, output_csv,
                             specie='human', num_workers=5, use_representative_domains=False, max_clusters=0,
                             filter_non_comparable=False, write_all_comparable=False):
    """Read a pair of leafcutter_ds output files and run the domain analysis, the
    same way analyze_ioe_file()/analyze_hadas_input() do for their formats."""
    df_junctions = leafcutter_read_input_files(con, significance_file, effect_sizes_file, specie=specie)
    analyze_junctions(con, df_junctions=df_junctions, output_path=output_csv, create_pdf=False,
                       num_workers=num_workers, use_representative_domains=use_representative_domains,
                       max_clusters=max_clusters, filter_non_comparable=filter_non_comparable,
                       write_all_comparable=write_all_comparable)


# --- internal2: LeafCutter results as a supplementary-table Excel -------------
# A collaborator's export of a leafcutter_ds run, with the cluster-significance
# and effect-size outputs already joined into a single sheet - one row per
# junction, under a title block:
#
#     | Supplementary Table 11                                                  |
#     | DTUs between TB and PTB placentas                                       |
#     |                                                                         |
#     | Cluster    | Splicing junction        | TB | PTB | deltaPSI | p.Adjust | Genes    |
#     | clu_1455_- | chr1:153626533-153627139 | .. | ..  | ..       | ..       | S100A13  |
#
# The table does not start at A1 (and in some files not in column A either), and
# the two PSI columns are named after the compared conditions, so the header row
# and its columns are found by label rather than by position. Only Cluster,
# Splicing junction and Genes are read: these rows are already the significant
# subset, and DOMAS applies no significance filtering of its own, exactly as for
# -format leafcutter.
_INTERNAL2_REQUIRED_LABELS = frozenset({'cluster', 'splicing junction', 'genes'})
# chr1:153626533-153627139. The leafcutter_ds files spell the same thing
# chr:start:end, so both separators are accepted.
_INTERNAL2_JUNCTION_RE = r'^(?P<chromosome>[^:\s]+):(?P<start_position>\d+)[-:](?P<end_position>\d+)$'


def _internal2_locate_table(df_raw, input_path):
    """Find the results table inside a sheet that starts with a title block.

    Returns (index of the header row, {lower-cased header label: column}).
    """
    for row_index in range(len(df_raw)):
        labels = {}
        for column, value in df_raw.iloc[row_index].items():
            if pd.notna(value):
                labels.setdefault(str(value).strip().lower(), column)
        if _INTERNAL2_REQUIRED_LABELS.issubset(labels):
            return row_index, labels
    raise ValueError(
        f"No results table found in {input_path}: expected a header row naming "
        f"{', '.join(sorted(_INTERNAL2_REQUIRED_LABELS))}")


def internal2_read_input_file(con, input_path, specie='human'):
    """Read an internal2 Excel file (LeafCutter differential-splicing results as a
    supplementary table) and return a per-junction DataFrame with the same columns
    hadas_read_input_file() produces, ready for analyze_junctions().

    The junction column supplies the coordinates, the cluster column the cluster -
    normalised to the chr:clu_<n> key leafcutter_read_input_files() builds, so the
    two readers produce the same identifiers for the same run - and the Genes
    column the symbols, resolved to a gene_ensembl_id via the DoChaP `Genes` table.
    """
    df_raw = pd.read_excel(input_path, header=None)
    header_row, labels = _internal2_locate_table(df_raw, input_path)
    body = df_raw.iloc[header_row + 1:]
    df = pd.DataFrame({
        'cluster': body[labels['cluster']],
        'junction': body[labels['splicing junction']],
        'genes': body[labels['genes']],
    })
    # Trailing blank rows below the table, and any row that names no junction,
    # carry nothing to analyse; they are not the unreadable-junction case warned
    # about below.
    df = df[df['cluster'].notna() & df['junction'].notna()].copy()

    parsed = df['junction'].astype(str).str.strip().str.extract(_INTERNAL2_JUNCTION_RE)
    unreadable = parsed['start_position'].isna()
    if unreadable.any():
        # These sheets are hand-assembled and do contain the occasional corrupted
        # cell (a stray spreadsheet reference pasted over a value), so one bad row
        # drops rather than failing a run of thousands.
        logger.warning("internal2: skipping %d row(s) with an unreadable junction, e.g. %s",
                       int(unreadable.sum()), df.loc[unreadable, 'junction'].head(3).tolist())
        df = df[~unreadable]
        parsed = parsed[~unreadable]
    if df.empty:
        raise ValueError(f"No readable junctions found in {input_path}")

    df['chromosome'] = parsed['chromosome']
    df['start_position'] = parsed['start_position'].astype(int)
    df['end_position'] = parsed['end_position'].astype(int)
    df['junction_name'] = (df['chromosome'] + ':' + df['start_position'].astype(str)
                           + ':' + df['end_position'].astype(str))
    # clu_<n>_<strand> -> chr:clu_<n>; the cluster column names no chromosome, so
    # it comes from the junction, as it does in the leafcutter_ds intron id.
    df['cluster_name'] = (df['chromosome'] + ':' + df['cluster'].astype(str).str.strip()
                          ).apply(_leafcutter_cluster_key)

    # The Genes column repeats the cluster's annotation on every one of its rows;
    # collecting it per cluster (rather than per row) keeps a cluster analysed once
    # per named gene even where the rows disagree.
    cluster_to_symbols = {}
    for cluster_name, genes in zip(df['cluster_name'], df['genes']):
        symbols = cluster_to_symbols.setdefault(cluster_name, [])
        for symbol in _leafcutter_gene_symbols(genes):
            if symbol not in symbols:
                symbols.append(symbol)

    df = _leafcutter_attach_genes(con, df, cluster_to_symbols, specie)
    logger.info("internal2: read %d junctions in %d clusters from %s",
                len(df), df['cluster_name'].nunique(), input_path)

    columns = ['gene_symbol', 'gene_ensembl_id', 'junction_name', 'chromosome',
               'start_position', 'end_position', 'specie', 'cluster_name']
    return df[columns].reset_index(drop=True)


def analyze_internal2_input(con, input_file, output_csv, specie='human', num_workers=5,
                            use_representative_domains=False, max_clusters=0,
                            filter_non_comparable=False, write_all_comparable=False):
    """Read an internal2 Excel file and run the domain analysis, the same way
    analyze_leafcutter_input() does for the leafcutter_ds files it exports."""
    df_junctions = internal2_read_input_file(con, input_file, specie=specie)
    analyze_junctions(con, df_junctions=df_junctions, specie=specie, output_path=output_csv, create_pdf=False,
                      num_workers=num_workers, use_representative_domains=use_representative_domains,
                      max_clusters=max_clusters, filter_non_comparable=filter_non_comparable,
                      write_all_comparable=write_all_comparable)



def read_junctions_csv(junctions_csv):
    """Read a plain junctions CSV (chromosome, gene_ensembl_id, start_position,
    end_position, cluster_name, ...) into a DataFrame ready for
    JunctionsAnalysis.analyze_junctions()."""
    return pd.read_csv(junctions_csv)


def _limit_clusters(df_junctions, max_clusters):
    """Keep only the first `max_clusters` clusters (sorted, for determinism). Used
    to cap the amount of work - e.g. the web GUI processes the first 100 clusters."""
    if not max_clusters or max_clusters <= 0 or df_junctions is None:
        return df_junctions
    col = 'cluster_name' if 'cluster_name' in df_junctions.columns else 'cluster'
    keep = sorted(df_junctions[col].dropna().unique())[:max_clusters]
    return df_junctions[df_junctions[col].isin(keep)].copy()


def analyze_junctions(con, df_junctions=None, junctions_csv=None, hadas_format=False, specie=None,
                        output_path='as_events_junctions_analysis.csv',
                        n=0, create_pdf=True, print_genes=None, num_workers=5,
                        use_representative_domains=False, max_clusters=0,
                        filter_non_comparable=False, write_all_comparable=False):
    """Read junctions (from a DataFrame, a plain CSV, or a hadas-format Excel file - exactly
    one of df_junctions/junctions_csv must be given) and run JunctionsAnalysis.analyze_junctions().

    max_clusters > 0 caps the analysis to the first N clusters (sorted)."""
    if df_junctions is None and junctions_csv is None:
        raise ValueError("Either df_junctions or junctions_csv must be provided.")
    if df_junctions is not None and junctions_csv is not None:
        raise ValueError("Only one of df_junctions or junctions_csv should be provided.")
    if junctions_csv is not None:
        df_junctions = hadas_read_input_file(con, junctions_csv) if hadas_format else read_junctions_csv(junctions_csv)

    df_junctions = _limit_clusters(df_junctions, max_clusters)

    analyzer = JunctionsAnalysis(con, logger_instance=logger)
    return analyzer.analyze_junctions(df_junctions=df_junctions, output_path=output_path, specie=specie,
                                      filter_transcript_count=n, create_pdf=create_pdf, print_genes=print_genes, num_workers=num_workers,
                                      use_representative_domains=use_representative_domains,
                                      filter_non_comparable=filter_non_comparable,
                                      write_all_comparable=write_all_comparable)

def analyze_ioe_file(con, ioe_file, output_csv, specie=None, num_workers=5, use_representative_domains=False, max_clusters=0,
                     filter_non_comparable=False, write_all_comparable=False):
    df_junctions = utils.ioe2junctions(ioe_file)
    gene_symbols_dict = utils.get_gene_symbols(con, df_junctions.gene_ensembl_id.unique().tolist())
    # add gene symbols to df_junctions
    df_junctions['gene_symbol'] = df_junctions['gene_ensembl_id'].map(gene_symbols_dict)
    analyze_junctions(con, df_junctions=df_junctions, specie=specie, output_path=output_csv, create_pdf=False, num_workers=num_workers,
                        use_representative_domains=use_representative_domains, max_clusters=max_clusters,
                        filter_non_comparable=filter_non_comparable,
                        write_all_comparable=write_all_comparable)


def analyze_rmats_input(con, rmats_dir, output_csv, specie=None, num_workers=5, use_representative_domains=False, max_clusters=0,
                        filter_non_comparable=False, write_all_comparable=False):
    """Read an rMATS-turbo output directory (the SE/A5SS/A3SS/MXE [Event].MATS.JC.txt files)
    and run the domain analysis. rMATS embeds the Ensembl GeneID and gene symbol
    in each event, so - unlike the leafcutter path - no symbol->ensembl lookup is
    needed. All events are analyzed (no significance filtering)."""
    df_junctions = utils.rmats2junctions(rmats_dir)
    analyze_junctions(con, df_junctions=df_junctions, specie=specie, output_path=output_csv, create_pdf=False,
                       num_workers=num_workers, use_representative_domains=use_representative_domains,
                       max_clusters=max_clusters, filter_non_comparable=filter_non_comparable,
                       write_all_comparable=write_all_comparable)


def analyze_voila_input(con, voila_tsv, output_csv, specie=None, num_workers=5, use_representative_domains=False, max_clusters=0,
                        filter_non_comparable=False, write_all_comparable=False):
    """Read a MAJIQ voila TSV (`voila tsv` output) and run the domain analysis.
    voila embeds the Ensembl Gene ID and gene name in each LSV, so no
    symbol->ensembl lookup is needed. All LSVs are analyzed (no filtering)."""
    df_junctions = utils.voila2junctions(voila_tsv)
    analyze_junctions(con, df_junctions=df_junctions, specie=specie, output_path=output_csv, create_pdf=False,
                       num_workers=num_workers, use_representative_domains=use_representative_domains,
                       max_clusters=max_clusters, filter_non_comparable=filter_non_comparable,
                       write_all_comparable=write_all_comparable)


def keep_min_transcript_clusters(df, examples_per_event=2):
    """Filters a DataFrame to keep all rows for the top N unique cluster_names per

    event_type that have the lowest number of transcripts.
    """
    # 1. Get the minimum transcript count per unique cluster
    cluster_mins = (
        df.groupby(["event_type", "cluster_name"])["num_transcripts"]
        .min()
        .reset_index()
    )

    # 2. Use a safer lambda structure that preserves columns during the groupby
    top_clusters = cluster_mins.groupby("event_type", group_keys=False).apply(
        lambda x: x.nsmallest(examples_per_event, "num_transcripts")
    )

    # 3. Safely pull unique names directly from the resulting column
    valid_cluster_names = top_clusters["cluster_name"].unique()

    # 4. Filter the original dataframe
    filtered_df = df[df["cluster_name"].isin(valid_cluster_names)].copy()

    return filtered_df

def analyze_ioe_files(con, input_path, pattern, output_csv, specie=None, examples_per_event=0, num_workers=5,
                       use_representative_domains=False, max_clusters=0, filter_non_comparable=False, write_all_comparable=False):
    dfs = []
    for file in os.listdir(input_path):
        if re.match(pattern, file):
            ioe_file = os.path.join(input_path, file)
            df_junctions = utils.ioe2junctions(ioe_file)
            dfs.append(df_junctions)
    if not dfs:
        logger.warning(f"No files matching pattern {pattern} were found in {input_path}.")
        return
    df_all_junctions = pd.concat(dfs, ignore_index=True)
    gene_symbols_dict = utils.get_gene_symbols(con, df_all_junctions.gene_ensembl_id.unique().tolist())
    df_all_junctions['gene_symbol'] = df_all_junctions['gene_ensembl_id'].map(gene_symbols_dict)
    gene_ensembl_ids = df_all_junctions['gene_ensembl_id'].unique().tolist()
    gene_transcript_counts = utils.get_genes_number_of_transcripts(con, gene_ensembl_ids)
    df_all_junctions['num_transcripts'] = df_all_junctions['gene_ensembl_id'].map(gene_transcript_counts)
    canonical_exon_counts = utils.get_canonical_exon_counts(con, gene_ensembl_ids)
    df_all_junctions['num_canonical_exons'] = df_all_junctions['gene_ensembl_id'].map(canonical_exon_counts)
    #df_all_junctions = df_all_junctions[df_all_junctions['num_transcripts'] >= 2]
    #df_all_junctions = df_all_junctions[df_all_junctions['num_canonical_exons'] <= 6]
    if examples_per_event > 0:
        # per event type, keep only examples_per_event unique cluster_name with minimum number of transcripts for the gene
        df_examples = keep_min_transcript_clusters(df_all_junctions, examples_per_event=examples_per_event)
        df_examples.to_csv('ioe_example_junctions.csv', index=False)
        analyze_junctions(con, df_junctions=df_examples, specie=specie, output_path=output_csv, create_pdf=False, num_workers=num_workers,
                            use_representative_domains=use_representative_domains, max_clusters=max_clusters,
                            filter_non_comparable=filter_non_comparable,
                            write_all_comparable=write_all_comparable)
    else:
        df_all_junctions.to_csv('ioe_all_junctions.csv', index=False)
        analyze_junctions(con, df_junctions=df_all_junctions, specie=specie, output_path=output_csv, create_pdf=False, num_workers=num_workers,
                            use_representative_domains=use_representative_domains, max_clusters=max_clusters,
                            filter_non_comparable=filter_non_comparable,
                            write_all_comparable=write_all_comparable)

def analyze_hadas_input(con, input_file, output_csv, print_genes=None, num_workers=5,
                         use_representative_domains=False, create_pdf=True, max_clusters=0,
                         filter_non_comparable=False, write_all_comparable=False):
    df_junctions = hadas_read_input_file(con, input_file)
    analyze_junctions(con, df_junctions=df_junctions, output_path=output_csv, create_pdf=create_pdf, print_genes=print_genes,
                        num_workers=num_workers, use_representative_domains=use_representative_domains, max_clusters=max_clusters,
                        filter_non_comparable=filter_non_comparable,
                        write_all_comparable=write_all_comparable)
    

