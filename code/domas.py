# To read excel file: pip install pandas openpyxl
import pandas as pd
import numpy  as np
import argparse
import sqlite3
from generate_gene_pdf import GeneVisualization



domain_name_columns = ['cdd', 'pfam', 'smart', 'tigr', 'interpro','CDD_id']

def get_gene_exons(conn, df_domains, gene_id):
    gene_id = gene_id.split('.')[0]  # remove version
    # Get all transcripts for this gene
    df_cur_transcripts = df_domains[df_domains.gene_ensembl_id == gene_id]
    
    if len(df_cur_transcripts) == 0:
        return None

    transcript_ids = np.unique(df_cur_transcripts.transcript_ensembl_id_version.values).tolist()
    placeholders = ','.join(['?' for _ in transcript_ids])
    exon_query = f'select * from Transcript_exon where transcript_ensembl_id in ({placeholders})'
    df_exons = pd.read_sql_query(exon_query, conn, params=transcript_ids)
    df_exons['aa_start'] = (df_exons['abs_start_CDS'] / 3).astype(int)
    df_exons['aa_end'] = (df_exons['abs_end_CDS'] / 3).astype(int)
    return df_exons

def find_aa_ranges_skipping_junction2(conn, df_exons, gene_id, start_pos, end_pos):
    # Keep only exons that overlap with the junction region [start_pos+1, end_pos-1]
    matching_transcript_ids = []
    df_start_exons = df_exons[abs(df_exons.genomic_start_tx - end_pos) <= 1]
    df_end_exons = df_exons[abs(df_exons.genomic_end_tx - start_pos) <= 1]
    initial_matching_transcript_ids = set(df_start_exons.transcript_ensembl_id.values) & set(df_end_exons.transcript_ensembl_id.values)
    for transcript_id in initial_matching_transcript_ids:
        start_rank = int(df_start_exons[df_start_exons.transcript_ensembl_id == transcript_id].rank.values[0])
        end_rank = int(df_end_exons[df_end_exons.transcript_ensembl_id == transcript_id].rank.values[0])
        if end_rank - start_rank == 1:
            matching_transcript_ids.append(transcript_id) 

    df_matching = df_exons[df_exons.transcript_ensembl_id.isin(matching_transcript_ids)]
    df_matching = df_matching[df_matching.abs_start_CDS > 0]
    
    df_matching_aa = df_matching.groupby('transcript_ensembl_id').agg({
        'aa_start': 'min',
        'aa_end': 'max'
    }).reset_index()

    df_matching_aa.drop_duplicates(inplace=True)
    return df_matching_aa


def find_aa_ranges_skipping_junction(conn, df_exons, gene_id, start_pos, end_pos):
    # Keep only exons that overlap with the junction region [start_pos+1, end_pos-1]
    skipping_mask = (abs(df_exons.genomic_start_tx - end_pos) <= 1) & (abs(df_exons.genomic_end_tx - start_pos) <= 1)
    
    df_skipping_exons = df_exons[skipping_mask]
    skipping_transcripts = np.unique(df_skipping_exons.transcript_ensembl_id.values)
    non_skipping_transcripts = list(set(df_exons.transcript_ensembl_id.values) - set(skipping_transcripts))
    
    df_exons = df_exons[df_exons.abs_start_CDS > 0]
    df_non_skipping_exons = df_exons[df_exons.transcript_ensembl_id.isin(non_skipping_transcripts)]
    df_non_skipping_exons = df_non_skipping_exons[(df_non_skipping_exons.genomic_start_tx <= end_pos) & (df_non_skipping_exons.genomic_end_tx >= start_pos)]
    
    # Group by transcript_id and get min(aa_start) and max(aa_end) for each transcript
    #df_skipping_transcript_aa_ranges = df_skipping_exons.groupby('transcript_ensembl_id').agg({
    #    'aa_start': 'min',
    #    'aa_end': 'max'
    #}).reset_index()
    df_not_skipping_transcript_aa_ranges = df_non_skipping_exons.groupby('transcript_ensembl_id').agg({
        'aa_start': 'min',
        'aa_end': 'max'
    }).reset_index()

    #df_skipping_transcript_aa_ranges.drop_duplicates(inplace=True)
    df_not_skipping_transcript_aa_ranges.drop_duplicates(inplace=True)

    return df_not_skipping_transcript_aa_ranges
   
 
def get_relevant_domains(df_domains, df_transcript_aa_ranges):
    relvant_columns = ['AA_start', 'AA_end'] + domain_name_columns + ['short_description']
    df_relvant_domains = pd.DataFrame(columns=relvant_columns)
    # Get domains for transcripts based on their AA ranges
    for _, row in df_transcript_aa_ranges.iterrows():
        transcript_id = row['transcript_ensembl_id']
        aa_start = row['aa_start']
        aa_end = row['aa_end']
    
        # Check if 'id' matches AND if coordinates overlap
        mask = (df_domains['transcript_ensembl_id_version'] == transcript_id) & \
            (df_domains['AA_end'] >= aa_start) & \
            (df_domains['AA_start'] <= aa_end)

        # Get the domains
        df_transcript_domains = df_domains.loc[mask, ['AA_start', 'AA_end'] + domain_name_columns + ['short_description']]
        df_relvant_domains = pd.concat([df_relvant_domains, df_transcript_domains], ignore_index=True)
    
    
    df_relvant_domains.drop_duplicates(inplace=True)
    df_relvant_domains = filter_domains_by_name_and_length(df_relvant_domains)
    return df_relvant_domains

def filter_domains_by_name_and_length(df_domains):
    # for rows with same value in domain_name_columns, keep only the one with the longest length (AA_end - AA_start)
    # do not compare name which are nan r empty string
    filtered_df_domains = df_domains.copy()
    filtered_df_domains['domain_length'] = filtered_df_domains['AA_end'] - filtered_df_domains['AA_start']
    filtered_df_domains['domain_name'] = filtered_df_domains[domain_name_columns].apply(
        lambda x: next((name for name in x if pd.notna(name) and name.strip() != '' and name.strip() != 'nan'), None), axis=1)
    filtered_df_domains = filtered_df_domains[filtered_df_domains['domain_name'].notna()]
    filtered_df_domains = filtered_df_domains.sort_values(by=['domain_name', 'domain_length'], ascending=[True, False])
    filtered_df_domains = filtered_df_domains.drop_duplicates(subset=['domain_name'], keep='first')
    filtered_df_domains = filtered_df_domains.drop(columns=['domain_length'])
    return filtered_df_domains


def get_gene_transcript_ids(conn, gene_ids):
    """
    Load transcript domains for a list of genes from the database.
    
    Args:
        conn: Database connection
        gene_ids: List of gene IDs
    
    Returns:
        DataFrame containing domains for all transcripts of the given genes
    """
    # Get transcripts for all genes in batches to avoid SQL query limits
    df_gene_transcripts = pd.DataFrame(columns=['gene_ensembl_id_version', 'transcript_ensembl_id_version'])
    gene_ids_no_version = [gene_id.split('.')[0] for gene_id in gene_ids]
    
    # Process in batches of 500 to avoid SQL IN clause limits
    batch_size = 500
    for i in range(0, len(gene_ids_no_version), batch_size):
        batch = gene_ids_no_version[i:i+batch_size]
        placeholders = ','.join(['?' for _ in batch])
        trans_query = f'select gene_ensembl_id, transcript_ensembl_id from Transcripts where gene_ensembl_id in ({placeholders})'
        cur_df = pd.read_sql_query(trans_query, conn, params=batch)
        cur_df.columns = ['gene_ensembl_id_version', 'transcript_ensembl_id_version']
        df_gene_transcripts = pd.concat([df_gene_transcripts, cur_df], ignore_index=True)
    
    print(f"Loaded {len(df_gene_transcripts)} transcripts for {len(gene_ids)} genes")
    
    transcript_ids = df_gene_transcripts.transcript_ensembl_id_version.values.tolist()
    return np.unique(transcript_ids).tolist()
    
def get_genes_df_transcripts(conn, gene_ids):
    """
    Load all transcripts (with their canonical flag) for a list of genes.

    Args:
        conn: Database connection
        gene_ids: List of gene Ensembl IDs

    Returns:
        DataFrame with one row per transcript, including transcript_ensembl_id,
        transcript_refseq_id, gene_ensembl_id and canonical columns.
    """
    dfs = []
    batch_size = 500
    for i in range(0, len(gene_ids), batch_size):
        batch = gene_ids[i:i + batch_size]
        placeholders = ','.join(['?' for _ in batch])
        query = f'SELECT * FROM Transcripts WHERE gene_ensembl_id IN ({placeholders})'
        dfs.append(pd.read_sql_query(query, conn, params=batch))
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def compare_domains(junctions, cluster_domains):
    if len(cluster_domains) == 1:
        print(f"No domains to compare. {junctions[0]} is the only junction in the cluster, skipping comparison.")
        return ['Nothing to compare'], ['no description']
    results = []
    descriptions = []

    for i in range(len(cluster_domains)):
        cur_results = []
        cur_descriptions = []
        cur_domains_df =  filter_domains_by_name_and_length(cluster_domains[i])
        cur_domains_df['domain_length'] = cur_domains_df['AA_end'] - cur_domains_df['AA_start'] + 1
        other_domains = [cluster_domains[j] for j in range(len(cluster_domains)) if j != i]
        other_domains_df = pd.concat(other_domains, ignore_index=True)
        other_domains_df = filter_domains_by_name_and_length(other_domains_df)
        other_domains_df['domain_length'] = other_domains_df['AA_end'] - other_domains_df['AA_start'] + 1
        # Compare domains between current junction and other junctions in the cluster
        cur_domains_names = set(cur_domains_df['domain_name'].values)
        other_domains_names = set(other_domains_df['domain_name'].values)   
        only_in_cur = cur_domains_names - other_domains_names
        only_in_other = other_domains_names - cur_domains_names
        in_both = cur_domains_names & other_domains_names
        for name in only_in_cur:
            cur_length = cur_domains_df[cur_domains_df['domain_name'] == name].domain_length.values[0]
            print(f"Domain {name}, length {cur_length}, is only in junction {junctions[i]} and not in other junctions in the cluster") 
            cur_results.append(f'only has {name}, {cur_length}')
            cur_descriptions.append(cur_domains_df[cur_domains_df['domain_name'] == name].short_description.values[0])
        for name in only_in_other:
            other_length = other_domains_df[other_domains_df['domain_name'] == name].domain_length.values[0]
            print(f"Domain {name}, length {other_length}, is only in other junctions in the cluster and not in junction {junctions[i]}")
            cur_results.append(f'Missing {name}, {other_length}')  
            cur_descriptions.append(other_domains_df[other_domains_df['domain_name'] == name].short_description.values[0])
        for name in in_both:
            cur_length = cur_domains_df[cur_domains_df['domain_name'] == name].domain_length.values[0]
            other_length = other_domains_df[other_domains_df['domain_name'] == name].domain_length.values[0]
            if cur_length == other_length:
                print(f"Domain {name} has the same length in junction {junctions[i]} and other junctions")
            elif cur_length > other_length:
                print(f"Domain {name} is longer in junction {junctions[i]} than in other junctions. {cur_length} vs {other_length}")
                cur_results.append(f'Longer {name}, {cur_length} vs {other_length}')
                cur_descriptions.append(cur_domains_df[cur_domains_df['domain_name'] == name].short_description.values[0])
            else:       
                print(f"Domain {name} is shorter in junction {junctions[i]} than in other junctions. {cur_length} vs {other_length}")
                cur_results.append(f'Shorter {name}, {cur_length} vs {other_length}')
                cur_descriptions.append(cur_domains_df[cur_domains_df['domain_name'] == name].short_description.values[0])
        results.append(';'.join(cur_results))
        cur_descriptions = [desc.replace(';', ' | ') for desc in cur_descriptions] # replace ';' in descriptions to avoid confusion with the separator
        descriptions.append(';'.join(cur_descriptions))
    return results, descriptions


def analyze_junctions_by_cluster(input_path, dochap_db_path, output_path='cluster_analysis.csv', args=None):
    """
    Main function that:
    1. Reads input file and groups by cluster
    2. For each cluster, finds gene and junction coordinates
    3. Finds transcripts that skip the junction
    4. Finds AA coordinates affected by junction skip
    5. Compares domains across transcripts in the cluster
    """
    # Connect to database
    conn = sqlite3.connect(dochap_db_path)
    
    # Read input file
    input_df = read_input_file(input_path)
    
    if input_df is None:
        print("Error reading input file")
        return None
    
    # Get transcript domains (bulk load for efficiency)
    print("Loading transcript domains from database...")
    gene_ids = np.unique(input_df.ensembl_h.values).tolist()
    transcript_ids = get_gene_transcript_ids(conn, gene_ids)
    df_domains = get_transcript_domains_db(conn, transcript_ids)
    df_results = pd.DataFrame(columns=['cluster', 'gene_name', 'junction_name', 'gene_id', 'start', 'end', 'comparison_results', 'domain_descriptions'])
    
    # Group by cluster
    grouped = input_df.groupby('cluster')
    count = 0
    total = len(grouped)
    for cluster_id, cluster_df in grouped:
        print(f"\n{'='*80}")
        print(f"Processing Cluster: {cluster_id}, {count+1}/{total}")
        print(f"{'='*80}")
        count += 1
        
        cluster_domains = []
        junctions = []
        for idx, row in cluster_df.iterrows():
            gene_id = row['ensembl_h']
            start_pos = row['start_position']
            end_pos = row['end_position']
            junction_name = row['h_junction']
           
            
            print(f"\n  Junction: {junction_name}")
            print(f"  Gene: {gene_id}")
            df_gene_domains = df_domains[df_domains.gene_ensembl_id == gene_id]
            # Step 3: Find transcripts that skip this junction
            df_exons = get_gene_exons(conn, df_domains, gene_id)
            if df_exons is None:
                print(f"Error:  No transcripts found for gene {gene_id}, skipping junction {junction_name}")
                continue
            junctions.append((junction_name, gene_id, start_pos, end_pos))
            df_non_skipping_aa_ranges = find_aa_ranges_skipping_junction(conn, df_exons, gene_id, start_pos, end_pos)
            print(f"  Found {len(df_non_skipping_aa_ranges)} transcripts that don't skip exons at this junction")
            df_non_skipping_domains = get_relevant_domains(df_gene_domains, df_non_skipping_aa_ranges)
            cluster_domains.append(df_non_skipping_domains)
        if len(cluster_domains) == 0: # for testing, if no junctions in the cluster,
            continue
        compare_results, domain_descriptions = compare_domains(junctions, cluster_domains)
        df_junctions= pd.DataFrame(junctions, columns=['junction_name', 'gene_id', 'start', 'end'])
        df_junctions['comparison_results'] = compare_results
        df_junctions['domain_descriptions'] = domain_descriptions
        df_junctions['cluster'] = cluster_id
        gene_id = df_junctions.gene_id.values[0].split(':')[0]
        df_gene = pd.read_sql_query(f'select * from Genes where gene_ensembl_id = "{gene_id}"', conn)
        gene_name = df_gene.gene_symbol.values[0] if len(df_gene) > 0 else gene_id
        df_junctions['gene_name'] = gene_name
        df_results = pd.concat([df_results, df_junctions], ignore_index=True)
        if (not args.gene_ids is None) and gene_id in [g.strip() for g in args.gene_ids.split(',')]:
            print_results(conn, df_junctions, gene_name)
    df_results.to_csv(output_path, index=False)
    conn.close()
   
def print_transcripts_by_cds_length(conn, gene_symbols, output_file=None, species=None):
    """
    For each gene symbol in gene_symbols, list its transcripts in descending
    order of CDS length, printing each transcript's id, cds_start, cds_end
    and cds_length.

    Args:
        conn: DoChaP database connection.
        gene_symbols: list of gene symbols to look up (case-insensitive).
        output_file: if given, write the listing to this path instead of
            printing it to stdout.
        species: if given (e.g. 'H_sapiens'), only consider genes from this
            species - otherwise every species with a matching symbol is
            included.
    """
    if not gene_symbols:
        return

    placeholders = ','.join(['?'] * len(gene_symbols))
    query = f"SELECT gene_ensembl_id, gene_symbol FROM Genes WHERE gene_symbol COLLATE NOCASE IN ({placeholders})"
    params = list(gene_symbols)
    if species:
        query += " AND specie = ?"
        params.append(species)
    df_genes = pd.read_sql_query(query, conn, params=params)

    lines = []
    found_symbols = set()
    for _, gene in df_genes.iterrows():
        gene_ensembl_id = gene['gene_ensembl_id']
        gene_symbol = gene['gene_symbol']
        found_symbols.add(gene_symbol.upper())

        df_transcripts = pd.read_sql_query(
            "SELECT transcript_ensembl_id, transcript_refseq_id, cds_start, cds_end "
            "FROM Transcripts WHERE gene_ensembl_id = ?",
            conn, params=[gene_ensembl_id],
        )
        df_transcripts = df_transcripts.dropna(subset=['cds_start', 'cds_end'])
        df_transcripts['cds_length'] = df_transcripts['cds_end'] - df_transcripts['cds_start']
        df_transcripts = df_transcripts.sort_values('cds_length', ascending=False)

        lines.append(f"Gene: {gene_symbol} ({gene_ensembl_id})")
        for _, t in df_transcripts.iterrows():
            transcript_id = t['transcript_ensembl_id'] if pd.notna(t['transcript_ensembl_id']) else t['transcript_refseq_id']
            lines.append(
                f"\t{transcript_id}\tcds_start={int(t['cds_start'])}\t"
                f"cds_end={int(t['cds_end'])}\tcds_length={int(t['cds_length'])}"
            )

    missing_symbols = [s for s in gene_symbols if s.upper() not in found_symbols]
    for symbol in missing_symbols:
        lines.append(f"Gene: {symbol} - NOT FOUND in database")

    output_text = "\n".join(lines)
    if output_file:
        with open(output_file, 'w') as f:
            f.write(output_text + "\n")
    else:
        print(output_text)


def print_results(conn, df_junctions, gene_name):
    viz = GeneVisualization(conn, gene_name)
    viz.create_pdf(
        gene_name + '_junction_comparison.pdf',
        protein_only=True,
        domains_only=True,
        df_junction=df_junctions
    )


    

    

def get_transcript_domains_db(con, transcript_ids):
    print('Starting getting domains from dochap')
    df_transcript = pd.read_sql_query('select * from Transcripts', con)
    df_transcript = df_transcript[df_transcript.transcript_ensembl_id.isin(transcript_ids)]
    df_protein = pd.read_sql_query('select * from Proteins', con)
    df_protein = df_protein[df_protein.transcript_ensembl_id.isin(transcript_ids)]
    # Drop rows with empty or NaN protein_ensembl_id
    df_protein = df_protein.dropna(subset=['protein_ensembl_id'])
    df_protein = df_protein[df_protein.protein_ensembl_id.str.strip() != '']
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
    merged_df = merged_df.drop(columns=['type_id', 'ext_id', 'name', 'other_name', 'description_x'])
    merged_df = merged_df.drop(columns=['transcript_refseq_id_x', 'tx_start', 'tx_end', 'cds_start', 'cds_end', 'exon_count'])
    merged_df = merged_df.drop(columns=['transcript_refseq_id_y','protein_refseq_id'])
    print('Done getting domains from dochap')
    print(f'df columns: {merged_df.columns}')
    return merged_df


REPRESENTATIVE_DOMAINS_COLUMNS = [
    'protein_ensembl_id_version', 'transcript_ensembl_id_version', 'protein_interpro_id',
    'gene_ensembl_id', 'canonical', 'AA_start', 'AA_end', 'short_description',
    'CDD_id', 'cdd', 'pfam', 'smart', 'tigr', 'interpro',
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


def get_representative_domains_db(con, transcript_ids):
    """
    Domains sourced from the RepresentativeDomains table (populated by
    DoChaP-db/InterProRepresentativeDomains.py) instead of DomainEvent/DomainType.

    Returns a DataFrame with the same columns as get_transcript_domains_db() so it's a
    drop-in replacement for build_domain_lookup(), but only contains rows for proteins
    that actually have a RepresentativeDomains entry - see get_domains_db() for combining
    it with the DomainEvent/DomainType fallback for proteins that don't.
    """
    print('Starting getting domains from RepresentativeDomains')
    df_transcript = pd.read_sql_query('select * from Transcripts', con)
    df_transcript = df_transcript[df_transcript.transcript_ensembl_id.isin(transcript_ids)]
    df_protein = pd.read_sql_query('select * from Proteins', con)
    df_protein = df_protein[df_protein.transcript_ensembl_id.isin(transcript_ids)]
    # Drop rows with empty or NaN protein_ensembl_id
    df_protein = df_protein.dropna(subset=['protein_ensembl_id'])
    df_protein = df_protein[df_protein.protein_ensembl_id.str.strip() != '']

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

    df_rep = df_rep[df_rep.protein_id.isin(interpro_ids)]
    df_rep = df_rep.dropna(subset=['start', 'end'])
    if df_rep.empty:
        return pd.DataFrame(columns=REPRESENTATIVE_DOMAINS_COLUMNS)

    merged_df = pd.merge(df_protein, df_transcript, on=['protein_ensembl_id', 'transcript_ensembl_id'])
    merged_df = pd.merge(merged_df, df_rep, left_on='protein_interpro_id', right_on='protein_id')

    domain_column = merged_df['domain_id'].map(_route_domain_id_to_column)
    for col in ('CDD_id', 'cdd', 'pfam', 'smart', 'tigr', 'interpro'):
        merged_df[col] = merged_df['domain_id'].where(domain_column == col)

    merged_df = merged_df.rename(columns={
        'protein_ensembl_id': 'protein_ensembl_id_version',
        'transcript_ensembl_id': 'transcript_ensembl_id_version',
        'start': 'AA_start',
        'end': 'AA_end',
        'domain_name': 'short_description',
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
    df_event_domains = get_transcript_domains_db(con, transcript_ids)
    if not use_representative_domains:
        return df_event_domains

    df_rep_domains = get_representative_domains_db(con, transcript_ids)
    proteins_with_rep_domains = set(df_rep_domains['protein_ensembl_id_version'].unique())
    df_fallback = df_event_domains[
        ~df_event_domains['protein_ensembl_id_version'].isin(proteins_with_rep_domains)
    ]
    return pd.concat([df_rep_domains, df_fallback], ignore_index=True)


def parse_args():
    parser = argparse.ArgumentParser(description="A script collecting  junction's domains.")
    parser.add_argument("-input", required=False, default="clusters_sum_table_H_vs_M_HN6.xlsx", type=str, help="Path to excel file containig junctions")
    parser.add_argument("-dochap", required=True, type=str, help="Path to dochap db")
    parser.add_argument("-output_csv", type=str, default="junctions.csv", help="Path to the output csv")
    parser.add_argument("-biomart_domains", type=bool, default=False, help="Path to the output csv")
    parser.add_argument("-gene_ids", type=str, default=None, help="Comma-separated list of gene IDs; if provided, print_results is called only for genes in this list")
    args = parser.parse_args()
    return args

# read only relevant columns from input file and
# parse h_junction into chromosome, start_position, end_position
human_relevant_columns_names = [
    'h_junction',
    'symbol_h', 
    'ensembl_h',
    'rank_h',
    'cluster',
]
mouse_relevant_columns_names = [
    'm_junction',
    'genes', 
    'rank_m',
    'cluster'
]
def hadas_read_input_file(con, input_path):
    try:
        df = pd.read_excel(input_path)
        print("Completed reading input file.")
        # create human df with relevant columns and parsed junction coordinates
        df_h = df[human_relevant_columns_names].copy()
        df_h[['chromosome', 'start_position', 'end_position']] = df_h['h_junction'].str.split(':', expand=True)
        df_h['start_position'] = df_h['start_position'].astype(int)
        df_h['end_position'] = df_h['end_position'].astype(int)
        df_h = df_h[(df_h.start_position >= 0) & (df_h.end_position >= 0)]
        df_h['specie'] = 'human'
        df_h.rename(columns={'symbol_h': 'gene_symbol', 'ensembl_h': 'gene_ensembl_id', 'rank_h': 'rank',
                             'h_junction': 'junction_name', 'cluster': 'cluster_name'}, inplace=True)
        df_h = df_h[['gene_symbol', 'gene_ensembl_id', 'junction_name', 'chromosome',
                      'start_position', 'end_position', 'specie', 'cluster_name', 'rank']]

        # create mouse df with relevant columns and parsed junction coordinates
        df_m = df[mouse_relevant_columns_names].copy()
        df_m[['chromosome', 'start_position', 'end_position']] = df_m['m_junction'].str.split(':', expand=True)
        df_m['start_position'] = df_m['start_position'].astype(int)
        df_m['end_position'] = df_m['end_position'].astype(int)
        df_m = df_m[(df_m.start_position >= 0) & (df_m.end_position >= 0)]
        df_m['specie'] = 'mouse'
        # get gene ensembl id from genes column, which is in the format "gene_symbol (gene_ensembl_id)"
        query = "SELECT gene_ensembl_id, gene_symbol FROM Genes WHERE specie = 'M_musculus' AND UPPER(gene_symbol) IN ({gene_symbols})"
        gene_symbols = df_m['genes'].unique().tolist()
        df_mouse_genes = pd.read_sql_query(query.format(gene_symbols=','.join(['?']*len(gene_symbols))), con, params=gene_symbols)
        df_mouse_genes['gene_symbol'] = df_mouse_genes['gene_symbol'].str.upper() 
        df_m = pd.merge(df_m, df_mouse_genes, left_on='genes', right_on='gene_symbol', how='left') 
        df_m.drop(columns=['genes'], inplace=True)     
        df_m.rename(columns={ 'm_junction': 'junction_name', 'cluster': 'cluster_name', 'rank_m': 'rank'}, inplace=True)
        df_m = df_m[['gene_symbol', 'gene_ensembl_id', 'junction_name', 'chromosome', 
                     'start_position', 'end_position', 'specie', 'cluster_name', 'rank']]
        merged_df = pd.concat([df_h, df_m], ignore_index=True)

        print("Completed parsing junctions.")
        return merged_df
    except Exception as e:
        print(f"Error reading input file: {e}")
        raise(e)

def main():
    pd.set_option('display.max_colwidth', None)
    args = parse_args()
    analyze_junctions_by_cluster(args.input, args.dochap, args.output_csv, args=args)


if __name__ == "__main__":
    main()





