import pandas as pd

def ioe2junctions(file_path):
    """
    Parses an IOE file and returns a DataFrame with the relevant data.

    Args:
        file_path (str): The path to the IOE file.
    """
    junctions = []
    df_ioe = pd.read_csv(file_path, sep='\t')
    count = 0
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
        if event_type == 'SE':
            # SE is giving the junction between the skipped exon and the flanking exons, so we need to add the junction between the flanking exons
            if len(current_junctions) != 2:
                continue
            start_short, end_short = current_junctions[0]
            start_long, end_long = current_junctions[0][0], current_junctions[1][1]
            junctions.append([chromosome, gene_ensembl_id, event_type, start_short, end_short, cluster_name])
            junctions.append([chromosome, gene_ensembl_id, event_type, start_long, end_long, cluster_name])
        else:
            for start, end in current_junctions:
                 junctions.append([chromosome, gene_ensembl_id, event_type, start, end, cluster_name])


    df_junctions = pd.DataFrame(junctions, columns=['chromosome', 'gene_ensembl_id', 'event_type', 'start_position', 'end_position', 'cluster_name'])
    return df_junctions

def get_gene_symbols(con, gene_ensembl_ids):
    """
    Retrieves gene symbols for a list of Ensembl gene IDs.

    Args:
        gene_ensembl_ids (list): A list of Ensembl gene IDs.
    Returns:
        dict: A dictionary mapping Ensembl gene IDs to gene symbols.
    """
    placeholders = ', '.join(['?'] * len(gene_ensembl_ids))
    query = f"SELECT gene_ensembl_id, gene_symbol FROM genes WHERE gene_ensembl_id IN ({placeholders})"
    df = pd.read_sql_query(query, con, params=gene_ensembl_ids)
    gene_symbol_dict = dict(zip(df['gene_ensembl_id'], df['gene_symbol']))
    return gene_symbol_dict


def get_genes_number_of_transcripts(con, gene_ensembl_ids):
    """
    Retrieves the number of transcripts for a list of Ensembl gene IDs.

    Args:
        gene_ensembl_ids (list): A list of Ensembl gene IDs.
    Returns:
        dict: A dictionary mapping Ensembl gene IDs to the number of transcripts.
    """
    placeholders = ', '.join(['?'] * len(gene_ensembl_ids))
    query = f"SELECT gene_ensembl_id, COUNT(COALESCE(transcript_ensembl_id, transcript_refseq_id)) AS num_transcripts FROM transcripts WHERE gene_ensembl_id IN ({placeholders}) GROUP BY gene_ensembl_id"
    df = pd.read_sql_query(query, con, params=gene_ensembl_ids)
    num_transcripts_dict = dict(zip(df['gene_ensembl_id'], df['num_transcripts']))
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
    placeholders = ','.join('?' * len(gene_ensembl_ids))
    query = f'''
        SELECT gene_ensembl_id, exon_count
            FROM Transcripts
            WHERE gene_ensembl_id IN ({placeholders})
              AND canonical != 0
        '''
    cur.execute(query, list(gene_ensembl_ids))
    for gene_ensembl_id, exon_count in cur.fetchall():
        result[gene_ensembl_id] = exon_count

    return result