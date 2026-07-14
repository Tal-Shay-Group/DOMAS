import glob
import logging
import os
import re
import sys
import sqlite3
from unittest import case
import pandas as pd
import numpy as np

from generate_gene_pdf import GeneVisualization
from junction_analisys import JunctionsAnalysis
import utils

logger = logging.getLogger(__name__)

def get_first_exon(df_e, strand):
    idx = max(df_e['order_in_transcript'].values) if (strand == '-') else 1
    return df_e[df_e['order_in_transcript'] == idx]
    

def get_last_exon(df_e, strand):
    idx = 1 if (strand == '-') else max(df_e['order_in_transcript'].values)
    return df_e[df_e['order_in_transcript'] == idx]

def check_mutually_exclusive_exons(gene_id, df_e1, df_e2, results):
    max1 = max(df_e1['order_in_transcript'].values)
    for cur_order in range(1,max1-2):
        df_cur_exon = df_e1[df_e1['order_in_transcript'] == cur_order]
        cur_start =  df_cur_exon.genomic_start_tx.values[0]
        cur_end =  df_cur_exon.genomic_end_tx.values[0]
        df_cur_exon2 = df_e2[(df_e2.genomic_start_tx == cur_start) & (df_e2.genomic_end_tx == cur_end)]
        if df_cur_exon2.empty:
            continue  # no matching first exon
        df_cur_exon2 = df_cur_exon2.iloc[0]
              
        last_exon = df_e1[df_e1['order_in_transcript'] == (cur_order + 2)]
        last_exon_start =  last_exon.genomic_start_tx.values[0]
        last_exon_end =  last_exon.genomic_end_tx.values[0]

        if df_e2[(df_e2.genomic_start_tx == last_exon_start) & 
                (df_e2.genomic_end_tx == last_exon_end) & 
                (df_e2.order_in_transcript == (df_cur_exon2.order_in_transcript + 2))].empty:
            continue  # no matching fourth exon

        middle_exon1 = df_e1[df_e1['order_in_transcript'] == (cur_order + 1)]
        middle_exon2 = df_e2[df_e2['order_in_transcript'] == (df_cur_exon2.order_in_transcript + 1)]
        if middle_exon1.empty or middle_exon2.empty:
            continue  # no middle exon
        check_overlap_middle = check_overlap(middle_exon1.genomic_start_tx.values[0], middle_exon1.genomic_end_tx.values[0],
                                            middle_exon2.genomic_start_tx.values[0], middle_exon2.genomic_end_tx.values[0])
        if not check_overlap_middle:
            t_id1 = df_e1.transcript_id.values[0]
            t_id2 = df_e2.transcript_id.values[0]
            results.append({'gene_id': gene_id, 'transcript1': t_id1, 'transcript2': t_id2, 'event': 'mutually_exclusive_exons', 
                            'rank1': cur_order+1, 'rank2': middle_exon2.order_in_transcript.values[0]})     
            logger.info(f'\tsplicing event: mutually_exclusive_exons {t_id1}, rank1={cur_order+1}, {t_id2}, rank2={middle_exon2.order_in_transcript.values[0]}')
            logger.info(f'\t\t({middle_exon1.genomic_start_tx.values[0]},{middle_exon1.genomic_end_tx.values[0]}), ({middle_exon2.genomic_start_tx.values[0]},{middle_exon2.genomic_end_tx.values[0]})')


'''
def check_mutually_exclusive_exons(df_e1, df_e2):
    # Do not overlap with each other, but their predecessor and next exon overlap with the same flanking exons
    min_order1 = df_e1.order_in_transcript.min()
    max_order1 = df_e1.order_in_transcript.max()
    min_order2 = df_e2.order_in_transcript.min()
    max_order2 = df_e2.order_in_transcript.max()

    # Candidate mutually-exclusive exons must be internal (have both predecessor and next exons).
    for exon1 in df_e1.itertuples(index=False):
        order1 = exon1.order_in_transcript
        if (order1 <= min_order1) or (order1 >= max_order1):
            continue

        prev1_df = df_e1[df_e1.order_in_transcript == (order1 - 1)]
        next1_df = df_e1[df_e1.order_in_transcript == (order1 + 1)]
        if prev1_df.empty or next1_df.empty:
            continue
        prev1 = prev1_df.iloc[0]
        next1 = next1_df.iloc[0]

        for exon2 in df_e2.itertuples(index=False):
            order2 = exon2.order_in_transcript
            if (order2 <= min_order2) or (order2 >= max_order2):
                continue

            # Mutually-exclusive exon candidates should not overlap.
            if check_overlap(exon1.genomic_start_tx, exon1.genomic_end_tx, exon2.genomic_start_tx, exon2.genomic_end_tx):
                continue

            prev2_df = df_e2[df_e2.order_in_transcript == (order2 - 1)]
            next2_df = df_e2[df_e2.order_in_transcript == (order2 + 1)]
            if prev2_df.empty or next2_df.empty:
                continue
            prev2 = prev2_df.iloc[0]
            next2 = next2_df.iloc[0]

            prev_overlap = check_overlap(prev1.genomic_start_tx, prev1.genomic_end_tx, prev2.genomic_start_tx, prev2.genomic_end_tx)
            next_overlap = check_overlap(next1.genomic_start_tx, next1.genomic_end_tx, next2.genomic_start_tx, next2.genomic_end_tx)
            if prev_overlap and next_overlap:
                return True

    return False
'''
'''
def check_alt_5_prime(df_e1, df_e2, strand):
    range1 =  range(df_e1['order_in_transcript'].max() -1, 0, -1) if (strand == '-') else range(2, df_e1['order_in_transcript'].max() + 1)
    range2 =  range(df_e2['order_in_transcript'].max() -1, 0, -1) if (strand == '-') else range(2, df_e2['order_in_transcript'].max() + 1)
    for rank1 in range1:
        exon1 = df_e1[df_e1['order_in_transcript'] == rank1].iloc[0]
        for rank2 in range2:
            exon2 = df_e2[df_e2['order_in_transcript'] == rank2].iloc[0]
            if (exon1.genomic_start_tx != exon2.genomic_start_tx) and check_overlap(exon1.genomic_start_tx, exon1.genomic_end_tx, exon2.genomic_start_tx, exon2.genomic_end_tx):
                logger.info(f'splicing event: alt_3_prime {df_e1.transcript_ensembl_id.values[0]}, {df_e2.transcript_ensembl_id.values[0]}')
                logger.info(f'\trank1={rank1}, rank2={rank2}, ({exon1.genomic_start_tx},{exon1.genomic_end_tx}), ({exon2.genomic_start_tx},{exon2.genomic_end_tx})')
'''
def check_alt_3_prime(gene_id, df_e1, df_e2, strand, results):
    check_alt_prime(gene_id, df_e1, df_e2, strand, prime=3, results=results)

def check_alt_5_prime(gene_id, df_e1, df_e2, strand, results):
    check_alt_prime(gene_id, df_e1, df_e2, strand, prime=5, results=results)

def check_alt_prime(gene_id, df_e1, df_e2, strand, prime=3, results=None):
    max1 = df_e1['order_in_transcript'].max()
    max2 = df_e2['order_in_transcript'].max()
    for idx, exon1 in df_e1.iterrows():
        rank1 = exon1.order_in_transcript
        if rank1 == 1 or rank1 == max1:
            continue # first or last alternative
        if (strand == '-' and prime == 3) or (strand == '+' and prime == 5):
            df_exon2 = df_e2[df_e2.genomic_start_tx == exon1.genomic_start_tx]
            if df_exon2.empty:
                continue
            df_exon2 = df_exon2.iloc[0]
            if df_exon2.genomic_end_tx == exon1.genomic_end_tx or df_exon2.order_in_transcript in [1,max2]:
                continue
           
            t_id1 = df_e1.transcript_id.values[0]
            t_id2 = df_e2.transcript_id.values[0]
            results.append({'gene_id': gene_id, 'transcript1': t_id1, 'transcript2': t_id2, 'event': f'alt_prime', 'rank1': rank1, 'rank2': df_exon2.order_in_transcript})
            logger.info(f'\tsplicing event: alt_prime {t_id1}, {t_id2}, rank1={rank1}, rank2={df_exon2.order_in_transcript}')
            logger.info(f'\t\t(exon1: {exon1.genomic_start_tx},{exon1.genomic_end_tx}), exon2: ({df_exon2.genomic_start_tx},{df_exon2.genomic_end_tx})')
        else:
            df_exon2 = df_e2[df_e2.genomic_end_tx == exon1.genomic_end_tx]
            if df_exon2.empty:
                continue
            df_exon2 = df_exon2.iloc[0]
            if df_exon2.genomic_start_tx == exon1.genomic_start_tx or df_exon2.order_in_transcript in [1,max2]:
                continue
            t_id1 = df_e1.transcript_id.values[0]
            t_id2 = df_e2.transcript_id.values[0]
            results.append({'gene_id': gene_id, 'transcript1': t_id1, 'transcript2': t_id2, 'event': f'alt_prime', 'rank1': rank1, 'rank2': df_exon2.order_in_transcript})
            logger.info(f'\tsplicing event: alt_prime {t_id1}, {t_id2}, rank1={rank1}, rank2={df_exon2.order_in_transcript}')
            logger.info(f'\t\t(exon1: {exon1.genomic_start_tx},{exon1.genomic_end_tx}), exon2: ({df_exon2.genomic_start_tx},{df_exon2.genomic_end_tx})')

                

def check_overlap(start1, end1, start2, end2):
    return (start1 <= end2) and (start2 <= end1)

def check_contain(start1, end1, start2, end2):
    return (start1 <= start2) and (end1 >= end2)

def check_skip_exon(gene_id, df_e1, df_e2, results):
    actual_check_skip_exon(gene_id, df_e1, df_e2, results)
    actual_check_skip_exon(gene_id, df_e2, df_e1, results)

def actual_check_skip_exon(gene_id, df_e1, df_e2, results):
    max1 = max(df_e1['order_in_transcript'].values)
    for cur_order in range(1,max1-2):
        df_cur_exon = df_e1[df_e1['order_in_transcript'] == cur_order]
        cur_start =  df_cur_exon.genomic_start_tx.values[0]
        cur_end =  df_cur_exon.genomic_end_tx.values[0]
        df_cur_exon2 = df_e2[(df_e2.genomic_start_tx == cur_start) & (df_e2.genomic_end_tx == cur_end)]
        if df_cur_exon2.empty:
            continue  # no matching first exon
        df_cur_exon2 = df_cur_exon2.iloc[0]
        
        next_exon = df_e1[df_e1['order_in_transcript'] == (cur_order + 2)]
        next_start =  next_exon.genomic_start_tx.values[0]
        next_end =  next_exon.genomic_end_tx.values[0]
        df_next_exon2 = df_e2[(df_e2.genomic_start_tx == next_start) & (df_e2.genomic_end_tx == next_end) & (df_e2.order_in_transcript == (df_cur_exon2.order_in_transcript + 1))]
        if df_next_exon2.empty:
            continue  # no matching third exon
        t_id1 = df_e1.transcript_id.values[0]
        t_id2 = df_e2.transcript_id.values[0]
        rank2_1 = df_cur_exon2.order_in_transcript
        rank2_2 = df_next_exon2.order_in_transcript.values[0]
        results.append({'gene_id': gene_id, 'transcript1': t_id1, 'transcript2': t_id2, 'event': 'skip_exon', 'rank1': cur_order+1, 'rank2': rank2_1, 'rank2_2': rank2_2})
        logger.info(f'\tsplicing event: skip_exon {t_id1}, rank={cur_order+1}, {t_id2}, ranks2={rank2_1}, {rank2_2}')

def check_alt_first_exon(gene_id, df_e1, df_e2, strand, results):
    df_first_exon1 = get_first_exon(df_e1, strand)
    df_first_exon2 = get_first_exon(df_e2, strand)
    first_exon1 = df_first_exon1.iloc[0]
    first_exon2 = df_first_exon2.iloc[0]
    if (first_exon1.genomic_start_tx != first_exon2.genomic_start_tx) or (first_exon1.genomic_end_tx != first_exon2.genomic_end_tx):
        results.append({'gene_id' : gene_id, 'transcript1': df_e1.transcript_id.values[0], 'transcript2': df_e2.transcript_id.values[0], 'event': 'alt_first_exon'})
        logger.info(f'\tsplicing event: alt_first_exon {df_e1.transcript_id.values[0]}, {df_e2.transcript_id.values[0]}')
        logger.info(f'\t\t({first_exon1.genomic_start_tx},{first_exon1.genomic_end_tx}), ({first_exon2.genomic_start_tx},{first_exon2.genomic_end_tx})')

def check_alt_last_exon(gene_id, df_e1, df_e2, strand, results):
    df_last_exon1 = get_last_exon(df_e1, strand)
    df_last_exon2 = get_last_exon(df_e2, strand)
    last_exon1 = df_last_exon1.iloc[0]
    last_exon2 = df_last_exon2.iloc[0]
    if (last_exon1.genomic_start_tx != last_exon2.genomic_start_tx) or (last_exon1.genomic_end_tx != last_exon2.genomic_end_tx):
        results.append({'gene_id': gene_id, 'transcript1': df_e1.transcript_id.values[0], 'transcript2': df_e2.transcript_id.values[0], 'event': 'alt_last_exon'})
        logger.info(f'\tsplicing event: alt_last_exon {df_e1.transcript_id.values[0]}, {df_e2.transcript_id.values[0]}')
        logger.info(f'\t\t({last_exon1.genomic_start_tx},{last_exon1.genomic_end_tx}), ({last_exon2.genomic_start_tx},{last_exon2.genomic_end_tx})')

def check_intron_retention(gene_id, df_e1, df_e2, strand, results):
    actual_check_intron_retention(gene_id, df_e1, df_e2, strand, results) 
    actual_check_intron_retention(gene_id, df_e2, df_e1, strand, results)                

def actual_check_intron_retention(gene_id, df_e1, df_e2, strand, results):
    max1 = max(df_e1['order_in_transcript'].values)
    for cur_order in range(1,max1+1):
        df_cur_exon = df_e1[df_e1['order_in_transcript'] == cur_order]
        cur_start =  df_cur_exon.genomic_start_tx.values[0]
        cur_end =  df_cur_exon.genomic_end_tx.values[0]
        df_exon_start2 = df_e2[(df_e2.genomic_start_tx == cur_start)]
        df_exon_end2 = df_e2[(df_e2.genomic_end_tx == cur_end)]
        if df_exon_start2.empty or df_exon_end2.empty or (df_exon_start2.order_in_transcript.values[0] == df_exon_end2.order_in_transcript.values[0]):
            continue  # no matching overlapping exons
        t_id1 = df_e1.transcript_id.values[0]
        t_id2 = df_e2.transcript_id.values[0]
        rank2_1 = df_exon_start2.order_in_transcript.values[0]
        rank2_2 = df_exon_end2.order_in_transcript.values[0]
        logger.info(f'splicing event: intron_retention {t_id1}, rank1={cur_order}, {t_id2}, ranks2={rank2_1}, {rank2_2}')
        results.append({'gene_id' : gene_id, 'transcript1': t_id1, 'transcript2': t_id2, 'event': 'intron_retention', 'rank1': cur_order, 'rank2': rank2_1, 'rank2_2': rank2_2})
        logger.info(f'\t\t ({cur_start},{cur_end})')

def compare_transcripts(gene_id, df_exons, t_id1, t_id2, strand, results):
    logger.info(f"Comparing {t_id1} to {t_id2}, strand={strand}")
    df_e1 = df_exons[(df_exons.transcript_ensembl_id == t_id1) | (df_exons.transcript_refseq_id == t_id1)].copy()
    df_e1['transcript_id'] = t_id1
    df_e2 = df_exons[(df_exons.transcript_ensembl_id == t_id2) | (df_exons.transcript_refseq_id == t_id2)].copy()
    df_e2['transcript_id'] = t_id2
    #check_skip_exon(gene_id, df_e1, df_e2, results)
    #check_alt_first_exon(gene_id, df_e1, df_e2, strand, results)
    #check_alt_last_exon(gene_id, df_e1, df_e2, strand, results)
    #check_intron_retention(gene_id, df_e1, df_e2, strand, results)
    check_mutually_exclusive_exons(gene_id, df_e1, df_e2, results)
    #check_alt_3_prime(gene_id, df_e1, df_e2, strand, results)
    #check_alt_5_prime(gene_id, df_e1, df_e2, strand, results)




def get_gene_alternative_splicing_events_by_name(con, gene_name, results=None):
    df_gene = pd.read_sql_query(f'select * from genes where gene_symbol == "{gene_name.upper()}"', con)
    strand = df_gene.strand.values[0]
    logger.info(df_gene.columns)
    gene_id = df_gene.iloc[0].gene_ensembl_id
    df_transcripts = pd.read_sql_query(f'select * from transcripts where gene_ensembl_id == "{gene_id}"', con)
    logger.info(df_transcripts.columns)
    get_gene_alternative_splicing_events_by_transcripts(con, df_transcripts, strand, results)

def get_gene_id_from_transcripts(df_transcripts):
    gene_id = df_transcripts.gene_ensembl_id.fillna(df_transcripts.gene_GeneID_id).values[0]
    return gene_id

def get_gene_alternative_splicing_events_by_transcripts(con, df_transcripts, strand, results):
    transcript_ids = df_transcripts.transcript_ensembl_id.fillna(df_transcripts.transcript_refseq_id).values.tolist()
    df_canonical = df_transcripts[df_transcripts.canonical != 0]
    if df_canonical.empty:
        gene_id = get_gene_id_from_transcripts(df_transcripts)
        results.append({'gene_id': gene_id, 'transcript1': None, 'transcript2': None, 'event': 'no_canonical_transcript'})
        logger.warning("No canonical transcript found for this gene. Skipping AS event analysis.")
        return
    canonical_id = df_canonical.transcript_ensembl_id.fillna(df_canonical.transcript_refseq_id).values[0]
    placeholders = ','.join(['?' for _ in transcript_ids])
    exon_query = f'select * from Transcript_exon where transcript_ensembl_id in ({placeholders}) or transcript_refseq_id in ({placeholders})'
    df_exons = pd.read_sql_query(exon_query, con, params=transcript_ids + transcript_ids)
    gene_id = get_gene_id_from_transcripts(df_transcripts)
    for transcript_id in transcript_ids:
        if transcript_id == canonical_id:
            continue
        compare_transcripts(gene_id, df_exons, canonical_id, transcript_id, strand, results)
    #compare_transcripts(df_exons, 'ENST00000440279.7', 'ENST00000488478.5', strand)

def sort_transcripts_by_exon_count(con, output_csv_path, only_canonical):
    if only_canonical:
        query = 'select * from transcripts where canonical != 0 order by exon_count asc'
    else:
        query = 'select * from transcripts order by exon_count asc'
    df_transcripts = pd.read_sql_query(query, con)
    df_transcripts.to_csv(output_csv_path, index=False)

def get_genes_with_multiple_transcripts(con):
    query = 'select gene_ensembl_id, count(*) as transcript_count from transcripts group by gene_ensembl_id having transcript_count > 1'
    df = pd.read_sql_query(query, con)
    return df.gene_ensembl_id.values.tolist()

# --- GTF Export Function ---
def export_genes_to_gtf(con, gene_ids, output_gtf_path):
    """
    Export genes, their transcripts, and exons for a list of gene_ids to GTF format.
    """
    # Get gene info
    placeholders = ','.join(['?'] * len(gene_ids))
    df_genes = pd.read_sql_query(f'SELECT * FROM Genes WHERE gene_ensembl_id IN ({placeholders})', con, params=gene_ids)
    # Get transcripts for these genes
    df_transcripts = pd.read_sql_query(f'SELECT * FROM Transcripts WHERE gene_ensembl_id IN ({placeholders})', con, params=gene_ids)
    transcript_ids = df_transcripts.transcript_ensembl_id.tolist()
    if not transcript_ids:
        logger.warning("No transcripts found for given gene_ids.")
        return
    placeholders_tx = ','.join(['?'] * len(transcript_ids))
    # Get exons for these transcripts
    df_exons = pd.read_sql_query(f'SELECT * FROM Transcript_Exon WHERE transcript_ensembl_id IN ({placeholders_tx})', con, params=transcript_ids)

    # Merge to get gene and transcript info for each exon
    df_exons = df_exons.merge(df_transcripts[['transcript_ensembl_id', 'gene_ensembl_id']], on='transcript_ensembl_id', how='left')
    df_exons = df_exons.merge(df_genes[['gene_ensembl_id', 'chromosome', 'strand', 'gene_symbol']], on='gene_ensembl_id', how='left')

    with open(output_gtf_path, 'w') as gtf:
        # Write gene features
        for _, gene in df_genes.iterrows():
            attr = f'gene_id "{gene.gene_ensembl_id}"; gene_name "{gene.gene_symbol}";'
            gtf.write(f'{gene.chromosome}\tcustom\tgene\t1\t-1\t.\t{gene.strand}\t.\t{attr}\n')
        # Write transcript features
        for _, tx in df_transcripts.iterrows():
            attr = f'gene_id "{tx.gene_ensembl_id}"; transcript_id "{tx.transcript_ensembl_id}";'
            gtf.write(f'chrNA\tcustom\ttranscript\t{tx.tx_start}\t{tx.tx_end}\t.\t.\t.\t{attr}\n')
        # Write exon features
        for _, exon in df_exons.iterrows():
            attr = f'gene_id "{exon.gene_ensembl_id}"; transcript_id "{exon.transcript_ensembl_id}"; exon_number "{exon.order_in_transcript}";'
            gtf.write(f'{exon.chromosome}\tcustom\texon\t{exon.genomic_start_tx}\t{exon.genomic_end_tx}\t.\t{exon.strand}\t.\t{attr}\n')
    logger.info(f"Exported GTF to {output_gtf_path}")

def combine_and_expand_ioe(file_pattern):
    all_files = glob.glob(file_pattern)
    if not all_files:
        return None

    df_list = []
    for file in all_files:
        temp_df = pd.read_csv(file, sep='\t')
        # Keep track of which event type file this came from (SE, A3, AL, etc.)
        temp_df['event_class'] = os.path.basename(file).split('_')[-2] 
        df_list.append(temp_df)

    combined_df = pd.concat(df_list, ignore_index=True)

    # Split the event_id by ':'
    # We use expand=True to create separate columns
    expanded = combined_df['event_id'].str.split(':', expand=True)

    # Define a generic naming scheme for the maximum possible columns
    # Based on your data, there are 7 columns
    col_names = [
        'gene_and_type', # e.g. ENSG...;A3
        'chrom',         # e.g. 1
        'coord_1',       # e.g. 1435821-1436926
        'coord_2', 
        'coord_3', 
        'coord_4', 
        'strand'
    ]

    # Assign names only up to the number of columns actually created
    expanded.columns = col_names[:len(expanded.columns)]

    # Further split the 'gene_and_type' to separate the Gene ID from the Event Type
    gene_type_split = expanded['gene_and_type'].str.split(';', expand=True)
    expanded['gene_id_clean'] = gene_type_split[0]
    expanded['specific_event_type'] = gene_type_split[1]

    # Combine back to the main dataframe
    final_df = pd.concat([combined_df, expanded], axis=1)
    return final_df

def get_all_as_events(con, max_transcripts_per_gene=-1):
    df_genes = pd.read_sql_query(f'select * from genes', con)
    df_transcripts = pd.read_sql_query(f'select * from transcripts', con)
    gene_count = df_transcripts.value_counts('gene_GeneID_id')
    df_transcripts['gene_count'] = df_transcripts['gene_GeneID_id'].map(gene_count)
    if max_transcripts_per_gene > 0:
        df_transcripts = df_transcripts[df_transcripts['gene_count'] <= max_transcripts_per_gene]
    
    gene_ids = df_transcripts.gene_GeneID_id.unique().tolist()
    all_results = []
    for idx, gene_id in enumerate(gene_ids):
        try:
            logger.info(f"Processing gene {idx+1}/{len(gene_ids)}: {gene_id}")
            strand = df_genes[df_genes.gene_GeneID_id == gene_id].strand.values[0]
            df_gene_transcripts = df_transcripts[df_transcripts.gene_GeneID_id == gene_id]
            get_gene_alternative_splicing_events_by_transcripts(con, df_gene_transcripts, strand, all_results)
            #if idx == 100:  # For testing, limit to first 10 genes
            #    break
        except Exception as e:
            logger.error(f"Error: processing gene {gene_id}: {e}")
    df_results = pd.DataFrame(all_results)
    df_results.to_csv('all_as_events.csv', index=False)

def find_clusters_with_n_transcripts(con, input_file, output_file, n=3):
    # get number of transcripts for each gene
    df_t = pd.read_sql_query(f'select * from transcripts', con)
    gene_count = df_t.value_counts('gene_GeneID_id')
    df_t['gene_count'] = df_t['gene_GeneID_id'].map(gene_count)
    ids_n = df_t[df_t['gene_count'] == n].gene_ensembl_id.unique().tolist()
    # get genes from clusters
    df_junctions = hadas_read_input_file(con, input_file)
    df_junctions_n = df_junctions[df_junctions['gene_ensembl_id'].isin(ids_n)]
    analyze_junctions2(con, df_junctions=df_junctions_n, output_path=output_file)

def get_df_junction_columns():
    columns =['chromosome', 'start_position', 'end_position', 'gene_symbol', 'gene_ensembl_id', 'cluster_name']
    return columns

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

    # Combine all chunks into one final DataFrame
    df_exons = pd.concat(df_list, ignore_index=True)
    return df_exons


def get_genes_df_transcripts(con, gene_ids):
    """
    Load all transcripts (with their canonical flag) for a list of genes.

    Args:
        con: Database connection
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
        dfs.append(pd.read_sql_query(query, con, params=batch))
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


def create_events_junctions(con, df_events, output_csv):
    all_junctions = []
    
    df_genes = pd.read_sql_query(f'select * from genes', con)
    gene_ids_x = df_events.gene_id.unique().tolist() # mix of gene_ensembl_id and gene_GeneID_id
    gene_placeholders = ','.join(['?']*len(gene_ids_x))
    df_genes = pd.read_sql_query(f'select * from genes where gene_ensembl_id in ({gene_placeholders}) or gene_GeneID_id in ({gene_placeholders})', con, params=gene_ids_x*2 )  
    gene_ids = df_genes.gene_GeneID_id.unique().tolist()
    df_transcripts = pd.read_sql_query(f'select * from transcripts where gene_GeneID_id in ({gene_placeholders})', con, params=gene_ids)
    transcript_ids = df_transcripts.transcript_ensembl_id.fillna(df_transcripts.transcript_refseq_id).tolist()  
    df_exons = get_exons_for_transcripts(con, transcript_ids)

    # Precompute transcript -> exons and transcript -> exon-by-rank maps for O(1) lookup in the loop.
    transcript_to_rows = {}
    #for pos, exon_row in enumerate(df_exons.itertuples(index=False)):
    #    t_ens = exon_row.transcript_ensembl_id
    #    t_ref = exon_row.transcript_refseq_id
    #    if pd.notna(t_ens):
    #        transcript_to_rows.setdefault(t_ens, []).append(pos)
    #    if pd.notna(t_ref) and t_ref != t_ens:
    #        transcript_to_rows.setdefault(t_ref, []).append(pos)


    df_exons['pos'] = range(len(df_exons))
    df_melted = df_exons.melt(
        id_vars=['pos'], 
        value_vars=['transcript_ensembl_id', 'transcript_refseq_id'], 
        value_name='transcript_id'
    )
    df_cleaned = df_melted.dropna(subset=['transcript_id'])
    df_cleaned = df_cleaned.drop_duplicates(subset=['transcript_id', 'pos'])
    transcript_to_rows = df_cleaned.groupby('transcript_id')['pos'].apply(list).to_dict()

    df_exons_mapped = df_cleaned.merge(
        df_exons[['order_in_transcript']].reset_index(), 
        left_on='pos', 
        right_on='index'
    )

    transcript_to_max_rank = df_exons_mapped.groupby('transcript_id')['order_in_transcript'].max().to_dict()
    transcript_to_exons = {
        transcript_id: df_exons.iloc[row_positions]
        for transcript_id, row_positions in transcript_to_rows.items()
    }
    ExonTuple = df_exons.itertuples(index=False)
    exon_tuples_list = list(ExonTuple) 
    transcript_to_rank_map = (
        df_exons_mapped.groupby('transcript_id')
            .apply(lambda x: {row.order_in_transcript: exon_tuples_list[int(row.pos)] for row in x.itertuples(index=False)})
            .to_dict()
    )
    #transcript_to_rank_map = {
    #    transcript_id: {ex.order_in_transcript: ex for ex in exons_df.itertuples(index=False)}
    #    for transcript_id, exons_df in transcript_to_exons.items()
    #}
    #transcript_to_max_rank = {
    #    transcript_id: exons_df['order_in_transcript'].max()
    #    for transcript_id, exons_df in transcript_to_exons.items()
    #}

    #transcript_to_gene_id = {}
    #for tx in df_transcripts.itertuples(index=False):
    #    if pd.notna(tx.transcript_ensembl_id):
    #        transcript_to_gene_id[tx.transcript_ensembl_id] = tx.gene_GeneID_id
    #    if pd.notna(tx.transcript_refseq_id):
    #       transcript_to_gene_id[tx.transcript_refseq_id] = tx.gene_GeneID_id
    
    df_tx_melted = df_transcripts.melt(
        id_vars=['gene_GeneID_id'],
        value_vars=['transcript_ensembl_id', 'transcript_refseq_id'],
        value_name='transcript_id'
    ).dropna(subset=['transcript_id'])


    transcript_to_gene_id = dict(zip(df_tx_melted['transcript_id'], df_tx_melted['gene_GeneID_id']))
    gene_meta = {
        g.gene_GeneID_id: g
        for g in df_genes.itertuples(index=False)
    }

    def get_exon(rank_map, rank):
        exon = rank_map.get(rank)
        if exon is None and isinstance(rank, float) and rank.is_integer():
            exon = rank_map.get(int(rank))
        return exon

    for idx, row in enumerate(df_events.itertuples(index=False)):
        if idx % 100 == 0:
            logger.info(f"Processing event {idx+1}/{len(df_events)}: gene_id={row.gene_id}, transcript1={row.transcript1}, transcript2={row.transcript2}, event={row.event}")
        junction1, junction2, junction3, junction4 = None, None, None, None
        gene_id = row.gene_id
        t_id1 = row.transcript1
        t_id2 = row.transcript2
        event = row.event
        rank1 = row.rank1 
        rank2 = row.rank2
        rank2_2_raw = getattr(row, 'rank2_2', np.nan)
        rank2_2 = None if pd.isna(rank2_2_raw) else rank2_2_raw

        df_exons1 = transcript_to_exons.get(t_id1)
        df_exons2 = transcript_to_exons.get(t_id2)
        if df_exons1 is None or df_exons2 is None:
            logger.warning(f"Skipping event for {gene_id} because transcript exons are missing.")
            continue
        if len(df_exons1) <= 1 or len(df_exons2) <= 1:
            logger.warning(f"Skipping event for {gene_id} because one of the transcripts has insufficient exons.")
            continue

        exons1_by_rank = transcript_to_rank_map[t_id1]
        exons2_by_rank = transcript_to_rank_map[t_id2]
        
        if event == 'no_canonical_transcript':
            logger.info(f"Gene {gene_id} has no canonical transcript. Skipping junction creation.")
            continue
        elif event == 'alt_first_exon':
            rank1 = 1
            rank2 = 1
        elif event == 'alt_last_exon':
            rank1 = transcript_to_max_rank[t_id1]
            rank2 = transcript_to_max_rank[t_id2]
        logger.debug(f"Creating junctions for event {event} between {t_id1} and {t_id2}, ranks: {rank1}, {rank2}, {rank2_2}")

        exon1 = get_exon(exons1_by_rank, rank1)
        exon2_1 = get_exon(exons2_by_rank, rank2)
        exon2_2 = None if rank2_2 is None else get_exon(exons2_by_rank, rank2_2)
        if exon1 is None or exon2_1 is None:
            logger.warning(f"Skipping event for {gene_id} because an expected exon rank is missing.")
            continue

        cur_gene_id = transcript_to_gene_id.get(t_id1)
        gene_cur_gene = gene_meta.get(cur_gene_id)
        if gene_cur_gene is None:
            logger.warning(f"Skipping event for {gene_id} because transcript-to-gene mapping is missing.")
            continue
        gene_symbol = gene_cur_gene.gene_symbol
        chromosome = gene_cur_gene.chromosome
        gene_ensembl_id = gene_cur_gene.gene_ensembl_id

        if event == 'mutually_exclusive_exons':    
            exon1_prev = get_exon(exons1_by_rank, rank1 - 1)
            exon1_next = get_exon(exons1_by_rank, rank1 + 1)
            exon2_1_prev = get_exon(exons2_by_rank, rank2 - 1)
            exon2_1_next = get_exon(exons2_by_rank, rank2 + 1)
            if exon1_prev is None or exon1_next is None or exon2_1_prev is None or exon2_1_next is None:
                logger.warning(f"Skipping mutually_exclusive_exons event for {gene_id} because neighboring exons are missing.")
                continue
            points1_1 = sorted([exon1.genomic_start_tx, exon1.genomic_end_tx, exon1_prev.genomic_start_tx, exon1_prev.genomic_end_tx])
            #points1_2 = sorted([exon1.genomic_start_tx, exon1.genomic_end_tx, exon1_next.genomic_start_tx, exon1_next.genomic_end_tx])
            points2_1 = sorted([exon2_1.genomic_start_tx, exon2_1.genomic_end_tx, exon2_1_prev.genomic_start_tx, exon2_1_prev.genomic_end_tx])
            #points2_2 = sorted([exon2_1.genomic_start_tx, exon2_1.genomic_end_tx, exon2_1_next.genomic_start_tx, exon2_1_next.genomic_end_tx])
            junction1 = (points1_1[1], points1_1[2])
            #junction3 = (points1_2[1], points1_2[2])
            junction2 = (points2_1[1], points2_1[2])
            #junction4 = (points2_2[1], points2_2[2])
        elif event == 'alt_prime':
            exon1_prev = get_exon(exons1_by_rank, rank1 - 1)
            exon2_prev = get_exon(exons2_by_rank, rank2 - 1)
            if exon1_prev is None or exon2_prev is None:
                logger.warning(f"Skipping alt_prime event for {gene_id} because the previous exon is missing.")
                continue
            points1 = sorted([exon1.genomic_start_tx, exon1.genomic_end_tx, exon1_prev.genomic_start_tx, exon1_prev.genomic_end_tx])
            points2 = sorted([exon2_1.genomic_start_tx, exon2_1.genomic_end_tx, exon2_prev.genomic_start_tx, exon2_prev.genomic_end_tx])
            junction1 = (points1[1], points1[2])
            junction2 = (points2[1], points2[2])
            if junction1 == junction2:
                exon1_next = get_exon(exons1_by_rank, rank1 + 1)
                exon2_next = get_exon(exons2_by_rank, rank2 + 1)
                if exon1_next is None or exon2_next is None:
                    logger.warning(f"Skipping alt_prime event for {gene_id} because the next exon is missing.")
                    continue
                points1_next = sorted([exon1.genomic_start_tx, exon1.genomic_end_tx, exon1_next.genomic_start_tx, exon1_next.genomic_end_tx])
                points2_next = sorted([exon2_1.genomic_start_tx, exon2_1.genomic_end_tx, exon2_next.genomic_start_tx, exon2_next.genomic_end_tx])
                junction1_next = (points1_next[1], points1_next[2])
                junction2_next = (points2_next[1], points2_next[2])
                if junction1_next == junction2_next:    
                    logger.warning(f"Skipping alt_prime event for {gene_id} because the resulting junctions are the same.")            
        elif event == 'skip_exon':
            exon1_next = get_exon(exons1_by_rank, rank1 + 1)
            if exon1_next is None or exon2_2 is None:
                logger.warning(f"Skipping skip_exon event for {gene_id} because required exons are missing.")
                continue
            points1 = sorted([exon1.genomic_start_tx, exon1.genomic_end_tx, exon1_next.genomic_start_tx, exon1_next.genomic_end_tx])
            points2 = sorted([exon2_1.genomic_start_tx, exon2_1.genomic_end_tx, exon2_2.genomic_start_tx, exon2_2.genomic_end_tx])
            junction1 = (points1[1], points1[2])
            junction2 = (points2[1], points2[2])
        elif event == 'alt_first_exon': 
            exon1_2 = get_exon(exons1_by_rank, rank1 + 1)
            if exon1_2 is None:
                logger.warning(f"Skipping alt_first_exon event for {gene_id} because the second exon is missing.")
                continue
            exon2_2 = get_exon(exons2_by_rank, rank2 + 1)
            if exon2_2 is None:
                logger.warning(f"Skipping alt_first_exon event for {gene_id} because the second exon is missing.")
                continue
            points1 = sorted([exon1.genomic_start_tx, exon1.genomic_end_tx, exon1_2.genomic_start_tx, exon1_2.genomic_end_tx])
            points2 = sorted([exon2_1.genomic_start_tx, exon2_1.genomic_end_tx, exon2_2.genomic_start_tx, exon2_2.genomic_end_tx])
            junction1 = (points1[1], points1[2])
            junction2 = (points2[1], points2[2])
            if junction1 == junction2:
                logger.warning(f"Skipping alt_first_exon event for {gene_id} because the resulting junctions are the same.")
                continue
        elif event == 'alt_last_exon':
            exon1_2 = get_exon(exons1_by_rank, rank1 - 1)
            exon2_2 = get_exon(exons2_by_rank, rank2 - 1)
            if exon1_2 is None or exon2_2 is None:
                logger.warning(f"Skipping alt_last_exon event for {gene_id} because the previous exon is missing.")
                continue
            points1 = sorted([exon1.genomic_start_tx, exon1.genomic_end_tx, exon1_2.genomic_start_tx, exon1_2.genomic_end_tx])
            points2 = sorted([exon2_1.genomic_start_tx, exon2_1.genomic_end_tx, exon2_2.genomic_start_tx, exon2_2.genomic_end_tx])
            junction1 = (points1[1], points1[2])
            junction2 = (points2[1], points2[2])
            if junction1 == junction2:
                logger.warning(f"Skipping alt_last_exon event for {gene_id} because the resulting junctions are the same.")
                continue
        elif event == 'intron_retention':
            if rank1 == 1 or rank1 == transcript_to_max_rank[t_id1]:   
                logger.warning(f"Skipping intron retention event for {gene_id} because it involves the first or last exon.")
                continue # skip if first or last exon
            continue # TBD: implement intron retention junction creation logic
        else:
            logger.error(f"Unknown event type: {event}.")
            exit(1)
        cluster_name = f"{gene_id}_{event}_{idx}" 
        for junction in [junction1, junction2, junction3, junction4]:
            if junction is not None:
                all_junctions.append([chromosome, junction[0], junction[1], gene_symbol, gene_ensembl_id, cluster_name])
        if idx % 10000 == 0:
            # intermediate save to avoid data loss in case of crash
            logger.info(f"Saving intermediate results after processing {idx+1} events.")    
            df_junctions = pd.DataFrame(all_junctions, columns=get_df_junction_columns()) 
            df_junctions.to_csv('as_events_junctions_intermediate.csv', index=False)
    df_junctions = pd.DataFrame(all_junctions, columns=get_df_junction_columns()) 
    # drop duplicates based on all columns except cluster_name
    df_junctions = df_junctions.drop_duplicates(subset=df_junctions.columns.drop(['cluster_name']), keep='first') 
    df_junctions.to_csv(output_csv, index=False)       

def read_junctions_csv(junctions_csv):
    """Read a plain junctions CSV (chromosome, gene_ensembl_id, start_position,
    end_position, cluster_name, ...) into a DataFrame ready for
    JunctionsAnalysis.analyze_junctions()."""
    return pd.read_csv(junctions_csv)


def analyze_junctions2(con, df_junctions=None, junctions_csv=None, hadas_format=False,
                        output_path='as_events_junctions_analysis.csv',
                        n=0, create_pdf=True, print_genes=None, num_workers=5,
                        use_representative_domains=False):
    """Read junctions (from a DataFrame, a plain CSV, or a hadas-format Excel file - exactly
    one of df_junctions/junctions_csv must be given) and run JunctionsAnalysis.analyze_junctions()."""
    if df_junctions is None and junctions_csv is None:
        raise ValueError("Either df_junctions or junctions_csv must be provided.")
    if df_junctions is not None and junctions_csv is not None:
        raise ValueError("Only one of df_junctions or junctions_csv should be provided.")
    if junctions_csv is not None:
        df_junctions = hadas_read_input_file(con, junctions_csv) if hadas_format else read_junctions_csv(junctions_csv)

    analyzer = JunctionsAnalysis(con, logger_instance=logger)
    return analyzer.analyze_junctions(df_junctions=df_junctions, output_path=output_path,
                                      filter_transcript_count=n, create_pdf=create_pdf, print_genes=print_genes, num_workers=num_workers,
                                      use_representative_domains=use_representative_domains)

def analyze_ioe_file(con, ioe_file, output_csv, num_workers=5, use_representative_domains=False):
    df_junctions = utils.ioe2junctions(ioe_file)
    gene_symbols_dict = utils.get_gene_symbols(con, df_junctions.gene_ensembl_id.unique().tolist())
    # add gene symbols to df_junctions
    df_junctions['gene_symbol'] = df_junctions['gene_ensembl_id'].map(gene_symbols_dict)
    analyze_junctions2(con, df_junctions=df_junctions, output_path=output_csv, create_pdf=False, num_workers=num_workers,
                        use_representative_domains=use_representative_domains)


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

def analyze_ioe_files(con, input_path, pattern, output_csv, examples_per_event=5, num_workers=5,
                       use_representative_domains=False):
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
        analyze_junctions2(con, df_junctions=df_examples, output_path=output_csv, create_pdf=False, num_workers=num_workers,
                            use_representative_domains=use_representative_domains)
    else:
        df_all_junctions.to_csv('ioe_all_junctions.csv', index=False)
        analyze_junctions2(con, df_junctions=df_all_junctions, output_path=output_csv, create_pdf=False, num_workers=num_workers,
                            use_representative_domains=use_representative_domains)

def analyze_hadas_input(con, input_file, output_csv, print_genes=None, num_workers=5,
                         use_representative_domains=False, create_pdf=True):
    df_junctions = hadas_read_input_file(con, input_file)
    analyze_junctions2(con, df_junctions=df_junctions, output_path=output_csv, create_pdf=create_pdf, print_genes=print_genes,
                        num_workers=num_workers, use_representative_domains=use_representative_domains)
    return
    gene_ids = df_junctions['gene_ensembl_id'].unique().tolist()
    gene_transcript_counts = utils.get_genes_number_of_transcripts(con, gene_ids)
    df_junctions['num_transcripts'] = df_junctions['gene_ensembl_id'].map(gene_transcript_counts)
    df_junctions.to_csv('hadas_junctions_with_transcript_counts.csv', index=False)
    df_junctions_filtered = df_junctions[df_junctions['num_transcripts'] <= 2]
    df_junctions_filtered.to_csv('hadas_junctions_filtered.csv', index=False)
    analyze_junctions2(con, df_junctions=df_junctions_filtered, output_path=output_csv, create_pdf=True, print_genes=print_genes)
    

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(message)s')
    logger.info("hello")
    #input_file = 'clusters_sum_table_H_vs_M_HN6.xlsx'
    #input_file = '/gpfs0/tals/users/melchio/clusters_sum_table_HN6.xlsx'
    #input_file = '/gpfs0/tals/users/melchio/short_H_vs_M_HN6.xlsx'
    input_file = 'short_H_vs_M_HN6.xlsx'
    #input_file = "hadas_prefered.xlsx"
    #print_genes = ['PFDN5', 'CD6','HNRNPH3','USP16', 'DOCK8', 'HNRNPK']
    print_genes = None
    #output_file = 'clusters_with_3_transcripts.csv'
    #dochap_path = '/gpfs0/tals/projects/Analysis/ariel/DoChap/DoChaP-db/xdDB_merged.sqlite'
    dochap_path = '/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-web/DB_merged.sqlite'
    con = sqlite3.connect(dochap_path)
    
    analyze_hadas_input(con, input_file, 'hadas_junctions_analysis.csv', print_genes=print_genes)
    exit(0)
    #analyze_ioe_files(con, '/Users/arielmelchior/documents/projects/DOMAS/external_data/', 'output_prefix_.*_strict.ioe', 'output_csv.csv', examples_per_event=0)
    analyze_ioe_files(con, '/gpfs0/tals/projects/Analysis/ariel2/DOMAS/external_data/', 'output_prefix_.*_strict.ioe', 'output_csv.csv', examples_per_event=0)
    #ioe_file = '/gpfs0/tals/projects/Analysis/ariel2/DOMAS/external_data/output_prefix_MX_strict.ioe'
    #analyze_ioe_file(con, ioe_file, 'ioe_mx_junctions_analysis.csv')
    exit(0 )
    get_all_as_events(con, 5)
    df_events = pd.read_csv('all_as_events.csv')
    # keep 2 rows oer each event type
    df_events = df_events[df_events['event'] != 'no_canonical_transcript'] # skip no_canonical_transcript for now since it doesn't fit well with junction analysis
    df_events = df_events.groupby('event').head(2)
    create_events_junctions(con, df_events, 'as_events_junctions2.csv')
    analyze_junctions2(con, junctions_csv='as_events_junctions2.csv', output_path='as_events_junctions_analysis2.csv')
    exit(0)
    #get_all_as_events(con)
    #main(con, 'PUF60')
    results = []
    get_gene_alternative_splicing_events_by_name(con, 'FGFR2', results)
    df_results = pd.DataFrame({
        'transcript1' : pd.Series([], dtype='str'),
        'transcript2' : pd.Series([], dtype='str'),
        'event' : pd.Series([], dtype='str'),
        'rank1' : pd.Series([], dtype='int64'),
        'rank2' : pd.Series([], dtype='int64'),
        'rank2_2' : pd.Series([], dtype='int64')
    })
    new_rows_df = pd.DataFrame(results)
    df_results = pd.concat([df_results, new_rows_df], ignore_index=True)
    df_results.to_csv('FGFR2_as_events.csv', index=False)
    exit(0)
    '''
    sort_transcripts_by_exon_count(con, 
                                   'transcripts_sorted_by_exon_count.csv', 
                                   only_canonical=True)
    df = pd.read_csv('transcripts_sorted_by_exon_count.csv')
    multiple_tx_genes = get_genes_with_multiple_transcripts(con)
    my_df = df[df.exon_count.isin([3,4])] # start with 3 exons
    my_df = my_df[my_df.gene_ensembl_id.str.startswith('ENSG')] # human
    my_df = my_df[my_df.gene_ensembl_id.isin(multiple_tx_genes)] # genes with multiple transcripts
    my_df = my_df[my_df['protein_ensembl_id'].notna()] # proper protein
    my_df = my_df[my_df['protein_refseq_id'].notna()] 
    gene_ids = my_df.gene_ensembl_id.values.tolist()
    logger.info(f"Selected {len(gene_ids)} genes for AS event analysis.")
    #gene_ids = gene_ids[:3]
    export_genes_to_gtf(con, gene_ids, 'dochap.gtf')
    os.system("python ../SUPPA/suppa.py generateEvents -i dochap.gtf -o supa_dochap -f ioe -e SE SS MX RI FL")
    as_event_df = combine_and_expand_ioe('supa_dochap_*.ioe')
    as_event_df.to_csv('combined_as_events.csv', index=False)
    '''
    # to do
    # sort transcripts by exon count
    # export 3 transcripts with the least exon count, to gtf format
    # install suppa2
    # Run suppa2 on the gtf
    






