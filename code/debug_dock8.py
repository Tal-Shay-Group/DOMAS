"""Debug script: run the normal pipeline but print domain-window details for DOCK8."""
import logging
import sqlite3
import pandas as pd
import junction_analisys as ja
import alternative_splicing

logging.basicConfig(level=logging.WARNING, format='%(levelname)s:%(message)s')

dochap_path = '/Users/arielmelchior/Documents/projects/DOMAS/DB_merged.sqlite'
con = sqlite3.connect(dochap_path)

# Patch find_relevant_domain_windows to print for DOCK8
_orig_frw = ja.find_relevant_domain_windows

def _debug_frw(transcript_exons, domain_lookup, canonical_transcript_id, transcript_id,
               canonical_junctions, transcript_junctions, junctions):
    t_domains, c_domains = _orig_frw(
        transcript_exons, domain_lookup, canonical_transcript_id, transcript_id,
        canonical_junctions, transcript_junctions, junctions
    )
    # We'll detect DOCK8 by checking if the transcript IDs belong to DOCK8
    # (we check the domain data for clues, or just always print)
    print(f"\n=== find_relevant_domain_windows: canonical={canonical_transcript_id}, transcript={transcript_id} ===")
    print(f"  canonical junctions: {[junctions[i] for i in canonical_junctions]}")
    print(f"  transcript junctions: {[junctions[i] for i in transcript_junctions]}")
    print(f"  t_domains (alt transcript, {len(t_domains)} rows):")
    if t_domains.empty:
        print("    (empty)")
    else:
        print(t_domains[['AA_start','AA_end','length','pfam','interpro','cdd','smart']].to_string())
    print(f"  c_domains (canonical, {len(c_domains)} rows):")
    if c_domains.empty:
        print("    (empty)")
    else:
        print(c_domains[['AA_start','AA_end','length','pfam','interpro','cdd','smart']].to_string())
    return t_domains, c_domains

ja.find_relevant_domain_windows = _debug_frw

# Run only for DOCK8
input_file = '/Users/arielmelchior/Documents/projects/DOMAS/hadas_prefered.xlsx'
df_junctions = alternative_splicing.hadas_read_input_file(con, input_file)
df_dock8 = df_junctions[df_junctions['gene_symbol'] == 'DOCK8']
print(f"DOCK8 junctions rows: {len(df_dock8)}")
print(df_dock8.to_string())

analyzer = ja.JunctionsAnalysis(con)
analyzer.analyze_junctions(df_junctions=df_dock8, output_path='/tmp/dock8_debug.csv', create_pdf=False)
