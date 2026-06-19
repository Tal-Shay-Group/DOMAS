"""
Analyze how each cross-reference source (interpro, pfam, cdd, smart, tigr,
CDD_id) contributes to domain annotations, for both:
  - DomainType: the domain *definitions*
  - DomainEvent: the domain *instances* (one row per occurrence in a protein)

For each source, reports:
  - coverage: how often the source has a value at all
  - unique: how often it's the *only* source with a value
  - multi: how often the source field itself lists more than one ID
  - winner: how often it's chosen as the display name (DOMAIN_NAME_PREFIX_PRIORITY)

Also reports a pairwise co-occurrence matrix (how often two sources are
both present on the same row).
"""

import sqlite3

import pandas as pd

SOURCE_COLUMNS = ['interpro', 'pfam', 'cdd', 'smart', 'tigr', 'CDD_id']
PREFIX_PRIORITY = ['IPR', 'pfam', 'cd', 'smart', 'tigr', 'CDD']


def _first_id(value):
    return str(value).split(';')[0].strip()


def source_presence(df, cols=SOURCE_COLUMNS):
    """Boolean DataFrame: True where each source column has a non-null value."""
    return pd.DataFrame({col: df[col].notna() & (df[col].astype(str).str.strip() != '') for col in cols})


def winning_source(df, presence, cols=SOURCE_COLUMNS, prefixes=PREFIX_PRIORITY):
    """Series of the winning source column name per row (or None if no match), vectorized."""
    first_ids = pd.DataFrame({col: df[col].astype(str).str.split(';').str[0].str.strip().str.lower() for col in cols})
    match = pd.DataFrame({
        col: presence[col] & first_ids[col].str.startswith(prefix.lower())
        for col, prefix in zip(cols, prefixes)
    })
    any_match = match.any(axis=1)
    winner = match.idxmax(axis=1)
    return winner.where(any_match, None)


def summarize(df, label, cols=SOURCE_COLUMNS):
    n = len(df)
    presence = source_presence(df, cols)
    num_sources = presence.sum(axis=1)
    winners = winning_source(df, presence, cols)

    rows = []
    for col in cols:
        coverage = presence[col].sum()
        unique = ((presence[col]) & (num_sources == 1)).sum()
        multi = (df[col].astype(str).str.contains(';')).sum()
        winner = (winners == col).sum()
        rows.append({
            'source': col,
            'coverage_n': coverage,
            'coverage_%': 100 * coverage / n,
            'unique_n': unique,
            'unique_%': 100 * unique / n,
            'multi_n': multi,
            'multi_%': 100 * multi / n if coverage else 0.0,
            'winner_n': winner,
            'winner_%': 100 * winner / n,
        })

    summary = pd.DataFrame(rows).set_index('source')

    cooccurrence = pd.DataFrame(
        [[ (presence[a] & presence[b]).sum() for b in cols] for a in cols],
        index=cols, columns=cols,
    )

    either_n = pd.DataFrame(
        [[ (presence[a] | presence[b]).sum() for b in cols] for a in cols],
        index=cols, columns=cols,
    )
    either_pct = 100 * either_n / n

    print(f"\n=== {label} (n={n}) ===")
    print(summary.round(2).to_string())
    print("\nCo-occurrence (rows & cols both present):")
    print(cooccurrence.to_string())
    print("\nEither present (rows | cols), count:")
    print(either_n.to_string())
    print("\nEither present (rows | cols), %:")
    print(either_pct.round(2).to_string())

    return summary, cooccurrence, either_n


SPECIES = ['D_rerio', 'H_sapiens', 'M_musculus', 'R_norvegicus', 'X_tropicalis']

# (label, filter applied to df_t, a DataFrame of Transcripts joined with Genes)
SCOPES = [("Entire DB", None)]
for _specie in SPECIES:
    SCOPES.append((_specie, lambda df_t, _specie=_specie: df_t['specie'] == _specie))
SCOPES += [
    ("H_sapiens canonical transcripts", lambda df_t: (df_t['specie'] == 'H_sapiens') & (df_t['canonical'] != 0)),
    ("H_sapiens non-canonical transcripts", lambda df_t: (df_t['specie'] == 'H_sapiens') & (df_t['canonical'] == 0)),
]

EVENTS_QUERY = """
    SELECT {cols}, e.type_id, e.protein_refseq_id, e.protein_ensembl_id,
           e.AA_start, e.AA_end
    FROM DomainEvent e
    JOIN DomainType dt ON e.type_id = dt.type_id
"""

TRANSCRIPTS_QUERY = """
    SELECT t.protein_refseq_id, t.protein_ensembl_id, t.canonical, g.specie
    FROM Transcripts t
    JOIN Genes g ON t.gene_ensembl_id = g.gene_ensembl_id
"""


def _events_contained_in(reference_df, other_df, tolerance=2):
    """
    Return the subset of other_df whose events are positionally contained
    within any event in reference_df on the same protein (±tolerance aa).

    Containment: other.AA_start >= ref.AA_start - tol
                 AND other.AA_end <= ref.AA_end + tol
    """
    protein_cols = ['protein_refseq_id', 'protein_ensembl_id']
    contained_idx = []

    ref_grouped = reference_df.groupby(protein_cols, dropna=False)

    for key, other_group in other_df.groupby(protein_cols, dropna=False):
        if key not in ref_grouped.groups:
            continue
        ref_group = ref_grouped.get_group(key)
        ref_starts = ref_group['AA_start'].values
        ref_ends = ref_group['AA_end'].values

        for row_idx, row in other_group.iterrows():
            if any(
                row['AA_start'] >= rs - tolerance and row['AA_end'] <= re + tolerance
                for rs, re in zip(ref_starts, ref_ends)
            ):
                contained_idx.append(row_idx)

    return other_df.loc[contained_idx]


def source_contribution_table(df_events, tolerance=2, cols=SOURCE_COLUMNS):
    """
    For each source X:
      1. Take source X events (rows where X is non-null).
      2. From all other events (X is null), remove those positionally
         contained within any source X event on the same protein (±tolerance aa).
      3. contribution_n  = len(source X events)
         effective_total = source X events + non-contained other events
         contribution_%  = contribution_n / effective_total * 100

    Returns a DataFrame with one row per source.
    """
    presence = source_presence(df_events, cols)
    total_n = len(df_events)
    rows = []

    for col in cols:
        src_mask = presence[col]
        src_events = df_events[src_mask]
        other_events = df_events[~src_mask]

        contained = _events_contained_in(src_events, other_events, tolerance)
        non_contained_others = other_events.drop(index=contained.index)

        effective_total = len(src_events) + len(non_contained_others)
        contribution_n = len(src_events)
        subsumed_n = len(contained)

        rows.append({
            'source': col,
            'source_events_n': contribution_n,
            'source_events_%_of_total': round(100 * contribution_n / total_n, 2),
            'subsumed_from_others_n': subsumed_n,
            'effective_total_n': effective_total,
            'contribution_%': round(100 * contribution_n / effective_total, 2) if effective_total else 0.0,
        })

    return pd.DataFrame(rows).set_index('source')


def analyze_scope(df_events_all, df_types_all, df_t, label, filter_fn):
    if filter_fn is None:
        df_events = df_events_all
    else:
        scoped_t = df_t[filter_fn(df_t)]
        refseq_ids = set(scoped_t['protein_refseq_id'].dropna())
        ensembl_ids = set(scoped_t['protein_ensembl_id'].dropna())
        df_events = df_events_all[
            df_events_all['protein_refseq_id'].isin(refseq_ids)
            | df_events_all['protein_ensembl_id'].isin(ensembl_ids)
        ]

    summarize(df_events, f"{label} - DomainEvent (instances)")

    contrib = source_contribution_table(df_events)
    print(f"\n--- {label} - Source Contribution (±2 aa containment removed) ---")
    print(contrib.to_string())

    type_ids = set(df_events['type_id'].unique())
    df_types = df_types_all[df_types_all['type_id'].isin(type_ids)]
    summarize(df_types, f"{label} - DomainType (definitions used)")


def main(db_path):
    con = sqlite3.connect(db_path)

    cols = ', '.join('dt.' + c for c in SOURCE_COLUMNS)
    df_events_all = pd.read_sql_query(EVENTS_QUERY.format(cols=cols), con)
    df_types_all = pd.read_sql_query(f"SELECT type_id, {', '.join(SOURCE_COLUMNS)} FROM DomainType", con)
    df_t = pd.read_sql_query(TRANSCRIPTS_QUERY, con)

    for label, filter_fn in SCOPES:
        print(f"\n{'#' * 70}\n# {label}\n{'#' * 70}")
        analyze_scope(df_events_all, df_types_all, df_t, label, filter_fn)


if __name__ == '__main__':
    main('/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-web/DB_merged.sqlite')
