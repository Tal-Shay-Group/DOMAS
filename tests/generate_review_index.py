"""Derive review_cases/INDEX.md and cases.json from the current reference outputs.

The picks are properties of the *output* - "a cluster where is_most_like_canonical
is set", "one containing a non-coding transcript" - so they stop being true whenever
the outputs move, which a refreshed fixture or a rebuilt DoChaP both do. Deriving
them keeps the index honest without anyone having to notice.

Writes two files:

  review_cases/INDEX.md    for reading
  review_cases/cases.json  the same picks, machine-readable, so
                           collect_review_pdfs.py needs no list of its own

Run after the references exist:

    python3 tests/generate_reference_outputs.py
    python3 tests/generate_review_pdfs.py
    python3 tests/generate_review_index.py
"""
import datetime
import glob
import json
import os
import sqlite3
import subprocess
import sys

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(TESTS_DIR, '..'))
sys.path.insert(0, os.path.join(REPO, 'code'))
sys.path.insert(0, TESTS_DIR)

import pandas as pd  # noqa: E402
import utils  # noqa: E402
from junction_analisys import (  # noqa: E402
    filter_representative_domains, _PRIMARY_ENTRY_TYPES, _SITE_ENTRY_TYPES,
)
from conftest import DEFAULT_DB_PATH  # noqa: E402

REF = os.path.join(REPO, 'tests', 'reference_outputs')
CASES = os.path.join(REPO, 'review_cases')
# Every reference case exists in both variants. `restrict_True` draws only the
# canonical transcript and the ones compared to it, which is what a case about domain
# resolution or intron retention wants. The A series needs `restrict_False`: its
# subject is which candidate was chosen, so the losing candidates must be visible.
IOE_CASE = 'ioe_csv__restrict_True__representative_True'
CAT_CASE = 'category_examples__restrict_True__representative_True'
IOE_CASE_ALL = 'ioe_csv__restrict_False__representative_True'
CAT_CASE_ALL = 'category_examples__restrict_False__representative_True'
# Live under review_cases/, not tests/reference_outputs/ - see cluster_facts(root=).
# Same split: the A pool draws every transcript, everything else only the compared.
A_POOL = 'leafcutter_all_transcripts'
LADDER_POOL = 'leafcutter'

SKIP_EVENTS = {
    'transcript_doesnt_have_features', 'no_unique_features', 'gene_not_in_db',
    'no_gene_specified', 'no_canonical_transcript', 'only_one_transcript',
    'no_canonical_features', 'feature_not_mapped',
}
S4_OUTCOMES = [
    'dropped domain', 'added_domain', 'longer', 'same', 'shorter', 'same_domains',
    'merged domain', 'split domain', 'increased_domain_number',
    'reduced_domain_number', 'shorter_domains', 'longer_domains',
]


# --- facts about the current outputs ---------------------------------------

def noncoding_map(con):
    """transcript id -> True when it carries no protein."""
    df = pd.read_sql_query(
        "SELECT transcript_ensembl_id t, transcript_refseq_id r, "
        "protein_ensembl_id p1, protein_refseq_id p2 FROM Transcripts", con)

    def blank(s):
        return s.isna() | s.astype(str).str.strip().isin(['', 'None', 'nan'])

    df['nc'] = blank(df.p1) & blank(df.p2)
    out = {}
    for row in df.itertuples():
        for key in (row.t, row.r):
            if key and str(key) not in ('nan', 'None', ''):
                out[str(key)] = bool(row.nc)
    return out


def gene_transcript_counts(con):
    """gene symbol -> how many transcripts DoChaP holds for it.

    The count the PDF title shows, and the one deciding whether a case is worth
    looking at: a two-transcript gene has at most one comparable candidate, so both
    flags fire on it unopposed and the case demonstrates nothing.
    """
    df = pd.read_sql_query(
        "SELECT g.gene_symbol, COUNT(*) n FROM Genes g JOIN Transcripts t "
        "  ON t.gene_ensembl_id = g.gene_ensembl_id OR t.gene_GeneID_id = g.gene_GeneID_id "
        "WHERE g.specie = 'H_sapiens' GROUP BY g.gene_symbol", con)
    return dict(zip(df.gene_symbol, df.n))


def gene_strands(con):
    """gene_ensembl_id -> strand. The D series is one case per (event type, strand),
    and the strand is not in the results - it is read from the database per gene."""
    return dict(pd.read_sql_query(
        "SELECT gene_ensembl_id, strand FROM Genes WHERE specie = 'H_sapiens'", con).values)


def cds_lengths(con, transcript_ids):
    """transcript id -> coding length in bases, the same quantity analyze() ranks on
    (the largest CDS-relative exon offset, not the genomic cds_end - cds_start span)."""
    ids = [t for t in transcript_ids if isinstance(t, str)]
    if not ids:
        return {}
    marks = ','.join('?' * len(ids))
    df = pd.read_sql_query(
        f"SELECT transcript_ensembl_id e, transcript_refseq_id r, MAX(abs_end_CDS) n "
        f"FROM Transcript_Exon "
        f"WHERE transcript_ensembl_id IN ({marks}) OR transcript_refseq_id IN ({marks}) "
        f"GROUP BY transcript_ensembl_id, transcript_refseq_id", con, params=ids + ids)
    out = {}
    for row in df.itertuples():
        for key in (row.e, row.r):
            if isinstance(key, str) and key in set(ids):
                out[key] = row.n
    return out


def cluster_facts(case, nc, root=REF, con=None, n_tx_of=None):
    """One row per cluster of `case`, with the properties the picks are chosen on.

    `root` is tests/reference_outputs for the reference cases and review_cases for
    the sets generate_review_pdfs.py builds; the A series draws from the latter, no
    reference fixture holding a gene with 3-5 transcripts.
    """
    path = os.path.join(root, case, 'results.csv')
    if not os.path.exists(path):
        return pd.DataFrame()
    df = pd.read_csv(path, low_memory=False)
    cds = cds_lengths(con, df.transcript_id.dropna().unique()) if con is not None else {}
    rows = []
    for cluster, group in df.groupby('cluster'):
        transcripts = group.transcript_id.dropna().unique()
        compared = group.loc[~group.event_type.isin(SKIP_EVENTS)]
        comparable = compared.transcript_id.dropna().unique()
        noncoding = [t for t in comparable if nc.get(str(t), False)]
        longest = compared.loc[compared.is_longest_cds == True, 'transcript_id'].dropna().unique()  # noqa: E712
        most = compared.loc[compared.is_most_like_canonical == True, 'transcript_id'].dropna().unique()  # noqa: E712
        # The candidate holding the longest CDS regardless of whether it codes -
        # when that is a non-coding transcript, step 1 of the priority is what
        # keeps the tag off it, which is the only way to see step 1 work.
        by_cds = sorted(comparable, key=lambda t: (cds.get(t, -1), t), reverse=True)
        rows.append({
            'case': case,
            'cluster': cluster,
            'gene': group.gene_symbol.iat[0],
            'n_tx': len(transcripts),
            'gene_n_tx': (n_tx_of or {}).get(group.gene_symbol.iat[0]),
            'n_comparable': len(comparable),
            'n_noncoding': sum(1 for t in transcripts if nc.get(str(t), False)),
            'nc_comparable': len(noncoding),
            'coding_comparable': len(comparable) - len(noncoding),
            'longest_id': longest[0] if len(longest) else None,
            'most_like_id': most[0] if len(most) else None,
            'rules_differ': bool(len(most) and len(longest) and most[0] != longest[0]),
            'nc_holds_longest_cds': bool(by_cds and nc.get(str(by_cds[0]), False)),
            'most_like': bool((group.is_most_like_canonical == True).any()),  # noqa: E712
            'outcomes': sorted(set(group.event_type) - SKIP_EVENTS),
            'as_type': str(cluster).split('_')[0],
        })
    return pd.DataFrame(rows)


# A duplicate-collapse removal reads two ways depending on how much the entries
# share. Above this fraction of the shorter one they are two calls of one region and
# collapsing is right; below it they are neighbouring repeat instances, and one is
# being lost.
TANDEM_OVERLAP = 0.20


def _duplicate_removals(raw, kept_idx):
    """Rows of `raw` the duplicate-collapse step removed, with how far each overlaps
    the kept entry of the same accession that displaced it.

    Identified from the filter's own output rather than by re-running its coverage
    phase: a row absent from the output while a LONGER kept row of the same accession
    overlaps it is exactly what that step removes.
    """
    out = []
    for i in raw.index:
        if i in kept_idx:
            continue
        for j in kept_idx:
            if raw.domain_id[j] != raw.domain_id[i]:
                continue
            overlap = min(raw.AA_end[i], raw.AA_end[j]) - max(raw.AA_start[i], raw.AA_start[j]) + 1
            if overlap <= 0:
                continue
            shorter = min(raw.AA_end[i] - raw.AA_start[i], raw.AA_end[j] - raw.AA_start[j]) + 1
            out.append({
                'domain_id': raw.domain_id[i],
                'dropped': f'{raw.AA_start[i]}-{raw.AA_end[i]}',
                'kept': f'{raw.AA_start[j]}-{raw.AA_end[j]}',
                'overlap_aa': int(overlap),
                'overlap_frac': overlap / shorter,
            })
            break
    return out


def domain_ladder_facts(con, pools):
    """Which branch of the InterPro tier ladder each transcript exercises.

    `pools` is a list of (case, root) pairs: the duplicate-collapse behaviours occur
    only in the LeafCutter review set, so one results.csv does not reach every case.
    """
    found = {}
    for case, root in pools:
        path = os.path.join(root, case, 'results.csv')
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, low_memory=False)
        ids = set(df.transcript_id.dropna()) | set(df.canonical_transcript_id.dropna())
        domains = utils.get_domains_db(con, ids, use_representative_domains=True)
        column = 'transcript_ensembl_id_version'
        gene_of = {}
        for row in df.itertuples():
            for t in (row.transcript_id, row.canonical_transcript_id):
                if isinstance(t, str):
                    gene_of.setdefault(t, row.gene_symbol)

        for t in sorted(ids):
            raw = domains[domains[column].astype(str).str.startswith(str(t).split('.')[0])]
            if raw.empty or 'type' not in raw.columns:
                continue
            kept = filter_representative_domains(raw)
            entry_type = raw['type'].fillna('')
            is_ipr = raw['domain_id'].astype(str).str.startswith('IPR')
            primary = int((is_ipr & entry_type.isin(_PRIMARY_ENTRY_TYPES)).sum())
            site = int((is_ipr & entry_type.isin(_SITE_ENTRY_TYPES)).sum())
            member = int((~is_ipr).sum())
            parent = int((is_ipr & ~entry_type.isin(_PRIMARY_ENTRY_TYPES)
                          & ~entry_type.isin(_SITE_ENTRY_TYPES)).sum())
            fact = {'gene': gene_of.get(t, '?'), 'transcript': t,
                    'raw_n': len(raw), 'kept_n': len(kept), 'case': case, 'root': root}
            if len(kept) < len(raw):
                found.setdefault('reduction', []).append(fact)
            if primary == 0 and member > 0 and len(kept):
                found.setdefault('member_only', []).append(fact)
            if primary and parent and member:
                found.setdefault('all_tiers', []).append(fact)
            if site:
                found.setdefault('site_ptm', []).append(fact)
            if raw['domain_id'].duplicated().any():
                found.setdefault('repeated', []).append(fact)
            for removal in _duplicate_removals(raw, set(kept.index)):
                key = ('tandem_collapsed' if removal['overlap_frac'] < TANDEM_OVERLAP
                       else 'duplicate_redundant')
                found.setdefault(key, []).append(dict(fact, **removal))
    return found


# --- choosing the picks ----------------------------------------------------

def smallest(frame, **conditions):
    """The fewest-transcript cluster satisfying every condition, or None."""
    sub = frame
    for column, wanted in conditions.items():
        sub = sub[sub[column] == wanted] if isinstance(wanted, bool) else sub[sub[column] >= wanted]
    sub = sub[sub.n_comparable >= 1]
    if sub.empty:
        return None
    return sub.sort_values(['n_tx', 'gene']).iloc[0].to_dict()


# A gene with one or two transcripts cannot show the selection working: a single
# comparable candidate takes both flags unopposed, and "longest CDS" means longest of
# one. This band is small enough to read at a sitting and large enough to choose in.
A_MIN_GENE_TX, A_MAX_GENE_TX = 3, 5


def richest(frame, **conditions):
    """The most-candidates cluster in the A band satisfying every condition, or None.

    Ties break on the smaller gene (fewer transcripts to read through) and then on
    gene symbol, so the pick does not depend on row order.
    """
    sub = frame[frame.gene_n_tx.between(A_MIN_GENE_TX, A_MAX_GENE_TX)
                & (frame.n_comparable >= 2)]
    for column, wanted in conditions.items():
        sub = sub[sub[column] == wanted] if isinstance(wanted, bool) else sub[sub[column] >= wanted]
    if sub.empty:
        return None
    return sub.sort_values(['n_comparable', 'gene_n_tx', 'gene'],
                           ascending=[False, True, True]).iloc[0].to_dict()


def build_cases(con):
    nc = noncoding_map(con)
    n_tx_of = gene_transcript_counts(con)
    ioe = cluster_facts(IOE_CASE, nc, con=con, n_tx_of=n_tx_of)
    cat = cluster_facts(CAT_CASE, nc, con=con, n_tx_of=n_tx_of)
    both = pd.concat([ioe, cat], ignore_index=True)
    # The same clusters, attributed to the every-transcript variant of each case.
    # Restricting changes which transcripts are DRAWN, never the analysis, so both
    # variants share a results.csv; only the case name differs, and that is what
    # resolves the directory a pick is read from.
    both_all = both.assign(case=both.case.map({IOE_CASE: IOE_CASE_ALL, CAT_CASE: CAT_CASE_ALL}))
    # The A series draws from the LeafCutter review set: it is the only pool with
    # genes in the A band (see A_MIN_GENE_TX). generate_review_pdfs.py keeps the
    # clusters; which one illustrates which rule is decided here, from the output.
    a_pool = cluster_facts(A_POOL, nc, root=CASES, con=con, n_tx_of=n_tx_of)
    if a_pool.empty:
        a_pool = both_all
    ladder = domain_ladder_facts(con, [(IOE_CASE, REF), (CAT_CASE, REF), (LADDER_POOL, CASES)])

    cases = []

    def add(code, label, pick, shows, source=None, pattern=None):
        if pick is None:
            cases.append({'code': code, 'label': label, 'shows': shows, 'missing': True})
            return
        # A pick from one of the review pools was written by generate_review_pdfs.py
        # into review_cases/, everything else by generate_reference_outputs.py.
        default_source = (os.path.join('review_cases', pick['case'])
                          if pick['case'] in (A_POOL, LADDER_POOL)
                          else os.path.join('tests', 'reference_outputs', pick['case']))
        entry = {
            'code': code, 'label': label, 'shows': shows,
            'gene': pick['gene'], 'cluster': pick['cluster'], 'case': pick['case'],
            'n_tx': int(pick['n_tx']), 'n_comparable': int(pick['n_comparable']),
            'n_noncoding': int(pick['n_noncoding']), 'outcomes': pick['outcomes'],
            'source': source or default_source,
            'pattern': pattern or f"*_{pick['gene']}_*.pdf",
        }
        if pick.get('gene_n_tx') is not None and pd.notna(pick['gene_n_tx']):
            entry['gene_n_tx'] = int(pick['gene_n_tx'])
        for key in ('longest_id', 'most_like_id'):
            if pick.get(key):
                entry[key] = pick[key]
        cases.append(entry)

    # A. choosing the comparable transcript. Each pick needs a real choice, so each
    # asks for at least two comparable candidates and for the property that makes the
    # rule visible - not merely for the flag to be set. With no such cluster in the
    # band, the fewest-transcript pick stands in so the case is not reported missing.
    add('A1', 'most_like_set', richest(a_pool, rules_differ=True) or smallest(both_all, most_like=True),
        'the outside-exon gate picks a DIFFERENT transcript than longest-CDS - '
        'the one case where is_most_like_canonical changes the answer')
    # nc_holds_longest_cds is excluded as A3's subject: with it, both land on the
    # same cluster (the tie-break prefers the smaller gene) and one PDF stands for two
    # different demonstrations.
    add('A2', 'most_like_unset',
        richest(a_pool, most_like=False, nc_holds_longest_cds=False)
        or smallest(both_all, most_like=False),
        'is_most_like_canonical is UNSET - no candidate passes the gate, so longest-CDS '
        'decides among coding candidates alone')
    add('A3', 'has_noncoding', richest(a_pool, nc_holds_longest_cds=True) or smallest(both_all, n_noncoding=1),
        'a non-coding candidate holds the longest CDS - step 1 excludes it from '
        'selection, so the tag goes to a coding transcript instead')
    add('A4', 'several_comparable', smallest(both_all, n_comparable=2),
        'several comparable transcripts, so the flags have to choose between them')

    # B. domain resolution
    def add_ladder(code, key, shows, hit, n_variants, detail=None):
        root_name = ('review_cases' if hit['root'] == CASES
                     else os.path.join('tests', 'reference_outputs'))
        cases.append({
            'code': code, 'label': key,
            'shows': f"{shows} ({hit['raw_n']} -> {hit['kept_n']} annotations)",
            'gene': hit['gene'], 'transcript': hit['transcript'], 'case': hit['case'],
            'source': os.path.join(root_name, hit['case']),
            'pattern': f"*_{hit['gene']}_*.pdf", 'n_variants': n_variants,
            **({'detail': detail} if detail else {}),
        })

    for code, key, shows in (
        ('B1', 'reduction', 'the tier ladder discards a redundant annotation'),
        ('B2', 'member_only', "demote, don't delete - a member-DB hit kept as the only evidence"),
        ('B3', 'all_tiers', 'all three tiers present in one transcript'),
        ('B4', 'site_ptm', 'Site/PTM entries removed - unconditional, unlike the tier rule'),
        ('B5', 'repeated', 'the same accession more than once'),
    ):
        hits = ladder.get(key, [])
        if not hits:
            cases.append({'code': code, 'label': key, 'shows': shows, 'missing': True})
            continue
        add_ladder(code, key, shows, hits[0], len(hits))

    # B6/B7 - the duplicate-collapse step, which no other case reaches. Both order by
    # how far the discarded entry overlaps the one displacing it, from opposite ends:
    # B6 takes the SMALLEST overlaps, neighbouring repeat instances where one is lost;
    # B7 the largest, two calls of one region where discarding one is right. One case
    # per gene - the same span recurs across a gene's transcripts.
    def add_duplicate_slate(prefix, key, hits, shows, count):
        seen = set()
        ranked = [h for h in hits if not (h['gene'] in seen or seen.add(h['gene']))]
        for n, hit in enumerate(ranked[:count]):
            code = f'{prefix}{chr(ord("a") + n)}' if count > 1 else prefix
            add_ladder(code, key, shows, hit, len(hits),
                       detail=(f"{hit['domain_id']} {hit['dropped']} dropped, "
                               f"{hit['kept']} kept, overlap {hit['overlap_aa']} aa "
                               f"({hit['overlap_frac'] * 100:.0f}% of the shorter)"))
        if not ranked:
            cases.append({'code': prefix, 'label': key, 'shows': shows, 'missing': True})

    add_duplicate_slate(
        'B6', 'tandem_collapsed',
        sorted(ladder.get('tandem_collapsed', []),
               key=lambda h: (h['overlap_aa'], h['gene'], h['transcript'])),
        'an instance of a repeat DISCARDED for overlapping its neighbour of the same '
        'accession - the domain count drops though both instances are real', count=5)
    add_duplicate_slate(
        'B7', 'duplicate_redundant',
        sorted(ladder.get('duplicate_redundant', []),
               key=lambda h: (-h['overlap_aa'], h['gene'], h['transcript'])),
        'the same rule applied where it IS right - two calls of one region, '
        'near-identical spans, one discarded', count=3)

    # C. intron retention
    ri = both[(both.as_type == 'RI') & (both.n_comparable >= 1)]
    ri_calls = ri[ri.outcomes.apply(lambda o: any(x not in ('no_domains_in_region',) for x in o))]
    add('C1', 'RI_domain_calls',
        ri_calls.sort_values('n_tx').iloc[0].to_dict() if not ri_calls.empty else None,
        'a retained intron producing real domain calls')
    ri_nc = ri[ri.n_noncoding > 0]
    add('C2', 'RI_noncoding',
        ri_nc.sort_values('n_tx').iloc[0].to_dict() if not ri_nc.empty else None,
        'a retained intron in a cluster with a non-coding transcript')

    # D. one case per input format x AS event type. The gene is read from the run,
    # not named here: generate_review_pdfs.py picks each rMATS cluster on what it
    # resolves, so a named gene would need re-editing whenever that choice moves.
    rmats_facts = cluster_facts('rmats', nc, root=CASES, con=con, n_tx_of=n_tx_of)
    strand_of = gene_strands(con)

    def rmats_pick(as_type, strand):
        for row in rmats_facts.sort_values(['n_comparable', 'gene'],
                                           ascending=[False, True]).itertuples():
            cluster = str(row.cluster)
            gene_id = cluster.split('_')[1] if '_' in cluster else None
            if cluster.split('_')[0] == as_type and strand_of.get(gene_id) == strand:
                return row
        return None

    for code, label, as_type, strand, shows in (
        ('D1a', 'A3SS_plus', 'A3SS', '+', 'plus-strand A3SS'),
        ('D1b', 'A3SS_minus', 'A3SS', '-', 'minus-strand A3SS - the boundary that varies with strand'),
        ('D2a', 'A5SS_plus', 'A5SS', '+', 'plus-strand A5SS'),
        ('D2b', 'A5SS_minus', 'A5SS', '-', 'minus-strand A5SS - the boundary that varies with strand'),
        ('D3a', 'SE_plus', 'SE', '+', 'plus-strand SE'),
        ('D3b', 'SE_minus', 'SE', '-', 'minus-strand SE'),
        ('D4a', 'MXE_plus', 'MXE', '+', 'plus-strand MXE'),
        ('D4b', 'MXE_minus', 'MXE', '-', 'minus-strand MXE'),
    ):
        pick = rmats_pick(as_type, strand)
        gene = pick.gene if pick is not None else None
        n_comparable = int(pick.n_comparable) if pick is not None else 0
        pattern = f'*_{as_type}_{gene}_*.pdf' if gene else f'*_{as_type}_*.pdf'
        note = '' if n_comparable else ' - NO comparable transcript, canonical only'
        entry = {
            'code': code, 'label': label, 'shows': shows + note,
            'source': os.path.join('review_cases', 'rmats'), 'pattern': pattern,
            'gene': gene, 'n_comparable': n_comparable,
            'missing': not (gene and glob.glob(os.path.join(CASES, 'rmats', pattern))),
        }
        if pick is not None:
            entry.update(cluster=str(pick.cluster), n_tx=int(pick.n_tx),
                         outcomes=pick.outcomes)
            if pick.gene_n_tx is not None and pd.notna(pick.gene_n_tx):
                entry['gene_n_tx'] = int(pick.gene_n_tx)
        cases.append(entry)

    cases.append({
        'code': 'D5', 'label': 'majiq_lsv', 'shows': 'MAJIQ per-LSV classification',
        'source': os.path.join('review_cases', 'majiq'), 'pattern': '*.pdf',
        'missing': not glob.glob(os.path.join(CASES, 'majiq', '*.pdf')),
    })

    # The multi-gene LeafCutter cluster, one case per gene it names, read from the run
    # for the same reason. Both halves must resolve a comparison, or the second-gene
    # case has nothing to show.
    lc_facts = cluster_facts(LADDER_POOL, nc, root=CASES, con=con, n_tx_of=n_tx_of)
    multigene = sorted({str(r.cluster).rsplit(':', 1)[0]: None
                        for r in lc_facts.itertuples()
                        if str(r.cluster).count(':') >= 2})
    halves = []
    for base in multigene:
        rows = [r for r in lc_facts.itertuples() if str(r.cluster).startswith(base + ':')]
        if len(rows) >= 2 and all(r.n_comparable >= 1 for r in rows):
            halves = sorted(rows, key=lambda r: str(r.cluster))
            break
    for n, row in enumerate(halves[:2]):
        cases.append({
            'code': f'D6{chr(ord("a") + n)}',
            'label': f'leafcutter_multigene_{row.gene}',
            'shows': ('one LeafCutter event, first of the two genes it names' if n == 0
                      else 'the same event, second of the two genes it names'),
            'source': os.path.join('review_cases', LADDER_POOL),
            'pattern': f'*_{row.gene}_*.pdf', 'gene': row.gene,
            'cluster': str(row.cluster), 'n_comparable': int(row.n_comparable),
            'missing': not glob.glob(os.path.join(CASES, LADDER_POOL, f'*_{row.gene}_*.pdf')),
        })
    if not halves:
        for code in ('D6a', 'D6b'):
            cases.append({'code': code, 'label': 'leafcutter_multigene',
                          'shows': 'a multi-gene LeafCutter event', 'missing': True})

    cases.append({
        'code': 'E1', 'label': 'longer_domains',
        'shows': 'longer_domains - equal counts >1 on both sides, compared side longer',
        'source': os.path.join('review_cases', LADDER_POOL), 'pattern': '*_DHX29_*.pdf',
        'missing': not glob.glob(os.path.join(CASES, LADDER_POOL, '*_DHX29_*.pdf')),
    })

    return cases, both, cat, ladder


# --- rendering -------------------------------------------------------------

def db_build_date(path):
    return datetime.datetime.fromtimestamp(os.path.getmtime(path)).strftime('%Y-%m-%d %H:%M')


def git_commit():
    try:
        return subprocess.run(['git', '-C', REPO, 'rev-parse', '--short', 'HEAD'],
                              capture_output=True, text=True).stdout.strip() or 'unknown'
    except Exception:
        return 'unknown'


def render(cases, clusters, cat, ladder, db_path):
    out = []
    w = out.append
    w('# Review cases — one example per decision the algorithm makes\n')
    w('> **Generated file — do not edit by hand.** Rebuild with')
    w('> `python3 tests/generate_review_index.py`.')
    w('>')
    w('> The picks are properties of the output, so they stop being true when the')
    w('> outputs move — which has happened both when a fixture was refreshed and when')
    w('> DoChaP was rebuilt. Deriving them is what keeps this file honest.\n')
    w(f'- reference outputs as of commit `{git_commit()}`')
    w(f'- DoChaP database built `{db_build_date(db_path)}`')
    w(f'- generated `{datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}`\n')
    w('The PDFs are not in git (`.gitignore` excludes `*.pdf`). Rebuild them with')
    w('`python3 tests/generate_reference_outputs.py` and')
    w('`python3 tests/generate_review_pdfs.py`, then collect them into one directory')
    w('with `python3 tests/collect_review_pdfs.py`.\n')
    w('---\n')

    groups = [('A', 'Choosing the comparable transcript'),
              ('B', 'Domain resolution (InterPro entry-type ladder)'),
              ('C', 'Intron retention'),
              ('D', 'Input format × AS event type'),
              ('E', 'Classification outcomes')]
    for letter, title in groups:
        w(f'## {letter}. {title}\n')
        w('| # | What it shows | Where | Detail |')
        w('| :-- | :--- | :--- | :--- |')
        for case in cases:
            if not case['code'].startswith(letter):
                continue
            if case.get('missing'):
                w(f"| {case['code']} | {case['shows']} | — | **no case available** |")
                continue
            where = f"`{case['pattern']}` in `{case['source']}`"
            bits = []
            if 'cluster' in case:
                bits.append(f"{case['gene']} `{case['cluster']}`")
                # the gene's transcript count, not the number of result rows: what
                # the PDF title shows, and what says whether there was a choice
                n_tx = case.get('gene_n_tx', case.get('n_tx'))
                if n_tx is not None:
                    bits.append(f"{n_tx} tx, {case['n_comparable']} comparable")
                else:
                    bits.append(f"{case['n_comparable']} comparable")
                if case.get('n_noncoding'):
                    bits.append(f"{case['n_noncoding']} non-coding")
                if case.get('longest_id'):
                    chosen = f"longest CDS `{case['longest_id']}`"
                    if case.get('most_like_id'):
                        same = case['most_like_id'] == case['longest_id']
                        chosen += (' (also most like canonical)' if same
                                   else f", most like canonical `{case['most_like_id']}`")
                    bits.append(chosen)
                if case.get('outcomes'):
                    bits.append(', '.join(case['outcomes'][:3]))
            elif 'transcript' in case:
                bits.append(f"{case['gene']} `{case['transcript']}`")
                if case.get('detail'):
                    bits.append(case['detail'])
                bits.append(f"{case['n_variants']} such transcripts")
            w(f"| {case['code']} | {case['shows']} | {where} | {' · '.join(bits) or '—'} |")
        w('')

    w('---\n')
    w('## Table S4 outcome coverage\n')
    w('| Outcome | Clusters in category_examples |')
    w('| :--- | ---: |')
    counts = {}
    if not cat.empty:
        for _, row in cat.iterrows():
            for outcome in row['outcomes']:
                counts[outcome] = counts.get(outcome, 0) + 1
    for outcome in S4_OUTCOMES:
        n = counts.get(outcome, 0)
        w(f"| {outcome} | {n if n else '**0 — see E1**'} |")
    w('')

    missing = [c for c in cases if c.get('missing')]
    w('---\n')
    w('## Gaps\n')
    if missing:
        for case in missing:
            w(f"- **{case['code']}** — {case['shows']}: no case in the current outputs.")
    else:
        w('- None: every case above resolved against the current outputs.')
    if not ladder.get('site_ptm'):
        w('- Site/PTM removal has no reviewable PDF here, though it is covered by unit')
        w('  tests. The entries exist in the database but no sampled cluster carries one.')
    w('- `no_gene_specified` cannot have a PDF: the figure is per-gene, and a cluster')
    w('  naming no gene has no transcripts, exons or domains to draw.')
    w('')
    return '\n'.join(out)


def main():
    if not os.path.exists(DEFAULT_DB_PATH):
        raise SystemExit(f'DoChaP database not found at {DEFAULT_DB_PATH}')
    con = sqlite3.connect(DEFAULT_DB_PATH)
    try:
        cases, clusters, cat, ladder = build_cases(con)
    finally:
        con.close()

    os.makedirs(CASES, exist_ok=True)
    with open(os.path.join(CASES, 'INDEX.md'), 'w') as handle:
        handle.write(render(cases, clusters, cat, ladder, DEFAULT_DB_PATH))
    with open(os.path.join(CASES, 'cases.json'), 'w') as handle:
        json.dump(cases, handle, indent=2)

    resolved = sum(1 for c in cases if not c.get('missing'))
    print(f'wrote review_cases/INDEX.md and cases.json — {resolved}/{len(cases)} cases resolved')
    for case in cases:
        if case.get('missing'):
            print(f"  no case for {case['code']} ({case['label']})")


if __name__ == '__main__':
    main()
