"""
Enrich a DOMAS results.csv with sequence-derived evidence, from the local
enrichment.sqlite (afdb_plddt / uniprot_feature / ensembl_sequence) + Pfam.

Adds, per (canonical_transcript, compared_transcript) comparison (the same
value repeated on that comparison's domain rows):

  canon_hmm            per-Pfam-family HMM coverage of the canonical protein
  alt_hmm              same for the compared protein (or SEQ-NOT-FOUND)
  changed_region       canonical AA span that differs (prefix/suffix trimmed)
  region_plddt_mean    mean AlphaFold pLDDT over that span (canonical model)
  region_plddt_pct70   % of that span with pLDDT > 70
  region_uniprot_sites UniProt features overlapping that span

These are the cross-transcript-stable signals the InterPro type filter can't
give (HMM family is uniform across transcripts; pLDDT/features are per residue).
Everything is read from local disk - no network.
"""
import os
import subprocess
import sqlite3
import tempfile

import pandas as pd

import build  # _parse_hmmsearch_domtbl (bitscores)
from junction_analisys import NON_COMPARISON_EVENTS
from utils import calc_spade_score, hmm_change_impact

_PRIMARY_FT_SKIP = {'CHAIN', 'SIGNAL', 'Chain', 'Signal'}
_IMPACT_RANK = {'none': 0, 'gain': 1, 'low': 1, 'moderate': 2, 'high': 3}


class Enricher:
    def __init__(self, enrichment_db, pfam_hmm, dochap_con):
        self.enr = sqlite3.connect(enrichment_db)
        self.pfam = pfam_hmm
        self.dochap = dochap_con

    # --- lookups -----------------------------------------------------------
    def ensp_for(self, transcript_ensembl_id):
        """transcript id (versioned) -> bare ENSP, or None (e.g. RefSeq NM_/XM_)."""
        row = self.dochap.execute(
            "SELECT protein_ensembl_id FROM Transcripts WHERE transcript_ensembl_id=?",
            (transcript_ensembl_id,)).fetchone()
        if not row or not row[0]:
            return None
        return row[0].split('.')[0]

    def uniprot_for(self, transcript_ensembl_id):
        row = self.dochap.execute(
            "SELECT p.protein_interpro_id FROM Transcripts t JOIN Proteins p "
            "ON t.protein_ensembl_id=p.protein_ensembl_id "
            "WHERE t.transcript_ensembl_id=?", (transcript_ensembl_id,)).fetchone()
        return row[0] if row and row[0] else None

    def seq(self, ensp):
        if not ensp:
            return None
        row = self.enr.execute(
            "SELECT seq FROM ensembl_sequence WHERE protein_ensembl_id=?", (ensp,)).fetchone()
        return row[0] if row else None

    def plddt(self, uniprot):
        if not uniprot:
            return None
        row = self.enr.execute(
            "SELECT plddt FROM afdb_plddt WHERE accession=?", (uniprot,)).fetchone()
        return [float(x) for x in row[0].split(',')] if row else None

    def features(self, uniprot, lo, hi):
        if not uniprot:
            return []
        q = ("SELECT ftype,start,end,note FROM uniprot_feature WHERE acc=? "
             "AND start IS NOT NULL AND NOT (end<? OR start>?) ORDER BY start")
        out = []
        for ftype, s, e, note in self.enr.execute(q, (uniprot, lo, hi)):
            if ftype in _PRIMARY_FT_SKIP:
                continue
            out.append(f"{ftype} {s}-{e}{(' ' + note) if note else ''}".strip())
        return out

    # --- HMM (batched: one hmmscan over all unique sequences) --------------
    def hmm_scan_many(self, id_to_seq):
        """{id: seq} -> {id: [(family, hmm_cov_pct), ...]}."""
        ids = [i for i, s in id_to_seq.items() if s]
        if not ids:
            return {}
        with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as fa:
            for i in ids:
                fa.write(f">{i}\n{id_to_seq[i]}\n")
            fa_path = fa.name
        out = fa_path + '.domtbl'
        subprocess.run(['hmmscan', '--cut_ga', '--domtblout', out, self.pfam, fa_path],
                       capture_output=True, check=True)
        hits = {i: [] for i in ids}
        with open(out) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                f = line.split()
                if len(f) < 23:
                    continue
                fam, hmm_len, qid = f[0], int(f[2]), f[3]
                hmm_from, hmm_to = int(f[15]), int(f[16])
                cov = round(100 * (hmm_to - hmm_from + 1) / hmm_len)
                if qid in hits:
                    hits[qid].append((fam, cov))
        os.remove(fa_path)
        os.remove(out)
        return hits

    def hmm_hits_many(self, id_to_seq):
        """{id: seq} -> {id: [(family, ali_start, ali_end, hmm_cov_pct), ...]}.

        Like hmm_scan_many but keeps the alignment coordinates on the protein,
        for drawing the HMM track and building the changed-element table."""
        ids = [i for i, s in id_to_seq.items() if s]
        if not ids:
            return {}
        with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as fa:
            for i in ids:
                fa.write(f">{i}\n{id_to_seq[i]}\n")
            fa_path = fa.name
        out = fa_path + '.domtbl'
        subprocess.run(['hmmscan', '--cut_ga', '--domtblout', out, self.pfam, fa_path],
                       capture_output=True, check=True)
        hits = {i: [] for i in ids}
        with open(out) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                f = line.split()
                if len(f) < 23:
                    continue
                fam, hmm_len, qid = f[0], int(f[2]), f[3]
                hmm_from, hmm_to = int(f[15]), int(f[16])
                ali_from, ali_to = int(f[17]), int(f[18])
                cov = round(100 * (hmm_to - hmm_from + 1) / hmm_len)
                if qid in hits:
                    hits[qid].append((fam, ali_from, ali_to, cov))
        os.remove(fa_path)
        os.remove(out)
        return hits

    @staticmethod
    def _fmt_hmm(hitlist):
        if hitlist is None:
            return 'SEQ-NOT-FOUND'
        if not hitlist:
            return ''
        return '; '.join(f'{fam}:{cov}' for fam, cov in hitlist)

    @staticmethod
    def changed_region(a, b):
        """Canonical AA span [lo,hi] that differs from the compared seq, by
        trimming the shared prefix and suffix. None if identical or missing."""
        if not a or not b:
            return None
        p = 0
        while p < min(len(a), len(b)) and a[p] == b[p]:
            p += 1
        s = 0
        while s < min(len(a), len(b)) - p and a[-1 - s] == b[-1 - s]:
            s += 1
        lo, hi = p + 1, len(a) - s
        return (lo, hi) if hi >= lo else None


def enrich_results(results_csv, out_csv, enrichment_db, pfam_hmm, dochap_con):
    df = pd.read_csv(results_csv)
    e = Enricher(enrichment_db, pfam_hmm, dochap_con)

    # Only enrich transcripts that were actually compared to the canonical.
    # Rows carrying a NON_COMPARISON event (no_canonical_junctions /
    # no_unique_junctions / transcript_doesnt_have_junctions / junction_not_mapped
    # / gene_not_in_db / ...) never produced a domain comparison, so there is
    # nothing to HMM-scan for them - skip so we don't scan unanalyzable proteins.
    analyzable = df[~df['event_type'].isin(NON_COMPARISON_EVENTS)]
    pairs = analyzable[['canonical_transcript_id', 'transcript_id']].dropna().drop_duplicates()
    pairs = [tuple(r) for r in pairs.to_numpy()]

    # gather unique sequences to scan once
    seqs = {}
    for canon, comp in pairs:
        for tx in (canon, comp):
            if tx not in seqs:
                seqs[tx] = e.seq(e.ensp_for(tx))
    hits = e.hmm_scan_many({tx: s for tx, s in seqs.items() if s})

    enr_cols = {}
    for canon, comp in pairs:
        a, b = seqs.get(canon), seqs.get(comp)
        uni = e.uniprot_for(canon)
        reg = e.changed_region(a, b)
        plddt_mean = plddt_pct = region = sites = ''
        if reg:
            lo, hi = reg
            region = f'{lo}-{hi}'
            pl = e.plddt(uni)
            if pl and len(pl) == len(a):
                seg = pl[lo - 1:hi]
                if seg:
                    plddt_mean = round(sum(seg) / len(seg))
                    plddt_pct = round(100 * sum(x > 70 for x in seg) / len(seg))
            sites = ' | '.join(e.features(uni, lo, hi))
        enr_cols[(canon, comp)] = {
            'canon_hmm': e._fmt_hmm(hits.get(canon, [] if a else None)),
            'alt_hmm': e._fmt_hmm(hits.get(comp, [] if b else None)),
            'changed_region': region,
            'region_plddt_mean': plddt_mean,
            'region_plddt_pct70': plddt_pct,
            'region_uniprot_sites': sites,
        }

    def row_enr(row, key):
        # non-comparable rows never carry enrichment, even when the same
        # transcript pair was scanned for a real comparison in another cluster.
        if row['event_type'] in NON_COMPARISON_EVENTS:
            return ''
        k = (row['canonical_transcript_id'], row['transcript_id'])
        return enr_cols.get(k, {}).get(key, '')

    for col in ('canon_hmm', 'alt_hmm', 'changed_region', 'region_plddt_mean',
                'region_plddt_pct70', 'region_uniprot_sites'):
        df[col] = df.apply(lambda r: row_enr(r, col), axis=1)

    df.to_csv(out_csv, index=False)
    return df


def _hmmsearch_matches(id_to_seq, pfam_hmm):
    """{id: seq} -> {id: [(pfam_acc, bitscore, coverage, ali_start, ali_end), ...]}
    via one hmmsearch pass (profiles -> sequences)."""
    ids = {i: s for i, s in id_to_seq.items() if s}
    if not ids:
        return {}
    with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as fa:
        for i, s in ids.items():
            fa.write(f">{i}\n{s}\n")
        fa_path = fa.name
    out = fa_path + '.domtbl'
    subprocess.run(['hmmsearch', '--cut_ga', '--domtblout', out, pfam_hmm, fa_path],
                   check=True, capture_output=True)
    res = {}
    for r in build._parse_hmmsearch_domtbl(out):
        # r = (prot, acc, name, ali_s, ali_e, env_s, env_e, hmm_f, hmm_t, hmm_len, bits, ie, cov)
        res.setdefault(r[0], []).append((r[1], r[10], r[12], r[3], r[4]))
    os.remove(fa_path)
    os.remove(out)
    return res


def add_scores(results_csv, out_csv, enrichment_db, pfam_hmm, dochap_con):
    """Add analysis-level score columns to a results table so the GUI reads them
    directly (no recompute in the PDF/GUI):

      spade_canonical / spade_compared : SPADE-style per-isoform domain-integrity
                                         score (utils.calc_spade_score, bitscore)
      impact                           : worst per-comparison change severity
                                         (utils.hmm_change_impact)

    Both are computed from a single hmmsearch over the analysed transcripts plus
    canonical AlphaFold pLDDT; non-comparable rows are left blank. Reads pfam_match
    if you later prefer the cache, but here scans live for the transcripts in play.
    """
    df = pd.read_csv(results_csv)
    e = Enricher(enrichment_db, pfam_hmm, dochap_con)
    analyzable = df[~df['event_type'].isin(NON_COMPARISON_EVENTS)]

    txs = set(analyzable['canonical_transcript_id'].dropna()) | set(analyzable['transcript_id'].dropna())
    seqs = {tx: e.seq(e.ensp_for(tx)) for tx in txs}
    matches = _hmmsearch_matches(seqs, pfam_hmm)      # {tx: [(acc,bits,cov,ali_s,ali_e)]}

    # SPADE per gene (rows of a gene share canonical_transcript_id)
    spade = {}
    for canon, grp in analyzable.groupby('canonical_transcript_id'):
        gene_txs = {canon} | set(grp['transcript_id'].dropna())
        doms = {tx: [(m[0], m[1]) for m in matches.get(tx, [])] for tx in gene_txs}
        for tx, info in calc_spade_score(doms).items():
            spade[tx] = info['spade_score']

    # impact per (canonical, compared): worst severity across changed Pfam families
    impact = {}
    for canon, comp in {(r['canonical_transcript_id'], r['transcript_id'])
                        for _, r in analyzable[['canonical_transcript_id', 'transcript_id']].dropna().iterrows()}:
        cm = {m[0]: m for m in matches.get(canon, [])}
        am = {m[0]: m for m in matches.get(comp, [])}
        pl = e.plddt(e.uniprot_for(canon))
        best_rank, best_label = 0, 'none'
        for acc in set(cm) | set(am):
            ccov = cm[acc][2] if acc in cm else None
            acov = am[acc][2] if acc in am else None
            cpl = None
            if acc in cm and pl:
                s, en = cm[acc][3], cm[acc][4]
                seg = pl[s - 1:en] if en <= len(pl) else []
                cpl = round(sum(seg) / len(seg), 1) if seg else None
            lab = hmm_change_impact(ccov, acov, cpl)
            if _IMPACT_RANK.get(lab, 0) > best_rank:
                best_rank, best_label = _IMPACT_RANK[lab], lab
        impact[(canon, comp)] = best_label

    def _spade(row, col):
        return '' if row['event_type'] in NON_COMPARISON_EVENTS else spade.get(row[col], '')

    df['spade_canonical'] = df.apply(lambda r: _spade(r, 'canonical_transcript_id'), axis=1)
    df['spade_compared'] = df.apply(lambda r: _spade(r, 'transcript_id'), axis=1)
    df['impact'] = df.apply(
        lambda r: '' if r['event_type'] in NON_COMPARISON_EVENTS
        else impact.get((r['canonical_transcript_id'], r['transcript_id']), ''), axis=1)
    df.to_csv(out_csv, index=False)
    return df
