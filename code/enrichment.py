"""
Enrich a DOMAS results.csv with sequence-derived evidence, from the local
enrichment.sqlite (afdb_plddt / uniprot_feature / ensembl_sequence) + Pfam.

Adds, per (canonical_transcript, compared_transcript) comparison (the same
value repeated on that comparison's domain rows):

  canon_hmm            per-Pfam-family HMM coverage of the canonical protein
  alt_hmm              same for the compared protein (or SEQ-NOT-FOUND)
  changed_region       AA span that differs (prefix/suffix trimmed); canonical
                       coords, or ins:<alt coords> for an alt-only insertion
  region_plddt_mean    mean AlphaFold pLDDT over that span (canonical model)
  region_plddt_pct70   % of that span with pLDDT > 70
  region_uniprot_sites UniProt features overlapping that span

These are the cross-transcript-stable signals the InterPro type filter can't
give (HMM family is uniform across transcripts; pLDDT/features are per residue).
Everything is read from local disk - no network.
"""
import argparse
import logging
import os
import subprocess
import sqlite3
import tempfile

import pandas as pd

import build  # _parse_hmmsearch_domtbl (bitscores)
from junction_analisys import NON_COMPARISON_EVENTS
from utils import (calc_spade_score, hmm_change_impact, insertion_impact,
                   impact_probability, fold_change_probability, fold_change_call,
                   fold_change_note, impact_note, missing_features,
                   NO_ACCESSION_NOTE)

_PRIMARY_FT_SKIP = {'CHAIN', 'SIGNAL', 'Chain', 'Signal'}
_IMPACT_RANK = {'none': 0, 'gain': 1, 'low': 1, 'moderate': 2, 'high': 3}
# UniProt functional-residue feature keys (losing one of these matters regardless
# of how much domain coverage/structure changed).
_FUNCTIONAL_FT = ('ACT_SITE', 'BINDING', 'SITE', 'MOTIF', 'MOD_RES', 'DISULFID',
                  'CROSSLNK', 'METAL', 'CA_BIND', 'DNA_BIND', 'NP_BIND', 'LIPID', 'CARBOHYD')


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

    def mean_plddt(self, uniprot):
        """Whole-structure mean pLDDT of the canonical AlphaFold model, or None.

        Used as the fold_change_prob reliability guard: DOMAS never folds the
        alternative isoform, so it cannot check that model's quality directly - but
        a protein AlphaFold models poorly overall is one whose isoform it also models
        poorly, and canonical mean pLDDT identifies that regime at AUC 0.873 on the
        benchmark. See utils.fold_change_note.
        """
        if not uniprot:
            return None
        row = self.enr.execute(
            "SELECT mean_plddt FROM afdb_plddt WHERE accession=?", (uniprot,)).fetchone()
        return row[0] if row and row[0] is not None else None

    def rsa(self, uniprot):
        """Per-residue relative solvent accessibility (0=buried .. 1=exposed) for
        a UniProt accession, from the canonical AlphaFold model. None if missing.
        Optional table (present only if afdb_rsa was provisioned)."""
        if not uniprot:
            return None
        try:
            row = self.enr.execute(
                "SELECT rsa FROM afdb_rsa WHERE accession=?", (uniprot,)).fetchone()
        except sqlite3.OperationalError:
            return None
        return [float(x) for x in row[0].split(',')] if row else None

    def pae_global(self, uniprot):
        """Whole-structure mean PAE (predicted aligned error) of the canonical
        AlphaFold model for a UniProt accession, or None. Higher = floppier /
        multi-domain (relative positions uncertain); the dominant fold-change
        feature. Optional table (present only if afdb_pae was provisioned)."""
        if not uniprot:
            return None
        try:
            row = self.enr.execute(
                "SELECT pae_global FROM afdb_pae WHERE accession=?", (uniprot,)).fetchone()
        except sqlite3.OperationalError:
            return None
        return row[0] if row and row[0] is not None else None

    def loeuf(self, uniprot):
        """gnomAD LOEUF (gene loss-intolerance; lower = more constrained) for a
        UniProt accession, or None. Optional table (present only if gnomAD was
        provisioned)."""
        if not uniprot:
            return None
        try:
            row = self.enr.execute(
                "SELECT loeuf FROM gene_constraint WHERE accession=?", (uniprot,)).fetchone()
        except sqlite3.OperationalError:
            return None
        return row[0] if row and row[0] is not None else None

    def am_pathogenicity(self, uniprot):
        """Per-residue mean AlphaMissense pathogenicity for a UniProt accession,
        or None. The am_pathogenicity table is optional (present only if
        AlphaMissense was provisioned)."""
        if not uniprot:
            return None
        try:
            row = self.enr.execute(
                "SELECT am FROM am_pathogenicity WHERE accession=?", (uniprot,)).fetchone()
        except sqlite3.OperationalError:
            return None
        if not row:
            return None
        return [float(x) if x else None for x in row[0].split(',')]

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

    def fold_state(self, uniprot, lo, hi):
        """Curated UniProt structure state of region [lo,hi]: 'folded' if it
        overlaps an experimental HELIX/STRAND, else 'disordered' if it overlaps a
        REGION annotated 'Disordered', else None. Experimental secondary
        structure wins over the (MobiDB-derived) disorder call."""
        if not uniprot:
            return None
        q = ("SELECT ftype, note FROM uniprot_feature WHERE acc=? "
             "AND start IS NOT NULL AND NOT (end<? OR start>?)")
        folded = disordered = False
        for ftype, note in self.enr.execute(q, (uniprot, lo, hi)):
            if ftype in ('HELIX', 'STRAND'):
                folded = True
            elif ftype == 'REGION' and note and 'isorder' in note.lower():
                disordered = True
        return 'folded' if folded else ('disordered' if disordered else None)

    def functional_sites(self, uniprot, lo, hi):
        """UniProt functional-residue features overlapping [lo,hi]."""
        if not uniprot:
            return []
        ph = ",".join("?" * len(_FUNCTIONAL_FT))
        q = (f"SELECT ftype,start,end,note FROM uniprot_feature WHERE acc=? "
             f"AND ftype IN ({ph}) AND start IS NOT NULL AND NOT (end<? OR start>?) ORDER BY start")
        out = []
        for ftype, s, e, note in self.enr.execute(q, (uniprot, *_FUNCTIONAL_FT, lo, hi)):
            out.append(f"{ftype} {s}-{e}{(' ' + note[:22]) if note else ''}".strip())
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
    def changed_region(canonical_seq, alt_seq):
        """Divergent span between the canonical and alternative protein sequences,
        found by trimming the shared prefix and shared suffix.

        Returns the change described in BOTH coordinate systems, because an
        insertion has residues only in the alt sequence (no canonical
        coordinate) and a deletion only in the canonical sequence:

            {'canon': (lo, hi) | None,   # 1-based span on the canonical sequence
             'alt':   (lo, hi) | None,   # 1-based span on the alternative sequence
             'kind':  'substitution' | 'insertion' | 'deletion'}

        None if either sequence is missing or the two are identical.
        """
        if not canonical_seq or not alt_seq:
            return None
        shorter_len = min(len(canonical_seq), len(alt_seq))
        shared_prefix = 0
        while (shared_prefix < shorter_len
               and canonical_seq[shared_prefix] == alt_seq[shared_prefix]):
            shared_prefix += 1
        shared_suffix = 0
        while (shared_suffix < shorter_len - shared_prefix
               and canonical_seq[-1 - shared_suffix] == alt_seq[-1 - shared_suffix]):
            shared_suffix += 1
        canon_lo, canon_hi = shared_prefix + 1, len(canonical_seq) - shared_suffix
        alt_lo, alt_hi = shared_prefix + 1, len(alt_seq) - shared_suffix
        canon_span = (canon_lo, canon_hi) if canon_hi >= canon_lo else None
        alt_span = (alt_lo, alt_hi) if alt_hi >= alt_lo else None
        if canon_span is None and alt_span is None:
            return None
        kind = ('insertion' if canon_span is None else
                'deletion' if alt_span is None else
                'substitution')
        return {'canon': canon_span, 'alt': alt_span, 'kind': kind}


def enrich_results(results_csv, out_csv, enrichment_db, pfam_hmm, dochap_con):
    df = pd.read_csv(results_csv)
    e = Enricher(enrichment_db, pfam_hmm, dochap_con)

    # Only enrich transcripts that were actually compared to the canonical.
    # Rows carrying a NON_COMPARISON event (no_canonical_features /
    # no_unique_features / transcript_doesnt_have_features / feature_not_mapped
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
        if reg and reg['canon']:
            # substitution / deletion: pLDDT and UniProt features are indexed on
            # the canonical sequence, so they only apply to a canonical span.
            lo, hi = reg['canon']
            region = f'{lo}-{hi}'
            if reg['alt']:
                # flag when the alt span adds net residues (expansion/insertion)
                net_added = (reg['alt'][1] - reg['alt'][0] + 1) - (hi - lo + 1)
                if net_added > 0:
                    region += f' +{net_added}alt'
            pl = e.plddt(uni)
            if pl and len(pl) == len(a):
                seg = pl[lo - 1:hi]
                if seg:
                    plddt_mean = round(sum(seg) / len(seg))
                    plddt_pct = round(100 * sum(x > 70 for x in seg) / len(seg))
            sites = ' | '.join(e.features(uni, lo, hi))
        elif reg and reg['kind'] == 'insertion':
            # alt-only residues: report the span in alt coordinates.
            alo, ahi = reg['alt']
            region = f'ins:{alo}-{ahi}'
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
      impact_prob / fold_change_prob   : calibrated functional / structural
                                         probabilities (utils, see their docstrings)
      fold_change_call                 : change / preserved / uncertain triage
      fold_change_note                 : why a structural call was withheld (canonical
                                         AlphaFold model too poor to trust); blank when
                                         the structure is well resolved
      functional_sites / region_am_mean: the UniProt and AlphaMissense evidence the
                                         above are computed FROM

    Both are computed from a single hmmsearch over the analysed transcripts plus
    canonical AlphaFold pLDDT; non-comparable rows are left blank. Reads pfam_match
    if you later prefer the cache, but here scans live for the transcripts in play.

    THESE COLUMNS ARE NOT INDEPENDENT EVIDENCE. `impact`, `impact_prob`,
    `region_am_mean` and `functional_sites` are largely one measurement shown four
    ways: impact_prob is ~78% rank-correlated with region_am (and beats raw
    AlphaMissense by only +0.013 AUC), while hmm_change_impact weights on that SAME
    AlphaMissense value using the same 0.564/0.34 thresholds and on the same functional
    sites. A row reading impact=high AND impact_prob=0.85 AND region_am_mean=0.72 with
    functional_sites populated is ONE piece of evidence displayed four times, not four
    agreeing signals - don't read the agreement as corroboration. The genuinely
    independent split is functional (impact_prob) vs structural (fold_change_prob):
    those two correlate at rho -0.05 across the benchmark.
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

    # impact per (canonical, compared): worst severity across changed Pfam
    # families, structure-weighted by UniProt fold-state (else pLDDT) and raised
    # when the changed region hits a UniProt functional residue.
    impact, func_col, am_col, prob_col, fold_col, note_col = {}, {}, {}, {}, {}, {}
    inote_col, plddt_col = {}, {}
    for canon, comp in {(r['canonical_transcript_id'], r['transcript_id'])
                        for _, r in analyzable[['canonical_transcript_id', 'transcript_id']].dropna().iterrows()}:
        cm = {m[0]: m for m in matches.get(canon, [])}
        am = {m[0]: m for m in matches.get(comp, [])}
        uni = e.uniprot_for(canon)
        pl = e.plddt(uni)
        # region-level signals over the changed span: UniProt fold-state +
        # functional sites, and mean AlphaMissense constraint (if provisioned).
        reg = Enricher.changed_region(seqs.get(canon), seqs.get(comp))
        fold_state, sites, region_am, region_rsa, buried_frac, identity = None, [], None, None, None, None
        added_label = 'none'
        if reg and reg['canon']:
            lo, hi = reg['canon']
            cs, alt_seq = seqs.get(canon), seqs.get(comp)
            if cs and alt_seq:                   # sequence identity (%) from the trimmed change
                shared = (lo - 1) + (len(cs) - hi)
                identity = 100.0 * shared / max(len(cs), len(alt_seq))
            fold_state = e.fold_state(uni, lo, hi)
            sites = e.functional_sites(uni, lo, hi)
            am_arr = e.am_pathogenicity(uni)
            if am_arr:
                seg = [v for v in am_arr[lo - 1:hi] if v is not None]
                region_am = round(sum(seg) / len(seg), 3) if seg else None
            rsa_arr = e.rsa(uni)                 # AlphaFold burial of the changed region
            if rsa_arr and len(rsa_arr) >= hi:
                seg_r = rsa_arr[lo - 1:hi]
                if seg_r:
                    region_rsa = round(sum(seg_r) / len(seg_r), 3)
                    buried_frac = round(sum(x < 0.25 for x in seg_r) / len(seg_r), 3)
        if reg and reg['alt']:
            # residues the alt isoform ADDS - a pure insertion, or the surplus of
            # a substitution whose alt span is longer than the canonical span -
            # have no canonical coordinate, so the loss-based per-Pfam loop below
            # misses them. Score the net-added segment by size + whether it lands
            # inside a Pfam domain on the alt sequence (disruptive) vs a
            # terminus/linker.
            alo, ahi = reg['alt']
            alt_span_len = ahi - alo + 1
            canon_span_len = (hi - lo + 1) if reg['canon'] else 0
            net_added = alt_span_len - canon_span_len
            if net_added > 0:
                inside_domain = any(
                    min(ahi, m[4]) - max(alo, m[3]) + 1 >= 0.5 * alt_span_len
                    for m in matches.get(comp, []))
                if reg['canon']:
                    f_lo, f_hi = lo, hi           # canonical residues being replaced/expanded
                else:
                    canon_len = len(seqs.get(canon) or '')
                    flank = alo - 1               # canonical residue before a pure insertion
                    f_lo, f_hi = max(1, flank), min(canon_len, flank + 1)
                junction_fold = (e.fold_state(uni, f_lo, f_hi) if f_hi >= f_lo else None)
                added_label = insertion_impact(net_added, inside_domain, junction_fold)
        hits = bool(sites)
        func_col[(canon, comp)] = ' | '.join(sites)
        am_col[(canon, comp)] = region_am if region_am is not None else ''
        # a functional residue in the changed region is at least 'moderate' even
        # if no Pfam coverage changes (the disordered-motif case); net-added
        # residues seed their own size/placement level (the loss-based loop can't).
        best_rank, best_label = (_IMPACT_RANK['moderate'], 'moderate') if hits else (0, 'none')
        if _IMPACT_RANK.get(added_label, 0) > best_rank:
            best_rank, best_label = _IMPACT_RANK[added_label], added_label
        max_cov_loss = 0.0
        for acc in set(cm) | set(am):
            ccov = cm[acc][2] if acc in cm else None
            acov = am[acc][2] if acc in am else None
            if ccov is not None:
                max_cov_loss = max(max_cov_loss, ccov - (acov or 0))
            cpl = None
            if acc in cm and pl:
                s, en = cm[acc][3], cm[acc][4]
                seg = pl[s - 1:en] if en <= len(pl) else []
                cpl = round(sum(seg) / len(seg), 1) if seg else None
            lab = hmm_change_impact(ccov, acov, cpl, fold_state=fold_state,
                                    hits_functional_site=hits, region_am=region_am,
                                    region_rsa=region_rsa)
            if _IMPACT_RANK.get(lab, 0) > best_rank:
                best_rank, best_label = _IMPACT_RANK[lab], lab
        impact[(canon, comp)] = best_label
        # two calibrated continuous companion scores, defined only where there is a
        # canonical changed region (not pure insertions): impact_prob = functional
        # (pathogenicity-relevance), fold_change_prob = structural (P(TM<0.5)).
        loeuf = e.loeuf(uni)
        if reg and reg['canon']:
            # Each score is withheld when ITS dominant feature is missing rather than
            # emitted off a median-imputed constant (see utils._DOMINANT); the note
            # records that, plus any secondary feature that was imputed.
            imiss = missing_features(region_am=region_am, loeuf=loeuf,
                                     max_cov_loss=max_cov_loss, buried_frac=buried_frac)
            # distinguish "this protein has no data" from "this transcript has no
            # accession to look data up by" - the second is a DoChaP mapping gap and
            # accounts for ~24% of otherwise-comparable human genes, so blaming the
            # enrichment tables for it would point at the wrong fix.
            inote_col[(canon, comp)] = NO_ACCESSION_NOTE if uni is None else impact_note(imiss)
            prob_col[(canon, comp)] = '' if region_am is None else impact_probability(
                region_am=region_am, loeuf=loeuf, max_cov_loss=max_cov_loss, buried_frac=buried_frac)
            # fold_change_prob is PAE-driven (E38): pae_global (whole canonical
            # structure) dominates, with identity, max_cov_loss and protein length.
            pae_global = e.pae_global(uni)
            protL = len(seqs.get(canon) or '') or None
            fmiss = missing_features(pae_global=pae_global, identity=identity,
                                     max_cov_loss=max_cov_loss, protL=protL)
            fold_col[(canon, comp)] = '' if pae_global is None else fold_change_probability(
                pae_global=pae_global, identity=identity,
                max_cov_loss=max_cov_loss, protL=protL)
            # reliability guard: fold_change_prob is calibrated on proteins AlphaFold
            # resolves well, and DOMAS cannot check the ISOFORM model's quality (it
            # never folds one). Canonical mean pLDDT is the proxy - below it the call
            # is downgraded to 'uncertain' and the note carries the reason.
            plddt_col[(canon, comp)] = e.mean_plddt(uni)
            note_col[(canon, comp)] = (NO_ACCESSION_NOTE if uni is None
                                       else fold_change_note(plddt_col[(canon, comp)], missing=fmiss))
        else:
            prob_col[(canon, comp)] = fold_col[(canon, comp)] = ''
            note_col[(canon, comp)] = 'not scored: pure insertion (no canonical changed region)'
            inote_col[(canon, comp)] = note_col[(canon, comp)]

    def _blank_nc(row, val):
        return '' if row['event_type'] in NON_COMPARISON_EVENTS else val

    df['spade_canonical'] = df.apply(lambda r: _blank_nc(r, spade.get(r['canonical_transcript_id'], '')), axis=1)
    df['spade_compared'] = df.apply(lambda r: _blank_nc(r, spade.get(r['transcript_id'], '')), axis=1)
    df['impact'] = df.apply(lambda r: _blank_nc(r, impact.get((r['canonical_transcript_id'], r['transcript_id']), '')), axis=1)
    df['impact_prob'] = df.apply(lambda r: _blank_nc(r, prob_col.get((r['canonical_transcript_id'], r['transcript_id']), '')), axis=1)
    df['fold_change_prob'] = df.apply(lambda r: _blank_nc(r, fold_col.get((r['canonical_transcript_id'], r['transcript_id']), '')), axis=1)
    df['fold_change_note'] = df.apply(lambda r: _blank_nc(r, note_col.get((r['canonical_transcript_id'], r['transcript_id']), '')), axis=1)
    df['impact_prob_note'] = df.apply(lambda r: _blank_nc(r, inote_col.get((r['canonical_transcript_id'], r['transcript_id']), '')), axis=1)
    # canonical mean pLDDT drives the one-sided guard inside fold_change_call: a
    # 'preserved' verdict on a poorly-resolved canonical becomes 'uncertain'
    _cpl = df.apply(lambda r: _blank_nc(r, plddt_col.get((r['canonical_transcript_id'], r['transcript_id']), None)), axis=1)
    df['fold_change_call'] = [fold_change_call(p, canonical_plddt=(v if v != '' else None)) or ''
                              for p, v in zip(df['fold_change_prob'], _cpl)]
    df['functional_sites'] = df.apply(lambda r: _blank_nc(r, func_col.get((r['canonical_transcript_id'], r['transcript_id']), '')), axis=1)
    df['region_am_mean'] = df.apply(lambda r: _blank_nc(r, am_col.get((r['canonical_transcript_id'], r['transcript_id']), '')), axis=1)
    df.to_csv(out_csv, index=False)
    return df


# ---------------------------------------------------------------------------
# Provisioning CLI.
#
# Building the enrichment database has nothing to do with analysing junctions,
# but it used to live in domas.py as a "setup mode" that short-circuited the
# analysis CLI: eight arguments, a branch that returned from argument validation
# before any of it ran, and a dispatch in main() that exited early. It belongs
# with the enrichment code it provisions.
#
#   python3 enrichment.py -download -build -data_dir enrichment_data
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Provision the local enrichment database (afdb_plddt / "
                    "uniprot_feature / uniprot_alias / ensembl_sequence + Pfam), "
                    "used to enrich a DOMAS results.csv with sequence-derived evidence.")
    parser.add_argument("-download", action="store_true",
                        help="Download the raw enrichment sources into -data_dir.")
    parser.add_argument("-build", action="store_true",
                        help="Build -enrichment_db from the raw files in -data_dir "
                             "(may be combined with -download to do both).")
    parser.add_argument("-data_dir", type=str, default="enrichment_data",
                        help="Directory for the raw source files (default: ./enrichment_data).")
    parser.add_argument("-enrichment_db", type=str, default=None,
                        help="Output sqlite path for -build (default: <data_dir>/enrichment.sqlite).")
    parser.add_argument("-species", type=str, default=None,
                        help="Comma-separated species subset for -download/-build "
                             "(default: all; choices: H_sapiens,M_musculus,R_norvegicus).")
    parser.add_argument("-only", type=str, default=None,
                        help="Comma-separated source subset for -download/-build "
                             "(default: all; choices: uniprot,afdb,ensembl,pfam).")
    parser.add_argument("-delete_raw", action="store_true",
                        help="After -build, delete the large raw downloads (keeps Pfam-A.hmm).")
    parser.add_argument("-force_download", action="store_true",
                        help="Re-download sources even if already present in -data_dir.")

    parser.add_argument("-enrich", type=str, default=None, metavar="RESULTS_CSV",
                        help="Enrich a DOMAS results.csv with sequence-derived evidence "
                             "(HMM coverage, changed region, AlphaFold pLDDT, UniProt "
                             "features) and write -out. Requires -dochap and a built "
                             "-enrichment_db.")
    parser.add_argument("-scores", action="store_true",
                        help="With -enrich, also add the analysis-level score columns "
                             "(SPADE-style domain integrity, change impact, fold-change "
                             "probability). Run as a second pass over the enriched table.")
    parser.add_argument("-out", type=str, default=None,
                        help="Output CSV for -enrich (default: the input with a "
                             "'.enriched.csv' suffix).")
    parser.add_argument("-dochap", type=str, default=None,
                        help="Path to the DoChaP sqlite db (required with -enrich).")
    parser.add_argument("-pfam_hmm", type=str, default=None,
                        help="Path to Pfam-A.hmm (default: <data_dir>/pfam/Pfam-A.hmm, "
                             "where -build leaves it).")

    args = parser.parse_args(argv)
    if not args.download and not args.build and not args.enrich:
        parser.error("nothing to do: give -download, -build and/or -enrich")
    if args.enrichment_db is None:
        args.enrichment_db = os.path.join(args.data_dir, "enrichment.sqlite")
    if args.pfam_hmm is None:
        args.pfam_hmm = os.path.join(args.data_dir, "pfam", "Pfam-A.hmm")

    if args.enrich:
        if not args.dochap:
            parser.error("-enrich requires -dochap")
        for label, path in (("-enrich input", args.enrich), ("-dochap", args.dochap),
                            ("-enrichment_db", args.enrichment_db), ("-pfam_hmm", args.pfam_hmm)):
            # Checked up front rather than failing partway through an hmmsearch that
            # can take minutes before it touches the missing file.
            if not os.path.exists(path):
                parser.error(f"{label} not found at {path}"
                             + ("  (run -download -build first)"
                                if label in ("-enrichment_db", "-pfam_hmm") else ""))
        if args.out is None:
            base = args.enrich[:-4] if args.enrich.endswith(".csv") else args.enrich
            args.out = f"{base}.enriched.csv"
    elif args.scores:
        parser.error("-scores only applies together with -enrich")
    return args


def run_setup(args):
    """Provision the enrichment database: -download fetches raw sources,
    -build parses them into -enrichment_db. Both may be given together."""
    import download

    species = [s.strip() for s in args.species.split(",")] if args.species else None
    only = [s.strip() for s in args.only.split(",")] if args.only else None
    if args.download:
        download.download_all(args.data_dir, species=species, only=only, force=args.force_download)
    if args.build:
        build.build_all(args.data_dir, args.enrichment_db, species=species, only=only,
                        delete_raw=args.delete_raw)


def run_enrich(args):
    """Enrich a results.csv, optionally following it with the score pass.

    The two passes are separate functions over the same table: enrich_results()
    adds the sequence-derived evidence columns, add_scores() the analysis-level
    scores derived from them. Chained here so -scores reads the enriched table
    rather than the raw one.
    """
    dochap_con = sqlite3.connect(args.dochap)
    try:
        logging.info("enriching %s -> %s", args.enrich, args.out)
        enrich_results(args.enrich, args.out, args.enrichment_db, args.pfam_hmm, dochap_con)
        if args.scores:
            logging.info("adding score columns to %s", args.out)
            add_scores(args.out, args.out, args.enrichment_db, args.pfam_hmm, dochap_con)
    finally:
        dochap_con.close()
    logging.info("wrote %s", args.out)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = parse_args()
    if args.download or args.build:
        run_setup(args)
    if args.enrich:
        run_enrich(args)


if __name__ == "__main__":
    main()
