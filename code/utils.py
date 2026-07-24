import math
import os

import pandas as pd


# rMATS-turbo per-event result files (junction-count variant). The JCEC files
# carry identical coordinates, so JC is enough for DOMAS's coordinate mapping.
_RMATS_EVENT_FILES = {
    'SE': 'SE.MATS.JC.txt',
    'A5SS': 'A5SS.MATS.JC.txt',
    'A3SS': 'A3SS.MATS.JC.txt',
    'MXE': 'MXE.MATS.JC.txt',
    'RI': 'RI.MATS.JC.txt',
}


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
    stored exon coordinates (626/626 exact across all five event types, zero
    off-by-one), so adding +1 would introduce a 1bp error rather than fix one. The
    matcher's 1bp tolerance is for exon- vs intron-boundary conventions, not to
    compensate for any base-offset here."""
    if event_type == 'SE':
        return [
            _ordered_pair(r['upstreamEE'], r['exonStart_0base']),  # inclusion: upstream -> exon
            _ordered_pair(r['exonEnd'], r['downstreamES']),        # inclusion: exon -> downstream
            _ordered_pair(r['upstreamEE'], r['downstreamES']),     # skipping
        ]
    if event_type == 'A5SS':
        return [
            _ordered_pair(r['longExonEnd'], r['flankingES']),      # long form
            _ordered_pair(r['shortEE'], r['flankingES']),          # short form
        ]
    if event_type == 'A3SS':
        return [
            _ordered_pair(r['flankingEE'], r['longExonStart_0base']),  # long form
            _ordered_pair(r['flankingEE'], r['shortES']),              # short form
        ]
    if event_type == 'MXE':
        return [
            _ordered_pair(r['upstreamEE'], r['1stExonStart_0base']),
            _ordered_pair(r['1stExonEnd'], r['downstreamES']),
            _ordered_pair(r['upstreamEE'], r['2ndExonStart_0base']),
            _ordered_pair(r['2ndExonEnd'], r['downstreamES']),
        ]
    if event_type == 'RI':
        return [
            _ordered_pair(r['upstreamEE'], r['downstreamES']),     # spliced (intron removed)
        ]
    return []


def rmats2junctions(rmats_dir):
    """Parse an rMATS-turbo output directory (the five [Event].MATS.JC.txt files)
    into a junctions DataFrame ready for JunctionsAnalysis.analyze_junctions().

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
            continue
        df = pd.read_csv(path, sep='\t')
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
            except (KeyError, ValueError, TypeError):
                continue  # malformed / incomplete row - skip it
            for start, end in event_junctions:
                junctions.append([chromosome, gene_ensembl_id, gene_symbol,
                                  event_type, start, end, cluster_name])

    if not junctions:
        raise ValueError(f"No rMATS MATS.JC.txt events found under {rmats_dir}")

    return pd.DataFrame(junctions, columns=[
        'chromosome', 'gene_ensembl_id', 'gene_symbol', 'event_type',
        'start_position', 'end_position', 'cluster_name'])


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

        event_junctions = _parse_coord_pairs(r['Junctions coords'])
        if ir_coords:
            event_junctions += _parse_coord_pairs(ir_coords)

        for start, end in event_junctions:
            junctions.append([chromosome, gene_ensembl_id, gene_symbol,
                              event_type, start, end, cluster_name])

    if not junctions:
        raise ValueError(f"No LSV junctions found in {tsv_path}")

    return pd.DataFrame(junctions, columns=[
        'chromosome', 'gene_ensembl_id', 'gene_symbol', 'event_type',
        'start_position', 'end_position', 'cluster_name'])


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
# to the training median. 5-fold CV AUC 0.768 and well-calibrated (predicted prob
# tracks observed pathogenic rate). Coefficients are on standardised features;
# region_am dominates (+0.80), gnomAD LOEUF is the non-circular add (-0.42, lower
# LOEUF = more constrained = more pathogenic), burial +0.24, coverage-loss ~0.
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
    """Calibrated probability (0-1) that a changed region is pathogenicity-relevant,
    from a logistic model over mean AlphaMissense (`region_am`), gnomAD gene
    constraint (`loeuf`), Pfam `max_cov_loss`, and AlphaFold `buried_frac`. This is
    the *continuous* companion to the categorical `hmm_change_impact`: it does not
    quantise, and it folds in the gene-level constraint. A missing feature is
    imputed to the training median (a neutral value), so the score degrades
    gracefully. See `_IMPACT_PROB_MODEL` for provenance and calibration.

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
# P(the alternative isoform has a different fold, TM-score < 0.5). Trained on
# 10,232 changed-regions labelled by the Genome Biology 2025 AlphaFold TM-scores,
# features standardised, missing imputed to the training median. Gene-grouped
# 5-fold CV AUC 0.777, accuracy 74% (@0.5), well-calibrated. This is a DIFFERENT
# question from impact_probability: burial flips sign - a buried changed region
# PRESERVES the fold here (coef -0.93) but harbours pathogenic variants there
# (+0.24). Regenerate with alphafold_benchmark/fit_foldchange.py.
_FOLD_CHANGE_MODEL = {
    'features': ['identity', 'buried_frac', 'mean_rsa', 'region_am', 'loeuf', 'max_cov_loss'],
    'median':   [75.0, 0.236, 0.4509, 0.483, 0.884, 11.0],
    'mean':     [66.3421, 0.2584, 0.4678, 0.481, 0.9085, 33.2563],
    'std':      [28.0247, 0.2339, 0.1819, 0.1855, 0.4535, 39.2996],
    'coef':     [-0.6662, -0.9294, 0.1828, 0.4646, 0.114, 0.537],
    'intercept': -0.8508,
}


def fold_change_probability(identity=None, buried_frac=None, mean_rsa=None,
                            region_am=None, loeuf=None, max_cov_loss=None):
    """Calibrated probability (0-1) that the alternative isoform adopts a DIFFERENT
    fold (AlphaFold TM-score < 0.5) - the structural companion to
    impact_probability (which is functional/pathogenicity). `identity` is the
    sequence identity (%) between the canonical and alternative protein. A missing
    feature is imputed to the training median. See `_FOLD_CHANGE_MODEL`.

    Returns a float in [0, 1], or None if every feature is missing.
    """
    m = _FOLD_CHANGE_MODEL
    vals = [identity, buried_frac, mean_rsa, region_am, loeuf, max_cov_loss]
    if all(v is None for v in vals):
        return None
    z = m['intercept']
    for i, v in enumerate(vals):
        x = m['median'][i] if v is None else v          # impute missing to the median
        z += m['coef'][i] * ((x - m['mean'][i]) / m['std'][i])
    return round(1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, z)))), 4)


def fold_change_call(fold_change_prob, low=0.4, high=0.6):
    """Triage a fold_change_probability into an actionable call. Cheap features
    top out at ~AUC 0.79 for fold change; the mid-probability band is genuinely
    ambiguous (~47% real change rate) and is the 'does the remainder refold' class
    that only actual folding resolves. So route it out rather than force a call:

        prob >= high  -> 'change'     (confident: different fold)
        prob <= low   -> 'preserved'  (confident: same fold)
        otherwise     -> 'uncertain'  (fold/inspect to resolve)

    Abstaining on 'uncertain' lifts accuracy on the called cases from ~75% (call
    all) to ~86% (call the confident half). Returns None if prob is None/blank.
    """
    if fold_change_prob is None or fold_change_prob == '':
        return None
    p = float(fold_change_prob)
    return 'change' if p >= high else 'preserved' if p <= low else 'uncertain'


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