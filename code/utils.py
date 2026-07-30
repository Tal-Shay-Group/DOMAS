import logging
import math
import os

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# rMATS-turbo per-event result files (junction-count variant). The JCEC files
# carry identical coordinates, so JC is enough for DOMAS's coordinate mapping.
#
# RI is deliberately absent. A retained-intron event yields exactly one junction
# (the spliced form - see _rmats_event_junctions()), and a single-junction cluster
# can never produce a comparable transcript: either the canonical transcript
# matches that junction, in which case every other transcript has no junction
# unique to it, or it doesn't, in which case the cluster stops at
# 'no_canonical_junctions'. Reading RI.MATS.JC.txt therefore only ever emitted
# unanalysable rows (verified: 880 rows over 40 clusters, none analysed), so the
# file is not read at all.
_RMATS_EVENT_FILES = {
    'SE': 'SE.MATS.JC.txt',
    'A5SS': 'A5SS.MATS.JC.txt',
    'A3SS': 'A3SS.MATS.JC.txt',
    'MXE': 'MXE.MATS.JC.txt',
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
    stored exon coordinates (626/626 exact across all event types, zero
    off-by-one), so adding +1 would introduce a 1bp error rather than fix one. The
    matcher's 1bp tolerance is for exon- vs intron-boundary conventions, not to
    compensate for any base-offset here.

    SE and MXE name their exons by *genomic* position - rMATS reports
    upstreamEE <= exonStart_0base and exonEnd <= downstreamES on both strands
    (verified 100% over the fixture) - so one set of pairs is correct either way.

    A5SS and A3SS instead name the long/short forms of one exon plus its flanking
    exon, so which genomic boundary distinguishes the two forms flips with strand:
    the alternative donor of an A5SS event is the exon's genomic END on the plus
    strand but its genomic START on the minus strand, and the flanking exon sits
    on the other side. Reading the plus-strand boundary regardless of strand made
    both forms yield the *same* junction, collapsing the event to a single
    junction - which can never produce a comparable transcript, so every
    minus-strand A5SS/A3SS event was silently unanalysable."""
    strand = str(r['strand']) if 'strand' in r else '+'
    if event_type == 'SE':
        return [
            _ordered_pair(r['upstreamEE'], r['exonStart_0base']),  # inclusion: upstream -> exon
            _ordered_pair(r['exonEnd'], r['downstreamES']),        # inclusion: exon -> downstream
            _ordered_pair(r['upstreamEE'], r['downstreamES']),     # skipping
        ]
    if event_type in ('A5SS', 'A3SS'):
        # The two alternative-splice-site types exchange roles under strand
        # reversal: a minus-strand A5SS has the same geometry as a plus-strand
        # A3SS. Built only in this branch - MXE has no longExon*/short*/flanking*
        # columns, and a KeyError here is swallowed by rmats2junctions().
        varying_boundary_is_upper = (event_type == 'A5SS') == (strand != '-')
        if varying_boundary_is_upper:
            return [
                _ordered_pair(r['longExonEnd'], r['flankingES']),          # long form
                _ordered_pair(r['shortEE'], r['flankingES']),              # short form
            ]
        return [
            _ordered_pair(r['flankingEE'], r['longExonStart_0base']),      # long form
            _ordered_pair(r['flankingEE'], r['shortES']),                  # short form
        ]
    if event_type == 'MXE':
        return [
            _ordered_pair(r['upstreamEE'], r['1stExonStart_0base']),
            _ordered_pair(r['1stExonEnd'], r['downstreamES']),
            _ordered_pair(r['upstreamEE'], r['2ndExonStart_0base']),
            _ordered_pair(r['2ndExonEnd'], r['downstreamES']),
        ]
    return []


def rmats2junctions(rmats_dir):
    """Parse an rMATS-turbo output directory (the SE/A5SS/A3SS/MXE [Event].MATS.JC.txt
    files; RI is skipped, see _RMATS_EVENT_FILES) into a junctions DataFrame ready
    for JunctionsAnalysis.analyze_junctions().

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
            logger.warning(f"rMATS {event_type} file not found at {path} - no {event_type} "
                           f"events will be analysed.")
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

    RI events are skipped, for the same reason RI.MATS.JC.txt is not read on the
    rMATS path (see _RMATS_EVENT_FILES): SUPPA writes a retained intron as
    `<gene>;RI:<chr>:<s1>:<e1>-<s2>:<e2>:<strand>`, in which only one token holds
    a junction, so the event yields a single junction - and a single-junction
    event can never produce a transcript differing from the canonical one within
    the event, so it is structurally unanalysable. Measured on
    events_RI_strict.ioe: 9,092 of 9,092 RI clusters had one junction, while
    every other SUPPA event type had none.

    Args:
        file_path (str): The path to the IOE file.
    """
    junctions = []
    df_ioe = pd.read_csv(file_path, sep='\t')
    count = 0
    skipped_ri = 0
    for row in df_ioe.itertuples():
        count += 1
        chromosome = row.seqname
        event_parts = row.event_id.split(';')
        gene_ensembl_id = event_parts[0]
        event_parts2 = event_parts[1].split(':')
        event_type = event_parts2[0]
        if event_type == 'RI':
            skipped_ri += 1
            continue
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


    if skipped_ri:
        logger.warning(f"{os.path.basename(file_path)}: skipped {skipped_ri} RI event(s) - a "
                       f"retained intron yields a single junction and cannot be analysed.")
    if not junctions:
        logger.warning(f"No analysable events found in {file_path}"
                       + (" (every event was RI)." if skipped_ri else "."))

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
# 10,227 changed-regions labelled by the Genome Biology 2025 AlphaFold TM-scores,
# features standardised, missing imputed to the training median.
#
# Feature set rebuilt around AlphaFold **PAE** (E38): `pae_global` (whole canonical
# structure mean predicted-aligned-error) is by far the strongest fold-change signal
# and dominates the model (std coef ~ +2). It lifts gene-grouped 5-fold CV from the
# old 6-feature AUC 0.777 to **AUC 0.894 / accuracy 0.815 / R2 0.659** (the full
# 10-feature model reaches 0.906; these four capture ~95% of the gain for a fraction
# of the provisioning). `identity` (trimmed-region sequence identity %) is the #2 term;
# `max_cov_loss` (max Pfam coverage lost) and `protL` (canonical protein length) add
# small real signal. Burial / region_am / loeuf were dropped: ~0 marginal here once
# PAE is in (they drive impact_probability, the functional axis - the fold-vs-function
# duality). Regenerate with alphafold_benchmark/fit_foldchange_pae.py.
_FOLD_CHANGE_MODEL = {
    'features': ['pae_global', 'identity', 'max_cov_loss', 'protL'],
    'median':   [17.1268, 75.0, 11.0, 388.0],
    'mean':     [16.5883, 66.3577, 33.2236, 379.1839],
    'std':      [7.1173, 28.0213, 39.2834, 134.3642],
    'coef':     [1.9553, -0.6826, 0.6514, -0.5351],
    'intercept': -1.2597,
}


def fold_change_probability(pae_global=None, identity=None, max_cov_loss=None, protL=None):
    """Calibrated probability (0-1) that the alternative isoform adopts a DIFFERENT
    fold (AlphaFold TM-score < 0.5) - the structural companion to
    impact_probability (which is functional/pathogenicity).

    Features (see `_FOLD_CHANGE_MODEL`):
      pae_global   - mean predicted-aligned-error over the canonical AlphaFold
                     structure (higher = floppier/multi-domain -> more fold change).
                     The dominant signal; provisioned via the afdb_pae table.
      identity     - sequence identity (%) between canonical and alternative protein.
      max_cov_loss - max Pfam coverage (%) lost in the changed region.
      protL        - canonical protein length (aa).
    A missing feature is imputed to the training median.

    Returns a float in [0, 1], or None if every feature is missing.
    """
    m = _FOLD_CHANGE_MODEL
    vals = [pae_global, identity, max_cov_loss, protL]
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

# ---------------------------------------------------------------------------
# DoChaP database readers.
#
# These live here rather than in alternative_splicing.py so that
# junction_analisys.py can import them directly: alternative_splicing.py
# imports junction_analisys, so reaching them from there required a
# function-local import to dodge the circular dependency.
# ---------------------------------------------------------------------------

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
    # domain_id + InterPro entry `type` are carried through so
    # junction_analisys.filter_representative_domains() can reduce the domain
    # set by curated type (Domain/Repeat vs Family/Homologous_superfamily)
    # instead of the geometric collapse heuristic. `type` is NULL for DBs built
    # before RepresentativeDomains gained the column (filter then no-ops).
    'domain_id', 'type',
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

    # `type` is absent in DBs built before RepresentativeDomains gained the
    # column; add it as NULL so downstream columns/filtering stay uniform.
    if 'type' not in df_rep.columns:
        df_rep['type'] = None

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
