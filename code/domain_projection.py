"""Project a canonical transcript's domain coordinates onto an alternative isoform.

Why this exists
---------------
DoChaP's ``RepresentativeDomains`` table is keyed on one column,
``protein_interpro_id`` (a UniProt accession), and nothing else. There is no
transcript or protein column in it. ``RepresentativeDomainsBuilder`` assigns a
gene's Swiss-Prot accession to every transcript UniProt's idmapping lists for
it, and its collision rule deliberately prefers that reviewed accession over any
isoform-specific TrEMBL one. The consequence is that several isoforms of
*different lengths* read back the SAME domain rows at the SAME amino-acid
coordinates - the canonical's.

Measured on DB_merged.sqlite (Aug 2025 build):

  * 530,604 proteins carry an accession, over 347,049 distinct accessions
  * 15,129 accessions are shared by proteins of differing lengths, covering
    92,314 proteins (~17% of all mapped proteins)
  * 85,362 stored domain rows fall past their protein's own C-terminus, on
    39,684 distinct proteins; 4,354 proteins carry a domain that *starts* past
    their end

Clipping to protein length (``utils._clip_domains_to_protein``) only catches the
C-terminal overruns. It cannot catch a domain that sits comfortably inside the
shorter protein but is encoded by an exon the shorter protein does not contain.
ARAP1 is exactly that case: the SAM domain at aa 3-70 is encoded entirely within
exon chr11:72726619-72727172, which ``ENST00000334211`` does not contain at all -
yet the stored rows place SAM on it, inside its 1205 aa, where no length check
can see the problem.

What this module does
---------------------
It re-derives each domain's position on the alternative isoform from the two
transcripts' own CDS structures, which ``Transcript_Exon`` records exactly:

  * ``abs_start_CDS`` / ``abs_end_CDS`` are 1-based, *contiguous* CDS nucleotide
    offsets; their maximum equals ``protein length * 3 + 3`` (the stop codon) -
    a relation ``check_db.py`` already asserts for canonical transcripts.
  * ``genomic_start_tx`` is BED-style (0-based). An exon's length is
    ``genomic_end_tx - genomic_start_tx``, NOT ``+ 1``; its 1-based closed
    interval is ``[genomic_start_tx + 1, genomic_end_tx]``. This is verifiable
    per transcript: every internal (fully coding) exon has an abs-CDS span
    exactly one shorter than ``end - start + 1``.

Mapping CDS nucleotides to genomic positions on both transcripts and
intersecting gives, for every canonical domain, the part of it the alternative
actually encodes - and in which reading frame. A stretch of shared DNA read in a
different frame yields a different peptide, so it does not preserve the domain;
that case is reported separately rather than silently counted as retained.

Scope
-----
Only isoforms that *inherited* the canonical's rows need projecting. A protein
whose ``protein_interpro_id`` differs from the canonical's got its own
isoform-specific accession, so its stored coordinates are genuinely its own and
must be left alone. ``needs_projection()`` encodes that test.
"""

from __future__ import annotations

from collections import namedtuple

import pandas as pd

__all__ = [
    'KEPT', 'TRUNCATED', 'FRAMESHIFTED', 'ABSENT',
    'CdsMap', 'build_cds_map', 'needs_projection',
    'ProjectedDomain', 'project_domain', 'project_domains',
]

# ---------------------------------------------------------------------------
# Projection outcomes
# ---------------------------------------------------------------------------

KEPT = 'kept'                  # every coding nucleotide retained, frame intact
TRUNCATED = 'truncated'        # part retained in frame, part lost
FRAMESHIFTED = 'frameshifted'  # DNA retained but read in a different frame
ABSENT = 'absent'              # the alternative encodes none of it

#: A contiguous run of CDS whose genomic<->CDS-nt relation is affine.
#: ``nt(g) = base + sign * g`` for any genomic position g in [glo, ghi],
#: where ``sign`` is +1 on the plus strand and -1 on the minus strand.
_Block = namedtuple('_Block', 'glo ghi base')


class CdsMap:
    """An affine, strand-aware map between a transcript's CDS nucleotides and
    genomic coordinates.

    The key property exploited by :func:`project_domain` is that within any
    (canonical block, alternative block) pair the difference between the two
    transcripts' CDS offsets at the same genomic position is *constant*::

        n_can(g) - n_alt(g) = base_can - base_alt

    so the reading-frame relationship between the two transcripts over a shared
    stretch of DNA is a single number, not a per-position computation.
    """

    __slots__ = ('blocks', 'sign', 'strand', 'cds_length')

    def __init__(self, blocks, sign, strand, cds_length):
        self.blocks = blocks
        self.sign = sign
        self.strand = strand
        self.cds_length = cds_length

    def __repr__(self):  # pragma: no cover - debugging aid
        return (f'CdsMap(strand={self.strand!r}, blocks={len(self.blocks)}, '
                f'cds_length={self.cds_length})')

    def nt_at(self, genomic_position):
        """CDS nucleotide offset at a genomic position, or None if not coding."""
        for block in self.blocks:
            if block.glo <= genomic_position <= block.ghi:
                return block.base + self.sign * genomic_position
        return None

    def genomic_intervals(self, nt_start, nt_end):
        """The genomic intervals encoding CDS nucleotides ``[nt_start, nt_end]``.

        Returns a list of ``(glo, ghi, base)`` with ``glo <= ghi``.
        """
        pieces = []
        for block in self.blocks:
            lo_nt = min(block.base + self.sign * block.glo,
                        block.base + self.sign * block.ghi)
            hi_nt = max(block.base + self.sign * block.glo,
                        block.base + self.sign * block.ghi)
            m1 = max(nt_start, lo_nt)
            m2 = min(nt_end, hi_nt)
            if m1 > m2:
                continue
            g1 = self.sign * (m1 - block.base)
            g2 = self.sign * (m2 - block.base)
            pieces.append((min(g1, g2), max(g1, g2), block.base))
        return pieces


def build_cds_map(df_exons, strand, cds_start=None, cds_end=None):
    """Build a :class:`CdsMap` from this transcript's ``Transcript_Exon`` rows.

    Parameters
    ----------
    df_exons : DataFrame
        Rows for ONE transcript, carrying ``genomic_start_tx``,
        ``genomic_end_tx``, ``abs_start_CDS`` and ``abs_end_CDS``.
    strand : {'+', '-', 1, -1}
        The gene's strand (``Genes.strand``).
    cds_start, cds_end : int, optional
        ``Transcripts.cds_start`` / ``cds_end``. Only consulted when the whole
        CDS lies in a single exon, where the exon boundaries alone cannot say
        where within it the coding part sits.

    Returns
    -------
    CdsMap or None
        None when the transcript has no coding exons, or when the single-exon
        case cannot be resolved. Returning None (rather than guessing) is
        deliberate: a caller that cannot project must fall back to *not* drawing
        inherited domains, never to drawing them at the wrong place.

    Notes
    -----
    Exon placement follows from the CDS offsets being contiguous. Reading the
    transcript 5'->3', every coding exon after the first is flush with its
    transcript-5' boundary, because its first CDS nucleotide continues directly
    from the previous exon's last one. The *first* coding exon is the exception:
    the 5'UTR occupies its transcript-5' side, so its coding block is flush with
    its transcript-3' boundary instead. The last coding exon needs no special
    case - its 3'UTR trails off the transcript-3' end, leaving the block flush
    at the 5' side like every other non-first exon.
    """
    sign = -1 if str(strand) in ('-', '-1') else 1

    if df_exons is None or len(df_exons) == 0:
        return None
    if 'abs_start_CDS' not in df_exons.columns:
        return None

    coding = df_exons[
        pd.to_numeric(df_exons['abs_start_CDS'], errors='coerce').fillna(0) > 0
    ].copy()
    if len(coding) == 0:
        return None

    coding['abs_start_CDS'] = coding['abs_start_CDS'].astype(int)
    coding['abs_end_CDS'] = coding['abs_end_CDS'].astype(int)
    coding = coding[coding['abs_end_CDS'] >= coding['abs_start_CDS']]
    if len(coding) == 0:
        return None
    coding = coding.sort_values('abs_start_CDS').reset_index(drop=True)

    blocks = []
    for position, exon in coding.iterrows():
        # genomic_start_tx is BED-style: the 1-based closed exon is
        # [genomic_start_tx + 1, genomic_end_tx].
        exon_lo = int(exon['genomic_start_tx']) + 1
        exon_hi = int(exon['genomic_end_tx'])
        a1 = int(exon['abs_start_CDS'])
        a2 = int(exon['abs_end_CDS'])
        n_cds = a2 - a1 + 1

        if n_cds > (exon_hi - exon_lo + 1):
            # The CDS span claims more nucleotides than the exon holds; the row
            # is inconsistent and any placement would be invented.
            return None

        first_coding = (position == 0)
        last_coding = (position == len(coding) - 1)

        if first_coding and last_coding:
            # Single-exon CDS: UTR on both sides, so neither boundary anchors it.
            if cds_start is None or cds_end is None:
                return None
            block_lo = min(int(cds_start), int(cds_end))
            block_hi = max(int(cds_start), int(cds_end))
            block_lo = max(block_lo, exon_lo)
            block_hi = min(block_hi, exon_hi)
            if block_hi - block_lo + 1 < n_cds:
                return None
            # Trim from the transcript-3' side, where the stop codon and 3'UTR sit.
            if sign < 0:
                block_lo = block_hi - n_cds + 1
            else:
                block_hi = block_lo + n_cds - 1
        elif first_coding:
            # Flush with the transcript-3' boundary; 5'UTR takes the other side.
            if sign < 0:
                block_lo, block_hi = exon_lo, exon_lo + n_cds - 1
            else:
                block_lo, block_hi = exon_hi - n_cds + 1, exon_hi
        else:
            # Flush with the transcript-5' boundary.
            if sign < 0:
                block_lo, block_hi = exon_hi - n_cds + 1, exon_hi
            else:
                block_lo, block_hi = exon_lo, exon_lo + n_cds - 1

        # nt(g) = base + sign * g, fixed by the block's transcript-5' end.
        g_at_a1 = block_hi if sign < 0 else block_lo
        base = a1 - sign * g_at_a1
        blocks.append(_Block(block_lo, block_hi, base))

    cds_length = int(coding['abs_end_CDS'].max())
    return CdsMap(blocks, sign, '-' if sign < 0 else '+', cds_length)


def needs_projection(protein_interpro_id, canonical_interpro_id, is_canonical):
    """True when this protein's stored domain rows were inherited, not its own.

    A protein that shares the canonical's accession is reading the canonical's
    rows out of ``RepresentativeDomains`` - there is only one set of rows per
    accession. A protein with a different accession has its own isoform-specific
    entry, so its coordinates are real and must be left untouched.
    """
    if is_canonical:
        return False
    if not protein_interpro_id or not canonical_interpro_id:
        return False
    return str(protein_interpro_id).strip() == str(canonical_interpro_id).strip()


#: The outcome of projecting one canonical domain onto one alternative isoform.
#:
#: ``aa_start``/``aa_end`` are the domain's position on the ALTERNATIVE protein
#: (None when absent). ``retained_fraction`` is the share of the domain's coding
#: nucleotides the alternative encodes in the canonical's reading frame.
ProjectedDomain = namedtuple(
    'ProjectedDomain',
    'status aa_start aa_end retained_fraction retained_nt frameshifted_nt domain_nt',
)


def project_domain(aa_start, aa_end, canonical_map, alt_map, alt_protein_length=None):
    """Project one domain from canonical amino-acid coordinates onto an isoform.

    Parameters
    ----------
    aa_start, aa_end : int
        The domain's 1-based inclusive position on the CANONICAL protein.
    canonical_map, alt_map : CdsMap
    alt_protein_length : int, optional
        Used only to clip the reported end, so a projected coordinate can never
        exceed the alternative's own length.
    """
    aa_start = int(aa_start)
    aa_end = int(aa_end)
    if aa_end < aa_start:
        aa_start, aa_end = aa_end, aa_start

    nt_start = (aa_start - 1) * 3 + 1
    nt_end = aa_end * 3
    domain_nt = nt_end - nt_start + 1

    if canonical_map is None or alt_map is None:
        return ProjectedDomain(ABSENT, None, None, 0.0, 0, 0, domain_nt)

    canonical_pieces = canonical_map.genomic_intervals(nt_start, nt_end)

    retained_nt = 0
    frameshifted_nt = 0
    alt_nt_lo = None
    alt_nt_hi = None

    for glo, ghi, base_can in canonical_pieces:
        for block in alt_map.blocks:
            olo = max(glo, block.glo)
            ohi = min(ghi, block.ghi)
            if olo > ohi:
                continue
            overlap_nt = ohi - olo + 1
            # Both transcripts are on the same strand, so the CDS-offset
            # difference over a shared stretch is a single constant.
            if (base_can - block.base) % 3 != 0:
                frameshifted_nt += overlap_nt
                continue
            retained_nt += overlap_nt
            n1 = block.base + alt_map.sign * olo
            n2 = block.base + alt_map.sign * ohi
            lo, hi = min(n1, n2), max(n1, n2)
            alt_nt_lo = lo if alt_nt_lo is None else min(alt_nt_lo, lo)
            alt_nt_hi = hi if alt_nt_hi is None else max(alt_nt_hi, hi)

    if retained_nt == 0:
        status = FRAMESHIFTED if frameshifted_nt > 0 else ABSENT
        return ProjectedDomain(status, None, None, 0.0,
                               retained_nt, frameshifted_nt, domain_nt)

    new_start = (alt_nt_lo - 1) // 3 + 1
    new_end = (alt_nt_hi - 1) // 3 + 1
    if alt_protein_length:
        new_end = min(int(new_end), int(alt_protein_length))
        new_start = min(int(new_start), int(alt_protein_length))
    if new_end < new_start:
        new_end = new_start

    status = KEPT if retained_nt == domain_nt else TRUNCATED
    return ProjectedDomain(status, int(new_start), int(new_end),
                           retained_nt / domain_nt,
                           retained_nt, frameshifted_nt, domain_nt)


def project_domains(df_domains, canonical_map, alt_map, alt_protein_length=None,
                    drop_absent=True):
    """Project a frame of canonical domains onto an alternative isoform.

    ``df_domains`` carries ``AA_start``/``AA_end`` as produced by
    ``generate_gene_pdf._representative_domains_to_domain_columns``. The returned
    frame keeps every other column untouched, rewrites those two, and adds:

      ``projection_status``    one of kept / truncated / frameshifted / absent
      ``retained_fraction``    share of the domain encoded in frame (0.0 - 1.0)
      ``canonical_AA_start``   the inherited coordinates, kept for reference
      ``canonical_AA_end``

    With ``drop_absent=True`` (the default) domains the isoform does not encode
    are removed, which is what a figure should draw. Pass False to keep them for
    reporting.
    """
    if df_domains is None or len(df_domains) == 0:
        return df_domains

    records = []
    for _, domain in df_domains.iterrows():
        projected = project_domain(domain['AA_start'], domain['AA_end'],
                                   canonical_map, alt_map, alt_protein_length)
        row = domain.copy()
        row['canonical_AA_start'] = domain['AA_start']
        row['canonical_AA_end'] = domain['AA_end']
        row['projection_status'] = projected.status
        row['retained_fraction'] = projected.retained_fraction
        row['AA_start'] = projected.aa_start
        row['AA_end'] = projected.aa_end
        records.append(row)

    result = pd.DataFrame(records).reset_index(drop=True)
    if drop_absent:
        keep = result['projection_status'].isin([KEPT, TRUNCATED])
        result = result[keep].reset_index(drop=True)
    return result
