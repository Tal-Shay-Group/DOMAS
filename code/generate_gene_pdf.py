"""
Generate gene visualization PDF similar to DoChap-web JavaScript visualization.
This module creates PDF reports showing gene transcripts with exons and protein domains.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle, FancyBboxPatch
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import sqlite3


# Color palette from DoChap-web
DOCHAP_COLORS = [
    "#FF4A46", "#1BE177", "#00CCFF", "#A30059", "#FFB500", "#006FA6", "#4FC601", "#FFDBE5",
    "#D16100", "#00C2A0", "#A079BF", "#C0B9B2", "#CC0744", "#549E79", "#B79762", "#B903AA",
    "#00846F", "#FF90C9", "#0AA6D8", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#BEC459",
    "#AA5199", "#0089A3", "#EEC3FF", "#8FB0FF", "#004D43", "#F4D749", "#997D87", "#3B5DFF",
    "#FF2F80", "#6B7900", "#FFAA92", "#A1C299", "#885578", "#B77B68", "#FAD09F", "#456D75",
    "#FF8A9A", "#0086ED", "#D157A0", "#00A6AA", "#B4A8BD", "#FF913F", "#636375", "#A3C8C9"
]

DOMAIN_LABEL_MAX_LEN = 14
PDF_RASTER_DPI = 110
DOMAIN_STRIP_COUNT = 80


class GeneVisualization:
    """Class to create gene visualization similar to DoChap-web."""
    
    def __init__(self, conn, gene_name):
        """
        Initialize gene visualization.
        
        Parameters:
        -----------
        conn : sqlite3.Connection
            Database connection to DoChap database
        gene_name : str
            Gene symbol to visualize
        """
        self.conn = conn
        self.gene_name = gene_name
        self.gene_data = None
        self.transcripts = []
        self.colors = {}
        self.color_index = 0
        
    def load_gene_data(self):
        """Load gene data from database."""
        # Get gene info
        gene_query = """
            SELECT * FROM Genes 
            WHERE gene_symbol = ? COLLATE NOCASE
        """
        df_gene = pd.read_sql_query(gene_query, self.conn, params=[self.gene_name])
        
        if len(df_gene) == 0:
            raise ValueError(f"Gene '{self.gene_name}' not found in database")
        
        self.gene_data = df_gene.iloc[0]
        
        # Get transcripts for this gene
        trans_query = """
            SELECT * FROM Transcripts 
            WHERE gene_ensembl_id = ?
            ORDER BY tx_start
        """
        df_transcripts = pd.read_sql_query(trans_query, self.conn, 
                                           params=[self.gene_data['gene_ensembl_id']])
        
        # Load exons and domains for each transcript
        for _, transcript in df_transcripts.iterrows():
            transcript_data = {
                'info': transcript,
                'exons': self._load_exons(transcript['transcript_ensembl_id']),
                'domains': self._load_domains(transcript['transcript_ensembl_id'])
            }
            self.transcripts.append(transcript_data)
        
        # Assign colors to exons
        self._assign_exon_colors()
        
    def _load_exons(self, transcript_id):
        """Load exons for a transcript."""
        exon_query = """
            SELECT * FROM Transcript_exon 
            WHERE transcript_ensembl_id = ?
            ORDER BY abs_start_CDS
        """
        df_exons = pd.read_sql_query(exon_query, self.conn, params=[transcript_id])
        # If any exons have abs_start_CDS = 0 (non-coding), sort those by genomic position
        # and put them at the beginning
        if (df_exons['abs_start_CDS'] == 0).any():
            non_coding = df_exons[df_exons['abs_start_CDS'] == 0].sort_values('genomic_start_tx')
            coding = df_exons[df_exons['abs_start_CDS'] > 0].sort_values('abs_start_CDS')
            df_exons = pd.concat([non_coding, coding], ignore_index=True)
        return df_exons
    
    def _load_domains(self, transcript_id):
        """Load protein domains for a transcript."""
        # Get protein ID
        protein_query = """
            SELECT protein_ensembl_id FROM Proteins 
            WHERE transcript_ensembl_id = ?
        """
        df_protein = pd.read_sql_query(protein_query, self.conn, params=[transcript_id])
        
        if len(df_protein) == 0:
            return pd.DataFrame()
        
        protein_id = df_protein.iloc[0]['protein_ensembl_id']
        
        # Get domain events
        domain_query = """
            SELECT de.*, dt.*
            FROM DomainEvent de
            JOIN DomainType dt ON de.type_id = dt.type_id
            WHERE de.protein_ensembl_id = ?
            ORDER BY de.AA_start
        """
        return pd.read_sql_query(domain_query, self.conn, params=[protein_id])
    
    def _assign_exon_colors(self):
        """Assign colors to unique exons across all transcripts based on genomic location."""
        seen_exons = {}
        
        for transcript in self.transcripts:
            for _, exon in transcript['exons'].iterrows():
                key = (int(exon['genomic_start_tx']), int(exon['genomic_end_tx']))
                
                if key not in seen_exons:
                    seen_exons[key] = DOCHAP_COLORS[self.color_index % len(DOCHAP_COLORS)]
                    self.color_index += 1
        
        self.colors = seen_exons
    
    def _get_exon_color(self, exon):
        """Get color for an exon, keyed on genomic location."""
        key = (int(exon['genomic_start_tx']), int(exon['genomic_end_tx']))

        color = self.colors.get(key)
        gray_like = {'#cccccc', '#808080', '#a9a9a9', 'gray', 'grey', 'darkgray', 'darkgrey'}
        if color is not None and str(color).strip().lower() not in gray_like:
            return color

        # Deterministic non-gray fallback.
        idx = abs(hash(key)) % len(DOCHAP_COLORS)
        return DOCHAP_COLORS[idx]

    def _format_domain_label(self, domain_row, compact_mode=False, max_len=DOMAIN_LABEL_MAX_LEN):
        """Build a compact domain label; append * when truncated or combined."""
        names = []
        for col in ['interpro', 'pfam', 'cdd', 'smart']:
            if col in domain_row and pd.notna(domain_row[col]):
                value = str(domain_row[col]).strip()
                if value and value.lower() != 'nan' and value not in names:
                    names.append(value)

        if not names:
            return None

        primary = names[0]
        needs_star = compact_mode or len(names) > 1 or len(primary) > max_len
        if len(primary) > max_len:
            primary = primary[:max_len]

        return f"{primary}*" if needs_star else primary
    
    def _format_number(self, num):
        """Format number with commas."""
        return f"{int(num):,}"

    def _format_axis_number(self, num):
        """Format axis numbers compactly to save vertical space on genomic scales."""
        value = float(num)
        abs_value = abs(value)
        if abs_value >= 1_000_000:
            return f"{value / 1_000_000:.4g}M"
        if abs_value >= 1_000:
            return f"{value / 1_000:.4g}K"
        return str(int(round(value)))

    def _is_negative_strand(self):
        """Return True when gene strand indicates reverse genomic orientation."""
        if self.gene_data is None or 'strand' not in self.gene_data:
            return False
        strand_value = self.gene_data['strand']
        if pd.isna(strand_value):
            return False

        strand_text = str(strand_value).strip().lower()
        return strand_text == '-' or strand_text.startswith('-') or strand_text == 'minus'

    def _transcript_produces_protein(self, transcript):
        """Return True when transcript has a protein identifier."""
        protein_id = transcript['info'].get('protein_ensembl_id')
        if protein_id is None or pd.isna(protein_id):
            return False
        protein_text = str(protein_id).strip()
        return bool(protein_text) and protein_text.lower() != 'nan'

    def _transcript_has_domains(self, transcript):
        """Return True when transcript protein has at least one domain row."""
        return len(transcript['domains']) > 0

    def _compute_domain_label_positions(self, label_items, axis_max, base_y, lane_step=0.08, lanes=4):
        """Assign below-domain label positions across lanes to reduce overlap."""
        if not label_items:
            return []

        # Keep deterministic ordering for stable output.
        ordered = sorted(label_items, key=lambda item: item['center'])
        lane_right_edges = [-1e12] * lanes
        placed = []

        for item in ordered:
            # Approximate label footprint in axis units from text length.
            est_half_width = max(item['width'] * 0.25, len(item['text']) * 2.2)
            left_edge = item['center'] - est_half_width
            right_edge = item['center'] + est_half_width

            lane_idx = None
            for candidate in range(lanes):
                if left_edge >= lane_right_edges[candidate]:
                    lane_idx = candidate
                    break

            if lane_idx is None:
                lane_idx = min(range(lanes), key=lambda idx: lane_right_edges[idx])

            lane_right_edges[lane_idx] = right_edge + max(3.0, axis_max * 0.004)
            placed.append({
                'center': item['center'],
                'text': item['text'],
                'label_y': base_y - lane_idx * lane_step,
            })

        return placed

    def _select_longest_labels_for_overlaps(self, domain_entries):
        """Pick one label per overlap cluster: the domain with the longest AA span."""
        if not domain_entries:
            return []

        ordered = sorted(domain_entries, key=lambda entry: entry['start'])
        selected = []

        cluster = [ordered[0]]
        cluster_end = ordered[0]['end']

        for entry in ordered[1:]:
            if entry['start'] <= cluster_end:
                cluster.append(entry)
                cluster_end = max(cluster_end, entry['end'])
            else:
                longest = max(cluster, key=lambda item: item['aa_span'])
                selected.append({
                    'center': longest['center'],
                    'width': longest['width'],
                    'domain': longest['domain'],
                    'has_overlap': len(cluster) > 1,
                })
                cluster = [entry]
                cluster_end = entry['end']

        longest = max(cluster, key=lambda item: item['aa_span'])
        selected.append({
            'center': longest['center'],
            'width': longest['width'],
            'domain': longest['domain'],
            'has_overlap': len(cluster) > 1,
        })

        return selected

    def _get_coding_exon_segments(self, transcript):
        """Return coding exon segments in amino-acid coordinate space."""
        segments = []
        for exon_idx, (_, exon) in enumerate(transcript['exons'].iterrows()):
            if exon['abs_start_CDS'] > 0 and exon['abs_end_CDS'] >= exon['abs_start_CDS']:
                start_nuc = int(exon['abs_start_CDS'])
                end_nuc = int(exon['abs_end_CDS'])
                # Convert nucleotide CDS coordinates to amino-acid axis.
                start_aa = (start_nuc - 1) / 3.0
                end_aa = end_nuc / 3.0
                color = self._get_exon_color(exon)
                if color == '#CCCCCC':
                    color = DOCHAP_COLORS[exon_idx % len(DOCHAP_COLORS)]
                segments.append({
                    'start_aa': start_aa,
                    'end_aa': end_aa,
                    'color': color,
                })
        return sorted(segments, key=lambda segment: segment['start_aa'])

    def _prepare_junction_display_df(self, df_junction):
        """Return a copy of junction metadata with stable idx/color columns."""
        if df_junction is None or len(df_junction) == 0:
            return None

        display_df = df_junction.copy().reset_index(drop=True)
        display_df['idx'] = np.arange(1, len(display_df) + 1)
        display_df['junction_color'] = [
            DOCHAP_COLORS[i % len(DOCHAP_COLORS)] for i in range(len(display_df))
        ]
        return display_df

    def _get_matching_junctions(self, transcript, df_junction):
        """Return matched junction rows with coordinates, idx, and display color.

        Draw criteria per transcript:
          1) Breakpoints sit at exon boundaries: left breakpoint at an exon end and right
              breakpoint at an exon start, with no fully contained exon in between.
        2) At least one breakpoint is inside an exon (not at exon boundary).
        """
        if df_junction is None or len(df_junction) == 0:
            return []
        if 'start' not in df_junction.columns or 'end' not in df_junction.columns:
            return []

        exon_ranges = []
        exon_starts = set()
        exon_ends = set()
        for _, exon in transcript['exons'].iterrows():
            exon_start = int(min(exon['genomic_start_tx'], exon['genomic_end_tx']))
            exon_end = int(max(exon['genomic_start_tx'], exon['genomic_end_tx']))
            exon_ranges.append((exon_start, exon_end))
            exon_starts.add(exon_start)
            exon_ends.add(exon_end)

        matches = []
        for _, junction in df_junction.iterrows():
            if pd.isna(junction['start']) or pd.isna(junction['end']):
                continue
            start = int(junction['start'])
            end = int(junction['end'])

            left = min(start, end)
            right = max(start, end)

            # Rule 1: boundary pair with zero complete exon(s) contained between breakpoints.
            left_is_exon_end = left in exon_ends
            right_is_exon_start = right in exon_starts
            boundary_pair = left_is_exon_end and right_is_exon_start
            contains_complete_exon = any(
                exon_start > left and exon_end < right
                for exon_start, exon_end in exon_ranges
            )
            criterion_boundary = boundary_pair and not contains_complete_exon

            # Rule 2: either breakpoint is strictly inside any exon (not first/last base).
            left_inside_exon = any(
                exon_start < left < exon_end for exon_start, exon_end in exon_ranges
            )
            right_inside_exon = any(
                exon_start < right < exon_end for exon_start, exon_end in exon_ranges
            )
            criterion_inside_exon = left_inside_exon or right_inside_exon

            if criterion_boundary or criterion_inside_exon:
                matches.append({
                    'start': start,
                    'end': end,
                    'idx': int(junction['idx']) if 'idx' in junction and pd.notna(junction['idx']) else None,
                    'color': junction['junction_color'] if 'junction_color' in junction else 'red',
                })

        matches.sort(key=lambda item: (item['start'], item['end'], item['idx'] if item['idx'] is not None else 0))
        return matches

    def _draw_junction_table(self, ax, df_junction):
        """Draw a compact first-page table for provided junction metadata."""
        ax.axis('off')
        if df_junction is None or len(df_junction) == 0:
            return

        table_df = df_junction.copy()
        if 'junction_color' in table_df.columns:
            junction_colors = list(table_df['junction_color'])
            table_df = table_df.drop(columns=['junction_color'])
        else:
            junction_colors = ['black'] * len(table_df)
        if 'idx' in table_df.columns:
            idx_column = list(table_df.columns).index('idx')
            front_cols = ['idx'] + [col for col in table_df.columns if col != 'idx']
            table_df = table_df[front_cols]
            idx_column = 0
        else:
            idx_column = None
        if 'cluster_name' in table_df.columns:
            table_df = table_df.drop(columns=['cluster_name'])
        if table_df.empty:
            return

        display_df = table_df.fillna('')
        column_names = list(display_df.columns)
        table = ax.table(
            cellText=display_df.astype(str).values,
            colLabels=column_names,
            loc='center',
            cellLoc='center',
        )
        table.auto_set_font_size(False)
        table.set_fontsize(7)
        table.scale(1, 1.15)

        width_map = {
            'idx': 0.035,
            'start': 0.095,
            'end': 0.095,
        }
        flexible_columns = [name for name in column_names if name not in width_map]
        remaining_width = max(0.10, 1.0 - sum(width_map.get(name, 0.0) for name in column_names))
        default_width = remaining_width / max(1, len(flexible_columns))

        for (row, col), cell in table.get_celld().items():
            column_name = column_names[col]
            cell.set_width(width_map.get(column_name, default_width))
            cell.set_linewidth(0.6)
            cell.set_edgecolor('#888888')
            if row == 0:
                cell.set_facecolor('#F1F1F1')
                cell.set_text_props(weight='bold')
            elif idx_column is not None and col == idx_column:
                cell.set_text_props(color=junction_colors[row - 1], weight='bold')

    def _draw_genomic_junctions(self, ax, junctions, exon_y, exon_height):
        """Draw colored bracket-like junction markers above the genomic exon track."""
        if not junctions:
            return

        baseline_top = exon_y + exon_height / 2
        for index, junction in enumerate(junctions):
            left = min(junction['start'], junction['end'])
            right = max(junction['start'], junction['end'])
            junction_top = min(0.95, baseline_top + 0.10 + index * 0.08)
            color = junction['color']
            ax.plot([left, left], [baseline_top, junction_top], color=color, linewidth=1.4, zorder=5)
            ax.plot([right, right], [baseline_top, junction_top], color=color, linewidth=1.4, zorder=5)
            ax.plot([left, right], [junction_top, junction_top], color=color, linewidth=1.4, zorder=5)
    
    def create_pdf(self, output_file='gene_visualization.pdf', transcripts_per_page=4,
                   protein_only=False, domains_only=False, df_junction=None):
        """
        Create PDF visualization of the gene, one page per transcripts_per_page transcripts.
        Each page has its own axis scales at the top so they never overlap with transcript rows.

        Parameters:
        -----------
        output_file : str
            Path to save the PDF.
        transcripts_per_page : int
            Number of transcript rows to draw per page.
        protein_only : bool
            If True, include only transcripts with a protein identifier.
        domains_only : bool
            If True, include only transcripts whose protein has at least one domain.
        df_junction : pandas.DataFrame | None
            Optional junction metadata table with at least `start` and `end` columns.
        """
        from matplotlib.backends.backend_pdf import PdfPages

        if self.gene_data is None:
            self.load_gene_data()

        junction_display_df = self._prepare_junction_display_df(df_junction)

        # Skip empty/invalid transcript entries
        valid_transcripts = [t for t in self.transcripts if len(t['exons']) > 0]
        if protein_only:
            valid_transcripts = [t for t in valid_transcripts if self._transcript_produces_protein(t)]
        if domains_only:
            valid_transcripts = [t for t in valid_transcripts if self._transcript_has_domains(t)]
        num_transcripts = len(valid_transcripts)
        if num_transcripts == 0:
            if protein_only and domains_only:
                print(
                    f"No transcripts with protein and domains found for gene {self.gene_name}"
                )
            elif domains_only:
                print(f"No domain-containing transcripts found for gene {self.gene_name}")
            elif protein_only:
                print(f"No protein-producing transcripts found for gene {self.gene_name}")
            else:
                print(f"No transcripts found for gene {self.gene_name}")
            return

        genomic_start = min(t['info']['tx_start'] for t in valid_transcripts)
        genomic_end   = max(t['info']['tx_end']   for t in valid_transcripts)

        # Calculate max AA axis from all transcripts once.
        max_protein_length = 0
        for transcript in valid_transcripts:
            protein_length = transcript['info'].get('protein_length')
            if protein_length is not None and not pd.isna(protein_length):
                max_protein_length = max(max_protein_length, float(protein_length))
            coding_segments = self._get_coding_exon_segments(transcript)
            if coding_segments:
                max_protein_length = max(max_protein_length,
                                         max(s['end_aa'] for s in coding_segments))
            if len(transcript['domains']) > 0:
                max_protein_length = max(max_protein_length,
                                         float(transcript['domains']['AA_end'].max()))
        if max_protein_length <= 0:
            max_protein_length = 1.0

        gene_title = (
            f"{self.gene_name} - {self.gene_data['specie']}  |  "
            f"chr{self.gene_data['chromosome']}:"
            f"{self._format_number(genomic_start)}-{self._format_number(genomic_end)}  |  "
            f"{num_transcripts} transcript(s)"
        )
        if junction_display_df is not None and 'cluster_name' in junction_display_df.columns:
            cluster_values = junction_display_df['cluster_name'].dropna()
            if len(cluster_values) > 0:
                gene_title += f"  |  {cluster_values.iloc[0]}"

        # Split into pages
        pages = [valid_transcripts[i:i + transcripts_per_page]
                 for i in range(0, num_transcripts, transcripts_per_page)]
        num_pages = len(pages)

        with PdfPages(output_file) as pdf:
            for page_idx, page_transcripts in enumerate(pages):
                n = len(page_transcripts)
                show_junction_table = page_idx == 0 and junction_display_df is not None and len(junction_display_df) > 0

                # Height layout per page:
                #   row 0  : page title
                #   row 1  : optional junction table on first page
                #   next   : axis scales
                #   next   : spacer
                #   rows.. : n * (exon/genomic row + protein/domain row + label row)
                if show_junction_table:
                    height_ratios = [0.55, 1.10, 0.75, 0.28] + [0.7, 1.2, 0.03] * n
                    scale_row = 2
                    transcript_start_row = 4
                else:
                    height_ratios = [0.55, 0.75, 0.28] + [0.7, 1.2, 0.03] * n
                    scale_row = 1
                    transcript_start_row = 3
                total_rows = len(height_ratios)

                fig = plt.figure(figsize=(11, 8.5))
                gs = GridSpec(
                    total_rows, 2,
                    figure=fig,
                    height_ratios=height_ratios,
                    width_ratios=[1.2, 1],
                    hspace=0.16,
                    wspace=0.35,
                    left=0.08,
                    right=0.95,
                    top=0.97,
                    bottom=0.04,
                )

                # ── title row ──────────────────────────────────────────────
                ax_title = fig.add_subplot(gs[0, :])
                ax_title.axis('off')
                page_label = f"  (page {page_idx + 1}/{num_pages})" if num_pages > 1 else ""
                ax_title.text(0.5, 0.5, gene_title + page_label,
                              ha='center', va='center',
                              fontsize=11, fontweight='bold',
                              transform=ax_title.transAxes)

                if show_junction_table:
                    ax_table = fig.add_subplot(gs[1, :])
                    self._draw_junction_table(ax_table, junction_display_df)

                # ── axis scale row (repeated on every page) ────────────────
                self._draw_genomic_scale(fig, gs[scale_row, 0], genomic_start, genomic_end)
                self._draw_protein_scale(fig, gs[scale_row, 1], max_protein_length)

                # ── transcript rows ────────────────────────────────────────
                for i, transcript in enumerate(page_transcripts):
                    row = transcript_start_row + i * 3

                    ax_genomic = fig.add_subplot(gs[row:row + 2, 0])
                    self._draw_genomic_view(ax_genomic, transcript,
                                            genomic_start, genomic_end,
                                            df_junction=junction_display_df)

                    ax_protein = fig.add_subplot(gs[row:row + 2, 1])
                    self._draw_protein_view(ax_protein, transcript, max_protein_length)

                    ax_label = fig.add_subplot(gs[row + 2, :])
                    ax_label.axis('off')
                    transcript_name = transcript['info']['transcript_ensembl_id']
                    protein_name = transcript['info'].get('protein_ensembl_id')
                    if protein_name is None or pd.isna(protein_name) or str(protein_name).strip() == '':
                        protein_name = 'N/A'
                    ax_label.text(
                        0.02, 0.98,
                        f"Transcript: {transcript_name}  |  Protein: {protein_name}",
                        fontsize=8, va='top', transform=ax_label.transAxes,
                    )

                pdf.savefig(fig, bbox_inches='tight', dpi=PDF_RASTER_DPI)
                plt.close(fig)

        print(f"PDF saved to: {output_file}  ({num_pages} page(s))")
    
    def _draw_genomic_scale(self, fig, gs_position, start, end):
        """Draw genomic coordinate scale with alternating above/below labels."""
        ax = fig.add_subplot(gs_position)
        left = min(start, end)
        right = max(start, end)
        if self._is_negative_strand():
            ax.set_xlim(right, left)
        else:
            ax.set_xlim(left, right)
        ax.set_ylim(-0.5, 1.2)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.set_xlabel('Genomic Position (bp)', fontsize=9, labelpad=18)
        ax.grid(axis='x', alpha=0.3, linestyle='--')

        # Compute tick positions then draw labels alternating above/below axis.
        locator = plt.MaxNLocator(nbins=8, prune='both')
        locator.set_axis(ax.xaxis)
        ticks = locator.tick_values(left, right)
        ticks = [t for t in ticks if left <= t <= right]
        ax.set_xticks(ticks)
        ax.set_xticklabels([])  # hide default labels
        ax.tick_params(axis='x', direction='out', length=4)

        for i, t in enumerate(ticks):
            label = f"{int(round(t)):,}"
            if i % 2 == 0:
                # above the axis line
                ax.text(t, 0.15, label, ha='center', va='bottom', fontsize=6,
                        transform=ax.get_xaxis_transform())
            else:
                # below the axis line
                ax.text(t, -0.35, label, ha='center', va='top', fontsize=6,
                        transform=ax.get_xaxis_transform())
        
    def _draw_protein_scale(self, fig, gs_position, max_length):
        """Draw protein coordinate scale in amino acids."""
        ax = fig.add_subplot(gs_position)
        ax.set_xlim(0, max_length)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.xaxis.tick_bottom()
        ax.xaxis.set_label_position('bottom')
        ax.set_xlabel('Amino Acid Position', fontsize=9, labelpad=3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, prune='both'))
        ax.tick_params(axis='x', labelsize=6.5, rotation=40)
    
    def _draw_genomic_view(self, ax, transcript, genomic_start, genomic_end, df_junction=None):
        """Draw genomic view with exons."""
        left = min(genomic_start, genomic_end)
        right = max(genomic_start, genomic_end)
        if self._is_negative_strand():
            ax.set_xlim(right, left)
        else:
            ax.set_xlim(left, right)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])
        
        exon_y = 0.5
        exon_height = 0.3
        non_cds_height = exon_height * 0.5

        cds_start_tx = transcript['info'].get('cds_start') if 'cds_start' in transcript['info'] else None
        cds_end_tx = transcript['info'].get('cds_end') if 'cds_end' in transcript['info'] else None
        has_transcript_cds = (
            cds_start_tx is not None
            and cds_end_tx is not None
            and pd.notna(cds_start_tx)
            and pd.notna(cds_end_tx)
            and int(cds_start_tx) > 0
            and int(cds_end_tx) > 0
        )
        if has_transcript_cds:
            cds_left = min(int(cds_start_tx), int(cds_end_tx))
            cds_right = max(int(cds_start_tx), int(cds_end_tx))
        else:
            cds_left = None
            cds_right = None
        
        # Draw baseline
        ax.plot([left, right], [exon_y, exon_y], 'k-', linewidth=1.5, zorder=1)
        
        # Draw exons
        for exon_idx, (_, exon) in enumerate(transcript['exons'].iterrows()):
            start = min(exon['genomic_start_tx'], exon['genomic_end_tx'])
            end = max(exon['genomic_start_tx'], exon['genomic_end_tx'])
            width = end - start
            exon_has_cds = (
                pd.notna(exon.get('abs_start_CDS'))
                and pd.notna(exon.get('abs_end_CDS'))
                and int(exon['abs_start_CDS']) > 0
                and int(exon['abs_end_CDS']) >= int(exon['abs_start_CDS'])
            )
            
            # Get exon color
            color = self._get_exon_color(exon)
            if color == '#CCCCCC':
                color = DOCHAP_COLORS[exon_idx % len(DOCHAP_COLORS)]
            
            if not has_transcript_cds:
                # No coding region for this transcript: draw exon as reduced-height non-CDS.
                if exon_has_cds:
                    coding_rect = Rectangle(
                        (start, exon_y - exon_height / 2),
                        width,
                        exon_height,
                        facecolor=color,
                        edgecolor='none',
                        linewidth=0,
                        zorder=4,
                    )
                    ax.add_patch(coding_rect)
                    continue
                gray_rect = Rectangle(
                    (start, exon_y - non_cds_height / 2),
                    width,
                    non_cds_height,
                    facecolor=color,
                    edgecolor='none',
                    alpha=0.35,
                    zorder=3,
                )
                ax.add_patch(gray_rect)
                continue

            coding_start = max(start, cds_left)
            coding_end = min(end, cds_right)

            if coding_end <= coding_start:
                # Exon is outside CDS: draw reduced-height non-CDS.
                if exon_has_cds:
                    coding_rect = Rectangle(
                        (start, exon_y - exon_height / 2),
                        width,
                        exon_height,
                        facecolor=color,
                        edgecolor='none',
                        linewidth=0,
                        zorder=4,
                    )
                    ax.add_patch(coding_rect)
                    continue
                gray_rect = Rectangle(
                    (start, exon_y - non_cds_height / 2),
                    width,
                    non_cds_height,
                    facecolor=color,
                    edgecolor='none',
                    alpha=0.35,
                    zorder=3,
                )
                ax.add_patch(gray_rect)
                continue

            # Draw non-CDS subsegments (5' and/or 3' UTR) at half height.
            if start < coding_start:
                left_non_cds = Rectangle(
                    (start, exon_y - non_cds_height / 2),
                    coding_start - start,
                    non_cds_height,
                    facecolor=color,
                    edgecolor='none',
                    alpha=0.35,
                    zorder=3,
                )
                ax.add_patch(left_non_cds)

            if coding_end < end:
                right_non_cds = Rectangle(
                    (coding_end, exon_y - non_cds_height / 2),
                    end - coding_end,
                    non_cds_height,
                    facecolor=color,
                    edgecolor='none',
                    alpha=0.35,
                    zorder=3,
                )
                ax.add_patch(right_non_cds)

            # Draw CDS part at full height with exon color.
            coding_rect = Rectangle(
                (coding_start, exon_y - exon_height / 2),
                coding_end - coding_start,
                exon_height,
                facecolor=color,
                edgecolor='none',
                linewidth=0,
                zorder=4,
            )
            ax.add_patch(coding_rect)

        matched_junctions = self._get_matching_junctions(transcript, df_junction)
        self._draw_genomic_junctions(ax, matched_junctions, exon_y, exon_height)

    def _draw_transcript_view(self, ax, transcript, max_protein_length):
        """Draw the transcript/protein rectangle view above the domain view."""
        ax.set_xlim(0, max_protein_length)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])

        coding_segments = self._get_coding_exon_segments(transcript)
        rect_y = 0.5
        rect_height = 0.38

        for segment in coding_segments:
            rect = Rectangle(
            (segment['start_aa'], rect_y - rect_height / 2),
            segment['end_aa'] - segment['start_aa'],
                rect_height,
                facecolor=segment['color'],
                edgecolor='black',
                linewidth=1.0,
                zorder=2,
            )
            ax.add_patch(rect)
    
    def _draw_protein_view(self, ax, transcript, max_protein_length):
        """Draw protein/domain view on the right (protein above domains)."""
        ax.set_xlim(0, max_protein_length)
        ax.set_ylim(0, 1)
        # Keep labels/borders as vector, but rasterize dense protein/domain fills.
        ax.set_rasterization_zorder(2.5)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_xticks([])

        protein_y = 0.72
        protein_height = 0.18
        domain_y = 0.32
        domain_height = 0.25
        min_domain_width = max(8.0, max_protein_length * 0.02)

        protein_length_aa = transcript['info'].get('protein_length')
        if protein_length_aa is None or pd.isna(protein_length_aa):
            protein_length_aa = max_protein_length
        protein_length_aa = min(max_protein_length, float(protein_length_aa))

        # Draw protein backbone rectangle and color it by contributing coding exons.
        coding_segments = self._get_coding_exon_segments(transcript)

        protein_bg = Rectangle(
            (0, protein_y - protein_height / 2),
            protein_length_aa,
            protein_height,
            facecolor='white',
            edgecolor='black',
            linewidth=1.2,
            zorder=1,
            rasterized=True,
        )
        ax.add_patch(protein_bg)

        for segment in coding_segments:
            seg_start = max(0, float(segment['start_aa']))
            seg_end = min(protein_length_aa, float(segment['end_aa']))
            if seg_end <= seg_start:
                continue
            seg_rect = Rectangle(
                (seg_start, protein_y - protein_height / 2),
                seg_end - seg_start,
                protein_height,
                facecolor=segment['color'],
                edgecolor='none',
                alpha=1.0,
                zorder=2,
                rasterized=True,
            )
            ax.add_patch(seg_rect)

        if len(transcript['domains']) == 0:
            ax.text(0.5, 0.22, 'No domains', transform=ax.transAxes,
                   ha='center', va='center', fontsize=10, style='italic', color='black')
            return
        
        # Sort domains by size (larger first)
        domains_sorted = transcript['domains'].sort_values('AA_end', ascending=False)
        compact_labels = len(domains_sorted) > 1
        
        # Draw each domain as ellipse with gradient
        label_candidates = []
        for domain_idx, (_, domain) in enumerate(domains_sorted.iterrows()):
            domain_start_aa = float(domain['AA_start'])
            domain_end_aa = float(domain['AA_end'])
            domain_width_raw = domain_end_aa - domain_start_aa
            domain_width = max(domain_width_raw, min_domain_width)
            domain_center = (domain_start_aa + domain_end_aa) / 2
            
            domain_name = self._format_domain_label(domain, compact_mode=compact_labels, max_len=DOMAIN_LABEL_MAX_LEN)
            
            overlapping_segments = []
            for segment in coding_segments:
                overlap_start = max(domain_start_aa, segment['start_aa'])
                overlap_end = min(domain_end_aa, segment['end_aa'])
                if overlap_start < overlap_end:
                    overlapping_segments.append({
                        'start_aa': overlap_start,
                        'end_aa': overlap_end,
                        'color': segment['color'],
                    })
            
            if len(overlapping_segments) == 0:
                # No exon overlap: use a stable non-gray palette fallback.
                fallback_color = DOCHAP_COLORS[domain_idx % len(DOCHAP_COLORS)]
                ellipse = Ellipse(
                    (domain_center, domain_y),
                    domain_width,
                    domain_height,
                    facecolor=fallback_color,
                    edgecolor='black',
                    linewidth=1.2,
                    alpha=0.8,
                    zorder=2,
                )
                ax.add_patch(ellipse)
            elif len(overlapping_segments) == 1:
                # Single exon, solid color
                color = overlapping_segments[0]['color']
                ellipse = Ellipse(
                    (domain_center, domain_y),
                    domain_width,
                    domain_height,
                    facecolor=color,
                    edgecolor='black',
                    linewidth=1.2,
                    alpha=0.8,
                    zorder=2,
                )
                ax.add_patch(ellipse)
            else:
                # Multiple exons: color the ellipse in amino-acid strips.
                num_strips = DOMAIN_STRIP_COUNT
                strip_width = domain_width / num_strips
                
                for i in range(num_strips):
                    strip_aa_pos = domain_start_aa + (i + 0.5) * strip_width
                    strip_color = overlapping_segments[0]['color']
                    
                    for segment in overlapping_segments:
                        if segment['start_aa'] <= strip_aa_pos <= segment['end_aa']:
                            strip_color = segment['color']
                            break
                    
                    # Calculate ellipse segment height (varies across width)
                    # Using ellipse equation: y = b * sqrt(1 - (x/a)^2)
                    x_from_center = strip_aa_pos - domain_center
                    a = domain_width / 2  # Semi-major axis
                    b = domain_height / 2  # Semi-minor axis
                    
                    if abs(x_from_center) < a:
                        # Height at this x position
                        relative_height = np.sqrt(1 - (x_from_center / a) ** 2)
                        strip_height = 2 * b * relative_height
                        
                        # Draw vertical strip
                        rect = Rectangle((strip_aa_pos - strip_width/2, domain_y - strip_height/2),
                                       strip_width * 1.1, strip_height,  # Slight overlap to avoid gaps
                                       facecolor=strip_color, edgecolor='none',
                                       alpha=1.0, zorder=2, rasterized=True)
                        ax.add_patch(rect)
                
                # Draw border ellipse on top
                ellipse_border = Ellipse((domain_center, domain_y), domain_width, domain_height,
                                       facecolor='none', edgecolor='black', linewidth=1.2, zorder=3)
                ax.add_patch(ellipse_border)
            
            # Collect candidates; label selection for overlaps is handled after drawing.
            label_candidates.append({
                'start': domain_start_aa,
                'end': domain_end_aa,
                'center': domain_center,
                'width': max(1.0, domain_width_raw),
                'aa_span': max(1.0, domain_end_aa - domain_start_aa + 1),
                'domain': domain,
            })

        label_items = []
        selected_labels = self._select_longest_labels_for_overlaps(label_candidates)
        for selected in selected_labels:
            domain_name = self._format_domain_label(
                selected['domain'],
                compact_mode=False,
                max_len=DOMAIN_LABEL_MAX_LEN,
            )
            if not domain_name:
                continue
            if selected['has_overlap'] and not domain_name.endswith('*'):
                domain_name = f"{domain_name}*"
            label_items.append({
                'center': selected['center'],
                'text': domain_name,
                'width': max(1.0, selected['width']),
            })

        label_base_y = domain_y - (domain_height / 2) - 0.05
        placed_labels = self._compute_domain_label_positions(
            label_items,
            max_protein_length,
            label_base_y,
            lane_step=0.08,
            lanes=4,
        )

        for label in placed_labels:
            ax.plot(
                [label['center'], label['center']],
                [domain_y - domain_height / 2, label['label_y'] + 0.01],
                color='gray',
                linewidth=0.5,
                zorder=4,
                alpha=0.8,
            )
            ax.text(
                label['center'],
                label['label_y'],
                label['text'],
                ha='center',
                va='top',
                fontsize=6.2,
                fontweight='bold',
                style='italic',
                zorder=5,
                rotation=32,
                color='black',
                clip_on=False,
            )


def generate_gene_pdf(gene_name, conn, output_file=None,
                      protein_only=False, domains_only=False):
    """
    Generate a PDF visualization for a gene similar to DoChap-web.
    
    Parameters:
    -----------
    gene_name : str
        Gene symbol (e.g., 'PUF60', 'BRCA1')
    conn : sqlite3.Connection
        Open SQLite connection to the DoChap database
    output_file : str, optional
        Path to save the PDF. If None, uses gene_name.pdf
    protein_only : bool, optional
        If True, include only transcripts that produce protein.
    domains_only : bool, optional
        If True, include only transcripts whose protein has at least one domain.
    
    Returns:
    --------
    str
        Path to the generated PDF file
    
    Example:
    --------
    >>> conn = sqlite3.connect('../DoChaP-web/DB_merged.sqlite')
    >>> generate_gene_pdf('PUF60', conn)
    >>> conn.close()
    """
    if output_file is None:
        output_file = f"{gene_name}_visualization.pdf"

    # Create visualization using provided DB connection
    viz = GeneVisualization(conn, gene_name)
    viz.create_pdf(
        output_file,
        protein_only=protein_only,
        domains_only=domains_only,
    )
    return output_file


if __name__ == "__main__":
    # Example usage
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate gene visualization PDF")
    parser.add_argument("gene_name", type=str, help="Gene symbol (e.g., PUF60)")
    parser.add_argument("-db", "--database", type=str, required=True,
                       help="Path to DoChap database")
    parser.add_argument("-o", "--output", type=str, default=None,
                       help="Output PDF file path")
    parser.add_argument("--protein-only", action="store_true",
                       help="Include only transcripts that produce protein")
    parser.add_argument("--domains-only", action="store_true",
                       help="Include only transcripts whose protein has domains")
    
    args = parser.parse_args()
    
    conn = sqlite3.connect(args.database)
    try:
        output = generate_gene_pdf(
            args.gene_name,
            conn,
            args.output,
            protein_only=args.protein_only,
            domains_only=args.domains_only,
        )
    finally:
        conn.close()
    print(f"Generated: {output}")
