import logging
import math

import pandas as pd

from generate_gene_pdf import GeneVisualization
import domas


logger = logging.getLogger(__name__)
            
class ClusterAnalysisResult:
    def __init__(self, cluster_name, gene_ensembl_id, gene_symbol, chromosome=None, as_event_type=None):
        self.cluster_name = cluster_name
        self.gene_ensembl_id = gene_ensembl_id
        self.gene_symbol = gene_symbol
        self.chromosome = chromosome
        self.as_event_type = as_event_type
        self.canonical_transcript_id = None
        self.junctions = []
        self.events = []
        

    def add_event(self, event, transcript_id=None, domain_name=None, canonical_domain_length=None, transcript_domain_length=None,
                  canonical_domains_number=None, transcript_domains_number=None):
        self.events.append((event, transcript_id, domain_name, canonical_domain_length, 
                            transcript_domain_length, canonical_domains_number, transcript_domains_number))

    def print_results(self, file_name='analysis_results.txt'):
        with open(file_name, 'a') as f:
            f.write(f"Cluster: {self.cluster_name}, Gene: {self.gene_symbol} ({self.gene_ensembl_id}), Chromosome: {self.chromosome}\n")
            f.write(f"\tCanonical Transcript ID: {self.canonical_transcript_id}\n")
            f.write(f"\tJunctions: {self.junctions}\n")
            f.write(f"\tEvents:\n")
            for event, transcript_id, domain_name, canonical_domain_length, transcript_domain_length, \
                canonical_domains_number, transcript_domains_number in self.events:
                msg = f"\t\t{event}: Transcript ID={transcript_id}, Domain Name={domain_name}, "
                msg += f"Canonical Length={canonical_domain_length}, Transcript Length={transcript_domain_length}, "
                msg += f"Canonical Domains Number={canonical_domains_number}, Transcript Domains Number={transcript_domains_number}\n"
                f.write(msg)

    def get_results_df(self):
        return pd.DataFrame(
            self.events,
            columns=['event', 'transcript_id', 'domain_name', 'c_domain_length', 't_domain_length',
                        'c_domains_number', 't_domains_number']
        )


class JunctionsAnalysis:
    def __init__(self, con, logger_instance=None, gene_visualization_cls=GeneVisualization):
        self.con = con
        self.logger = logger_instance or logger
        self.gene_visualization_cls = gene_visualization_cls

    @staticmethod
    def _get_domain_name(row):
        domain_name_columns = ['interpro', 'pfam', 'cdd', 'smart', 'tigr', 'CDD_id']
        for col in domain_name_columns:
            val = row[col]
            if pd.notna(val) and str(val).strip() not in ['', 'nan', 'None']:
                return val
        return None
    def _find_junction_first_last_exons(self, df_transcript_exons, min_bp, max_bp):
        t_last_exon = df_transcript_exons.loc[
            df_transcript_exons[df_transcript_exons['genomic_start_tx'] <= (max_bp + 1)]['genomic_start_tx'].idxmax()
        ]
        t_first_exon = df_transcript_exons.loc[
            df_transcript_exons[df_transcript_exons['genomic_end_tx'] >= (min_bp - 1)]['genomic_end_tx'].idxmin()
        ]
        return t_first_exon, t_last_exon
    
    def _min_skip_nan(self, values):
        return min((v for v in values if not math.isnan(v)), default=None)
    
    def _max_skip_nan(self, values):
        return max((v for v in values if not math.isnan(v)), default=None)
    
    def _get_min_max_aa(self, first_exon, last_exon):
        min_aa = min(first_exon['abs_start_CDS'], last_exon['abs_start_CDS']) // 3
        max_aa = max(last_exon['abs_end_CDS'], first_exon['abs_end_CDS']) // 3
        return min_aa, max_aa

    def _find_min_max_bp_for_domains(self, df_transcript_exons,domains_in_region):
        domains_in_region = domains_in_region[domains_in_region.AA_start !=0].copy()
        min_domain_aa = (domains_in_region['AA_start'].min()) * 3
        max_domain_aa = (domains_in_region['AA_end'].max()) * 3

        relevant_exons = df_transcript_exons[(df_transcript_exons['abs_start_CDS'] <= max_domain_aa) & (df_transcript_exons['abs_end_CDS'] >= max_domain_aa)]
        min_domains_bp = relevant_exons['genomic_start_tx'].min() if not relevant_exons.empty else float('nan')
        max_domains_bp = relevant_exons['genomic_end_tx'].max() if not relevant_exons.empty else float('nan')
        return min_domains_bp, max_domains_bp
    
    def _find_relevant_domains(self, map_t2e, cluster_domains, canonical_transcript_id, transcript_id, canonical_junctions, transcript_junctions, curr_cluster_result):
        # get max/min of the junction coordinates for both the canonical and transcript junctions
        junction_indices = canonical_junctions | transcript_junctions

        min_bp = min(curr_cluster_result.junctions[idx][0] for idx in junction_indices)
        max_bp = max(curr_cluster_result.junctions[idx][1] for idx in junction_indices)
                
        
        t_map_t2e = map_t2e[transcript_id]
        c_map_t2e = map_t2e[canonical_transcript_id]
        # find the exons in the canonical and transcript that are overlap to the max/min coordinates
        t_first_exon, t_last_exon = self._find_junction_first_last_exons(t_map_t2e, min_bp, max_bp)
        c_first_exon, c_last_exon = self._find_junction_first_last_exons(c_map_t2e, min_bp, max_bp)
        
        # find the domains that overlap to the region between the min/max coordinates in both the canonical and transcript
        t_min_aa1, t_max_aa1 = self._get_min_max_aa(t_first_exon, t_last_exon)
        c_min_aa1, c_max_aa1 = self._get_min_max_aa(c_first_exon, c_last_exon)
        

        df_t_domains = cluster_domains[
            (cluster_domains.transcript_ensembl_id == transcript_id)
            | (cluster_domains.transcript_refseq_id == transcript_id)
        ]
        df_c_domains = cluster_domains[
            (cluster_domains.transcript_ensembl_id == canonical_transcript_id)
            | (cluster_domains.transcript_refseq_id == canonical_transcript_id)
        ]
        t_domains_in_region1 = df_t_domains[
            (df_t_domains['AA_end'] >= t_min_aa1) & (df_t_domains['AA_start'] <= t_max_aa1)
        ].copy()
        c_domains_in_region1 = df_c_domains[
            (df_c_domains['AA_end'] >= c_min_aa1) & (df_c_domains['AA_start'] <= c_max_aa1)
        ].copy()
        # find the max bp coordinates for the domains in the canonical and transcript
        t_min_bp, t_max_bp = self._find_min_max_bp_for_domains(t_map_t2e, t_domains_in_region1)
        c_min_bp, c_max_bp = self._find_min_max_bp_for_domains(c_map_t2e, c_domains_in_region1)
        common_min_bp = self._min_skip_nan([min_bp, t_min_bp, c_min_bp])
        common_max_bp = self._max_skip_nan([max_bp, t_max_bp, c_max_bp])
        if common_min_bp is not None and common_max_bp is not None:
            t_first_exon2, t_last_exon2 = self._find_junction_first_last_exons(t_map_t2e, common_min_bp, common_max_bp)
            t_min_aa2, t_max_aa2 = self._get_min_max_aa(t_first_exon2, t_last_exon2)
            t_domains_in_region2 = df_t_domains[(df_t_domains['AA_end'] >= t_min_aa2) & (df_t_domains['AA_start'] <= t_max_aa2)].copy()
            c_first_exon2, c_last_exon2 = self._find_junction_first_last_exons(c_map_t2e, common_min_bp, common_max_bp)
            c_min_aa2, c_max_aa2 = self._get_min_max_aa(c_first_exon2, c_last_exon2)   
            c_domains_in_region2 = df_c_domains[(df_c_domains['AA_end'] >= c_min_aa2) & (df_c_domains['AA_start'] <= c_max_aa2)].copy()
        else:
            t_domains_in_region2 = t_domains_in_region1
            c_domains_in_region2 = c_domains_in_region1
        
        t_domains_in_region2['length'] = t_domains_in_region2['AA_end'] - t_domains_in_region2['AA_start'] + 1
        c_domains_in_region2['length'] = c_domains_in_region2['AA_end'] - c_domains_in_region2['AA_start'] + 1
        return t_domains_in_region2, c_domains_in_region2

    def _total_covered_length(self, df, idxs, start_col, end_col):
        """
        Total length covered by the union of [start, end] intervals (inclusive),
        merging overlapping intervals so overlaps aren't double-counted.
        """
        if not idxs:
            return None
        intervals = sorted((df.loc[i, start_col], df.loc[i, end_col]) for i in idxs)
        total = 0
        cur_start, cur_end = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_end:  # overlapping intervals -> merge
                cur_end = max(cur_end, e)
            else:
                total += cur_end - cur_start + 1
                cur_start, cur_end = s, e
        total += cur_end - cur_start + 1
        return total
    
    def _compare_domains(self, cluster_domains, map_t2e, canonical_transcript_id, transcript_id, canonical_junctions, transcript_junctions, curr_cluster_result):
        
        t_domains_in_region, c_domains_in_region = self._find_relevant_domains(map_t2e, cluster_domains, canonical_transcript_id, transcript_id, canonical_junctions, transcript_junctions, curr_cluster_result)
        name_cols = ['interpro', 'pfam', 'cdd', 'smart', 'tigr', 'CDD_id']
        prefixes_by_priority = ['IPR', 'pfam', 'cd', 'smart', 'tigr', 'CDD']
        start_col = 'AA_start'
        end_col = 'AA_end'
        length_col = 'length'

        """
        Compare domains of a transcript against the canonical transcript.

        Domains are grouped into "identity groups": two domains belong to the
        same group if they share at least one non-empty name (pfam/interpro/
        smart/etc), transitively.

        Returns one row per group, classified as:
        - C=0, T>0  -> 'new domain'
        - C>0, T=0  -> 'dropped domain'
        - C=1, T=1  -> 'same' / 'longer' / 'shorter' (by length)
        - C=1, T>1  -> 'split domain'
        - C>1, T=1  -> 'merged domain'
        - C>1, T>1, C==T -> 'same_domains' / 'longer_domains' / 'shorter_domains' (by total length)
        - C>1, T>1, C<T  -> 'increased_domain_number'
        - C>1, T>1, C>T  -> 'reduced_domain_number'

        @return: DataFrame with one row per group:
                status, group_id, canonical_idx (list), transcript_idx (list),
                canonical_domain_number, transcript_domain_number,
                canonical_total_length, transcript_total_length, names
        """
        canonical_domain_names = {idx: self._name_set(row, name_cols) for idx, row in c_domains_in_region.iterrows()}
        transcript_domain_names = {idx: self._name_set(row, name_cols) for idx, row in t_domains_in_region.iterrows()}

        items = ([(('C', idx), names) for idx, names in canonical_domain_names.items()]
                + [(('T', idx), names) for idx, names in transcript_domain_names.items()])

        groups = self._group_by_shared_names(items)

        results = []
        group_id = 0

        for members in groups:
            group_id += 1
            c_idxs = sorted([idx for kind, idx in members if kind == 'C'],
                            key=lambda i: c_domains_in_region.loc[i, start_col])
            t_idxs = sorted([idx for kind, idx in members if kind == 'T'],
                            key=lambda i: t_domains_in_region.loc[i, start_col])

            canonical_domains_number, transcript_domains_number = len(c_idxs), len(t_idxs)
            canonical_total_domains_length = self._total_covered_length(c_domains_in_region, c_idxs, start_col, end_col)
            transcript_total_domains_length = self._total_covered_length(t_domains_in_region, t_idxs, start_col, end_col)

            names = set()
            for i in c_idxs:
                names |= canonical_domain_names[i]
            for i in t_idxs:
                names |= transcript_domain_names[i]

            if canonical_domains_number == 0:
                status = 'new domain'
            elif transcript_domains_number == 0:
                status = 'dropped domain'
            elif canonical_domains_number == 1 and transcript_domains_number == 1:
                status = self._classify_pair(transcript_total_domains_length, canonical_total_domains_length)
            elif canonical_domains_number == 1 and transcript_domains_number > 1:
                status = 'split domain'
            elif canonical_domains_number > 1 and transcript_domains_number == 1:
                status = 'merged domain'
            else:  # canonical_domains_number > 1 and transcript_domains_number > 1
                if canonical_domains_number == transcript_domains_number:
                    status = self._classify_pair(transcript_total_domains_length, canonical_total_domains_length) + '_domains'
                elif canonical_domains_number < transcript_domains_number:
                    status = 'increased_domain_number'
                else:
                    status = 'reduced_domain_number'
            domain_name = self._choose_name(names, prefixes_by_priority)
            curr_cluster_result.add_event(
                    status,
                    transcript_id=transcript_id,
                    domain_name=domain_name,
                    canonical_domain_length=canonical_total_domains_length,
                    transcript_domain_length=transcript_total_domains_length,
                    canonical_domains_number=canonical_domains_number,
                    transcript_domains_number=transcript_domains_number,
                )
            

    def _choose_name(self, names: set[str], prefixes: list[str]) -> str | None:
        for prefix in prefixes:
            for s in names:
                if s.lower().startswith(prefix.lower()):
                    return s
        return next(iter(names), None)


    def compare_domains1(self, cluster_domains, map_t2e, canonical_transcript_id, transcript_id, canonical_junctions, transcript_junctions, curr_cluster_result):
        t_domains_in_region, c_domains_in_region = self._find_relevant_domains(map_t2e, cluster_domains, canonical_transcript_id, transcript_id, canonical_junctions, transcript_junctions, curr_cluster_result)
        all_domains = pd.concat([t_domains_in_region, c_domains_in_region], ignore_index=True)
        if all_domains.empty: 
            self.logger.warning(
                f"No domains found in the relevant region for transcript {transcript_id} and canonical transcript {canonical_transcript_id}. Skipping domain comparison."
            )
            curr_cluster_result.add_event('no_domains_in_region', transcript_id=transcript_id)
            return

        common_columns = ['CDD_id', 'cdd', 'pfam', 'smart', 'tigr', 'interpro', 'length']
        all_domains['is_shared'] = all_domains.duplicated(subset=common_columns, keep=False)
        shared_domains = all_domains[all_domains['is_shared']].drop_duplicates(subset=common_columns)

        for _, row in shared_domains.iterrows():
            curr_cluster_result.add_event(
                'same_domain',
                transcript_id=transcript_id,
                domain_name=self._get_domain_name(row),
                canonical_domain_length=row['length'],
                transcript_domain_length=row['length'],
            )

        unique_t_domains = all_domains[~all_domains['is_shared']]
        canonical_domain_names = unique_t_domains[unique_t_domains['canonical'] == '1'].apply(self._get_domain_name, axis=1).tolist()
        transcript_domain_names = unique_t_domains[unique_t_domains['canonical'] == '0'].apply(self._get_domain_name, axis=1).tolist()

        for canonical_domain in canonical_domain_names:
            if canonical_domain not in transcript_domain_names:
                c_row = unique_t_domains[
                    (unique_t_domains['canonical'] != 0) # AMC
                    & (unique_t_domains.apply(lambda row: self._get_domain_name(row) == canonical_domain, axis=1))
                ].iloc[0]
                curr_cluster_result.add_event(
                    'dropped_domain',
                    transcript_id=transcript_id,
                    domain_name=canonical_domain,
                    canonical_domain_length=c_row['length'],
                )
            else:
                c_row = unique_t_domains[
                    (unique_t_domains['canonical'] != 0) # AMC
                    & (unique_t_domains.apply(lambda row: self._get_domain_name(row) == canonical_domain, axis=1))
                ].iloc[0]
                t_row = unique_t_domains[
                    (unique_t_domains['canonical'] == 0) # AMC
                    & (unique_t_domains.apply(lambda row: self._get_domain_name(row) == canonical_domain, axis=1))
                ].iloc[0]
                length_diff = t_row['length'] - c_row['length']
                if length_diff < 0:
                    curr_cluster_result.add_event(
                        'shorter_domain',
                        transcript_id=transcript_id,
                        domain_name=canonical_domain,
                        canonical_domain_length=c_row['length'],
                        transcript_domain_length=t_row['length'],
                    )
                elif length_diff > 0:
                    curr_cluster_result.add_event(
                        'longer_domain',
                        transcript_id=transcript_id,
                        domain_name=canonical_domain,
                        canonical_domain_length=c_row['length'],
                        transcript_domain_length=t_row['length'],
                    )
                else:
                    self.logger.error(
                        f"Unexpected case where domain lengths are the same but domain is not marked as shared for transcript {transcript_id} and canonical domain {canonical_domain}. This may indicate an issue with the domain comparison logic."
                    )
                    

        for transcript_domain in transcript_domain_names:
            if transcript_domain not in canonical_domain_names:
                row = unique_t_domains[
                    (unique_t_domains['canonical'] == 0) # AMC
                    & (unique_t_domains.apply(lambda row: self._get_domain_name(row) == transcript_domain, axis=1))
                ].iloc[0]
                length = row['length']
                curr_cluster_result.add_event(
                    'novel_transcript_domain',
                    transcript_id=transcript_id,
                    domain_name=transcript_domain,
                    canonical_domain_length=0,
                    transcript_domain_length=length,
                )

    def analyze_junctions(self, junctions_csv='as_events_junctions.csv', output_path='as_events_junctions_analysis.csv', df_junctions=None, 
                          hadas_format=False, n=0):
        debug_gene_names = ['PUF60']
        results = []
        if df_junctions is None and junctions_csv is None:
            raise ValueError("Either df_junctions or junctions_csv must be provided.")
        elif df_junctions is not None and junctions_csv is not None:
            raise ValueError("Only one of df_junctions or junctions_csv should be provided.")
        elif junctions_csv is not None:
            if hadas_format:
                df_junctions = domas.read_hadas_input_file(junctions_csv)
            else:
                df_junctions = pd.read_csv(junctions_csv)
        else:
            df_junctions = df_junctions.copy()
        if 'cluster_name' not in df_junctions.columns and 'cluster' in df_junctions.columns:
            df_junctions = df_junctions.rename(columns={'cluster': 'cluster_name'})

        df_junctions = df_junctions[df_junctions['gene_ensembl_id'].str.startswith('ENSG', na=False)]
        # check all needed columns are in df_junctions
        required_columns = ['gene_ensembl_id', 'start_position', 'end_position', 'cluster_name']
        for col in required_columns:
            if col not in df_junctions.columns:
                raise ValueError(f"Column '{col}' is required in df_junctions but not found.")
        #df_junctions = df_junctions[df_junctions['gene_symbol'].isin(debug_gene_names)]
        df_t = pd.read_sql_query('select * from transcripts', self.con)
        if n > 0:
            gene_count = df_t.value_counts('gene_GeneID_id')
            df_t['gene_count'] = df_t['gene_GeneID_id'].map(gene_count)
            ids_n = df_t[df_t['gene_count'] == n].gene_GeneID_id.unique().tolist()
            df_junctions_n = df_junctions[df_junctions['gene_ensembl_id'].isin(ids_n)]
        else:
            df_junctions_n = df_junctions

        gene_ids = df_junctions_n.gene_ensembl_id.unique().tolist()
        transcript_ids = domas.get_gene_transcript_ids(self.con, gene_ids)
        df_transcripts = domas.get_genes_df_transcripts(self.con, gene_ids)
        canonical_transcripts_ids = set(
            df_transcripts[df_transcripts.canonical != 0]
            .transcript_ensembl_id.fillna(df_transcripts.transcript_refseq_id)
            .values.tolist()
        )
        transcript_ids_combined = df_transcripts.transcript_ensembl_id.fillna(df_transcripts.transcript_refseq_id).values.tolist()
        transcript_ids = set(transcript_ids_combined) - set(['', None, 'nan', 'None'])
        df_domains = domas.get_transcript_domains_db(self.con, transcript_ids)

        # Local import avoids circular imports with aleternative_splicing module.
        from aleternative_splicing import get_exons_for_transcripts

        df_exons = get_exons_for_transcripts(self.con, transcript_ids)
        cluster_groups = df_junctions_n.groupby('cluster_name')
        total = len(cluster_groups)
        count = 0

        for cluster_name, cluster_df in cluster_groups:
            count += 1
            if count % 100 == 0:
                self.logger.info(f"Analyzing cluster {count}/{total}: {cluster_name}")

            gene_ensembl_id = cluster_df.gene_ensembl_id.values[0]
            gene_symbol = cluster_df.gene_symbol.values[0]
            event_type = cluster_df.event_type.values[0] if 'event_type' in cluster_df.columns else None
            curr_cluster_result = ClusterAnalysisResult(cluster_name, gene_ensembl_id, gene_symbol, as_event_type=event_type)
            curr_cluster_result.junctions = list(zip(cluster_df['start_position'], cluster_df['end_position']))
            results.append(curr_cluster_result)

            df_gene_transcript = df_transcripts[df_transcripts.gene_ensembl_id == gene_ensembl_id]
            gene_transcript_ids = set(
                df_gene_transcript.transcript_ensembl_id.fillna(df_gene_transcript.transcript_refseq_id).values.tolist()
            )
            gene_canoincal_set = canonical_transcripts_ids.intersection(gene_transcript_ids)
            if not gene_canoincal_set:
                curr_cluster_result.add_event('no_canonical_transcript')
                self.logger.warning(f"No canonical transcript found for cluster {cluster_name}. Skipping analysis.")
                continue

            if len(gene_transcript_ids) == 1:
                curr_cluster_result.add_event('only_one_transcript')
                self.logger.info(f"Only one transcript found for cluster {cluster_name}.")
                continue
            curr_cluster_result.canonical_transcript_id, = gene_canoincal_set
            cluster_domains = df_domains[
                df_domains.transcript_ensembl_id.isin(gene_transcript_ids)
                | df_domains.transcript_refseq_id.isin(gene_transcript_ids)
            ]

            map_t2j = {transcript_id: [] for transcript_id in gene_transcript_ids}
            map_j2t = {i: [] for i in range(len(curr_cluster_result.junctions))}
            map_t2e = {transcript_id: [] for transcript_id in gene_transcript_ids}
            for transcript_id in gene_transcript_ids:
                df_transcript_exons = df_exons[
                    (df_exons.transcript_ensembl_id == transcript_id)
                    | (df_exons.transcript_refseq_id == transcript_id)
                ]
                map_t2e[transcript_id] = df_transcript_exons
                for idx, (start_position, end_position) in enumerate(curr_cluster_result.junctions):
                    start_exon = df_transcript_exons[abs(df_transcript_exons.genomic_start_tx - end_position) <= 1]
                    end_exon = df_transcript_exons[abs(df_transcript_exons.genomic_end_tx - start_position) <= 1]
                    if (
                        (len(start_exon) == 1)
                        and (len(end_exon) == 1)
                        and (abs(start_exon.order_in_transcript.values[0] - end_exon.order_in_transcript.values[0]) == 1)
                    ):
                        map_t2j[transcript_id].append(idx)
                        map_j2t[idx].append(transcript_id)

            for idx, transcripts in map_j2t.items():
                if not transcripts:
                    self.logger.warning(
                        f"Junction {curr_cluster_result.junctions[idx]} in cluster {cluster_name} does not map to any transcript. This may indicate an issue with the junction definition or exon mapping."
                    )
                    curr_cluster_result.add_event('junction_not_mapped', idx)

            canonical_junctions = set(map_t2j.get(curr_cluster_result.canonical_transcript_id, []))
            if not canonical_junctions:
                curr_cluster_result.add_event('no_canonical_junctions')
                self.logger.warning(f"No canonical junctions found for cluster {cluster_name}. Skipping analysis.")
                continue

            for transcript_id, transcript_junctions in map_t2j.items():
                if transcript_id == curr_cluster_result.canonical_transcript_id:
                    continue
                if not transcript_junctions:
                    self.logger.warning(
                        f"Transcript {transcript_id} in cluster {cluster_name} does not have any junctions. This may indicate an issue with the exon mapping for the canonical transcript."
                    )
                    curr_cluster_result.add_event('transcript_doesnt_have_junctions', transcript_id=transcript_id)
                    continue
                transcript_unique_junctions = set(transcript_junctions) - canonical_junctions
                if not transcript_unique_junctions:
                    self.logger.info(
                        f"Transcript {transcript_id} in cluster {cluster_name} does not have any unique junctions compared to the canonical transcript. Skipping this transcript for comparison."
                    )
                    curr_cluster_result.add_event('no_unique_junctions', transcript_id=transcript_id)
                    continue
                self._compare_domains(
                    cluster_domains,
                    map_t2e,
                    curr_cluster_result.canonical_transcript_id,
                    transcript_id,
                    canonical_junctions,
                    set(transcript_junctions),
                    curr_cluster_result,
                )

        count = 0
        for cluster_result in results:
            count += 1
            #cluster_result.print_results(file_name=output_path)
            try:
                df_cluster_junctions = pd.DataFrame(cluster_result.junctions, columns=['start_position', 'end_position'])
                viz = self.gene_visualization_cls(self.con, cluster_result.gene_symbol)
                if cluster_result.as_event_type is None:
                    file_name = f'{cluster_result.gene_symbol}_{count}_junction_comparison.pdf'
                else:
                    file_name = f'{cluster_result.as_event_type}_{cluster_result.gene_symbol}_{count}_junction_comparison.pdf'
                    
                viz.create_pdf(
                    file_name,
                    protein_only=False,
                    domains_only=False,
                    df_junction=df_cluster_junctions,
                    df_results=cluster_result.get_results_df(),
                )
            except ValueError as e:
                print(f"Warning: Skipping PDF generation for {cluster_result.gene_symbol}: {e}")
        # write results to csv
        df_results_columns = ['gene_symbol', 'event_type', 'canonical_transcript_id', 'transcript_id', 'domain_name', 
                              'c_domain_length', 't_domain_length', 'c_domains_number', 't_domains_number']
        all_results = []
        for cluster_result in results:  
            for event, transcript_id, domain_name, canonical_domain_length, transcript_domain_length, \
                canonical_domains_number, transcript_domains_number in cluster_result.events:
                all_results.append([
                    cluster_result.gene_symbol,
                    event,
                    cluster_result.canonical_transcript_id,
                    transcript_id,
                    domain_name,
                    canonical_domain_length,
                    transcript_domain_length,
                    canonical_domains_number,
                    transcript_domains_number
                ])
        df_all_results = pd.DataFrame(all_results, columns=df_results_columns)
        df_all_results.to_csv(output_path, index=False)
        return results




    def _name_set(self, row, name_cols):
        """Return the set of non-empty/non-null identifier names for a domain row."""
        names = set()
        for col in name_cols:
            val = row[col]
            if val is None:
                continue
            if isinstance(val, float) and np.isnan(val):
                continue
            if str(val).strip() in ('', 'None', 'nan'):
                continue
            names.add(val)
        return names


    def _group_by_shared_names(self, items_with_names):
        """
        items_with_names: list of (key, set_of_names)
        Groups keys together if their name-sets overlap (transitively).
        Returns: list of lists of keys.
        """
        groups = []  # list of {'keys': [...], 'names': set(...)}
        for key, names in items_with_names:
            matching = [g for g in groups if g['names'] & names] # match even if one name is identicl
            if not matching:
                groups.append({'keys': [key], 'names': set(names)})
            else:
                merged = matching[0]
                for other in matching[1:]:
                    merged['keys'].extend(other['keys'])
                    merged['names'] |= other['names']
                    groups.remove(other)
                merged['keys'].append(key)
                merged['names'] |= names
        return [g['keys'] for g in groups]


    def _classify_pair(self, t_len, c_len):
        if t_len == c_len:
            return 'same'
        return 'longer' if t_len > c_len else 'shorter'

