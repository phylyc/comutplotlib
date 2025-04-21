from functools import reduce
import os
import pandas as pd
import re

from comutplotlib.gistic import Gistic, join_gistics
from comutplotlib.maf import MAF, join_mafs
from comutplotlib.seg import SEG, join_segs
from comutplotlib.sif import SIF, join_sifs
from comutplotlib.snv import SNV
from comutplotlib.cnv import CNV
from comutplotlib.meta import Meta


class ComutData(object):

    uninteresting_effects = [
        # MAF.utr5, MAF.utr3, MAF.flank5, MAF.flank3,
        MAF.synonymous, MAF.silent, MAF.igr, MAF.intron
    ]

    def __init__(
        self,

        maf_paths: list[str] = None,
        maf_pool_as: dict | None = None,

        seg_paths: list[str] = None,
        gistic_paths: list[str] = None,

        mutsig_paths: list[str] = None,

        sif_paths: list[str] = None,
        meta_data_rows: list[str] = (),
        meta_data_rows_per_sample: list[str] = (),

        drop_empty_columns: bool = False,

        by: str = MAF.patient,
        column_order: tuple[str] = None,
        index_order: tuple[str] = None,
        column_sort_by: tuple[str] = None,

        interesting_gene: str = None,
        interesting_gene_comut_percent_threshold: float = None,
        interesting_genes: set = None,
        ground_truth_genes: dict[str, list[str]] = None,
        snv_interesting_genes: set = None,
        cnv_interesting_genes: set = None,
        total_recurrence_threshold: float = None,
        snv_recurrence_threshold: int = 5,

        low_amp_threshold: int | float = 1,
        mid_amp_threshold: int | float = 1.5,
        high_amp_threshold: int | float = 2,
        baseline: int | float = 0,
        low_del_threshold: int | float = -1,
        mid_del_threshold: int | float = -1.5,
        high_del_threshold: int | float = -2,
        show_low_level_cnvs: bool = False,
    ):
        self.maf = join_mafs([MAF.from_file(path_to_file=maf) for maf in maf_paths]) if maf_paths is not None else MAF()
        self.maf.pool_annotations(pool_as=maf_pool_as, inplace=True)
        self.maf.select(selection={MAF.effect: self.uninteresting_effects}, complement=True, inplace=True)

        self.seg = join_segs([SEG.from_file(path_to_file=seg) for seg in seg_paths]) if seg_paths is not None else SEG()
        self.gistic = join_gistics([Gistic.from_file(path_to_file=gistic) for gistic in gistic_paths]) if gistic_paths is not None else Gistic()
        if not show_low_level_cnvs:
            self.gistic.data = (
                self.gistic.data
                .replace([low_del_threshold, low_amp_threshold], baseline)
                .replace([high_amp_threshold], mid_amp_threshold)
                .replace([high_del_threshold], mid_del_threshold)
            )

        self.mutsig = pd.concat([pd.read_csv(path_to_file, index_col=0, sep="\t") for path_to_file in mutsig_paths]) if mutsig_paths is not None else None

        self.sif = join_sifs([SIF.from_file(path_to_file=sif) for sif in sif_paths]) if sif_paths is not None else SIF()
        self.sif.add_annotations(inplace=True)

        self.snv = None
        self.cnv = None
        self.tmb = None
        self.meta = None

        self.col_order = column_order
        self.columns = None

        self.meta_data_rows = meta_data_rows
        self.meta_data_rows_per_sample = meta_data_rows_per_sample

        self.drop_empty_columns = drop_empty_columns

        self.interesting_gene = interesting_gene
        self.interesting_gene_comut_percent_threshold = interesting_gene_comut_percent_threshold
        self.snv_interesting_genes = set(snv_interesting_genes) if snv_interesting_genes is not None else set()
        self.cnv_interesting_genes = set(cnv_interesting_genes) if cnv_interesting_genes is not None else set()
        self.interesting_genes = set(interesting_genes) if interesting_genes is not None else set()
        self.ground_truth_genes = ground_truth_genes
        self.total_recurrence_threshold = total_recurrence_threshold
        self.snv_recurrence_threshold = snv_recurrence_threshold

        self.low_amp_threshold = low_amp_threshold
        self.mid_amp_threshold = mid_amp_threshold
        self.high_amp_threshold = high_amp_threshold
        self.baseline = baseline
        self.low_del_threshold = low_del_threshold
        self.mid_del_threshold = mid_del_threshold
        self.high_del_threshold = high_del_threshold

        self.by = by
        self.col_order = column_order
        self.idx_order = index_order
        self.column_sort_by = column_sort_by
        self.columns = None
        self.genes = None

    def get_dimensions(self):
        n_genes = len(self.genes) if self.genes is not None else 0
        n_samples = len(self.columns) if self.columns is not None else 0
        n_meta = len(self.meta.rows) if self.meta.rows is not None else 0
        return n_genes, n_samples, n_meta

    def preprocess(self):
        self.columns = self.get_columns()
        self.snv = SNV(maf=self.maf, by=self.by)
        self.cnv = CNV(seg=self.seg, gistic=self.gistic, baseline=self.baseline,
                       low_amp_threshold=self.low_amp_threshold, mid_amp_threshold=self.mid_amp_threshold, high_amp_threshold=self.high_amp_threshold,
                       low_del_threshold=self.low_del_threshold, mid_del_threshold=self.mid_del_threshold, high_del_threshold=self.high_del_threshold)
        self.meta = Meta(sif=self.sif, by=self.columns.name, rows=self.meta_data_rows, rows_per_sample=self.meta_data_rows_per_sample)
        self.reindex_data()
        self.genes = self.get_genes()
        self.tmb = self.get_tmb()
        self.reindex_data()
        if self.drop_empty_columns:
            self._drop_empty_columns()
        self.sort_genes()
        self.reindex_data()
        self.sort_columns()
        self.reindex_data()

    def save(self, out_dir, name):
        # Write columns and genes to file, each comma separated
        # This can be read into the bash script using
        # --column_order $(cat <filename>.columns.txt) --index_order $(cat <filename>.genes.txt)
        with open(os.path.join(out_dir, f"{name}.columns.txt"), "w+") as f:
            f.write(",".join(self.columns))
        with open(os.path.join(out_dir, f"{name}.genes.txt"), "w+") as f:
            f.write(",".join(self.genes))

    def _drop_empty_columns(self):
        not_empty = self.snv.has_snv.any(axis=0) | ~self.cnv.isna.all(axis=0)
        self.columns = self.columns[not_empty]
        self.reindex_data()

    def reindex_data(self):
        if self.genes is not None:
            self.snv.reindex(index=self.genes)
            self.cnv.reindex(index=self.genes)
        else:
            genes = self.snv.df.index.union(self.cnv.df.index)
            self.snv.reindex(index=genes)
            self.cnv.reindex(index=genes)

        if self.columns is not None:
            self.snv.reindex(columns=self.columns)
            self.cnv.reindex(columns=self.columns)
            self.meta.reindex(columns=self.columns)
            if self.mutsig is not None:
                self.mutsig = self.mutsig.reindex(index=self.columns)
            if self.tmb is not None:
                self.tmb = self.tmb.reindex(index=self.columns)
        else:
            columns = self.snv.df.columns.union(self.cnv.df.columns).union(self.meta.df.columns)
            self.snv.reindex(columns=columns)
            self.cnv.reindex(columns=columns)
            self.meta.reindex(columns=columns)
            if self.mutsig is not None:
                self.mutsig = self.mutsig.reindex(index=columns)
            if self.tmb is not None:
                self.tmb = self.tmb.reindex(index=columns)

    def get_columns(self):
        if not self.sif.empty:
            entity = self.sif
        elif not self.maf.empty:
            entity = self.maf
        else:
            entity = self.gistic

        if self.by == MAF.sample:
            columns = pd.Index(entity.samples, name=SIF.sample)
        elif self.by == MAF.patient:
            columns = pd.Index(entity.patients, name=SIF.patient)
        else:
            raise ValueError("The argument 'by' needs to be one of {" + f"{MAF.sample}, {MAF.patient}" + "} but received " + f"{self.by}.")
        return columns

    def get_genes(self):
        if self.interesting_gene is not None:
            has_mut = self.snv.has_snv | self.cnv.has_high_cnv | self.cnv.has_mid_cnv
            has_mut_in_gene = has_mut.loc[self.interesting_gene]
            total_mut_in_gene = has_mut_in_gene.sum()
            good_genes = has_mut.apply(
                lambda g: (g & has_mut_in_gene).sum() / total_mut_in_gene >= self.interesting_gene_comut_percent_threshold,
                axis=1
            )
            self.snv_interesting_genes = set(good_genes.loc[good_genes].index)
            self.cnv_interesting_genes = self.snv_interesting_genes
            self.interesting_genes = self.snv_interesting_genes | self.cnv_interesting_genes
            self.snv_recurrence_threshold = int(total_mut_in_gene * self.interesting_gene_comut_percent_threshold) + 1

        self.interesting_genes |= self.snv_interesting_genes | self.cnv_interesting_genes

        # if ground_truth_genes is not None:
        #     for _, ground_truth_gene_list in ground_truth_genes.items():
        #         self.interesting_genes |= set(ground_truth_gene_list)

        if self.total_recurrence_threshold is not None:
            recurrence = self.get_total_recurrence()["high"]
            recurrence = recurrence.loc[recurrence >= self.total_recurrence_threshold]
            self.interesting_genes = {g for g in self.interesting_genes if g in recurrence.index}

        return pd.Index(self.interesting_genes, name=MAF.gene_name)

    def get_total_recurrence(self):
        has_high_mut = (self.snv.has_snv | self.cnv.has_high_cnv | self.cnv.has_mid_cnv).fillna(False)
        has_low_mut = (self.snv.has_snv | self.cnv.has_high_cnv | self.cnv.has_mid_cnv | self.cnv.has_low_cnv).fillna(False)
        has_mut = pd.concat([
            has_high_mut.astype(int).sum(axis=1).to_frame("high"),
            has_low_mut.astype(int).sum(axis=1).to_frame("low"),
        ], axis=1) / len(self.columns)
        return has_mut

    def get_total_recurrence_overall(self):
        has_high_mut = (self.snv.has_snv | self.cnv.has_high_cnv | self.cnv.has_mid_cnv).fillna(False)
        has_low_mut = (self.snv.has_snv | self.cnv.has_high_cnv | self.cnv.has_mid_cnv | self.cnv.has_low_cnv).fillna(False)
        has_mut = {
            "high": has_high_mut.any(axis=0).astype(int).sum(),
            "low": has_low_mut.any(axis=0).astype(int).sum(),
        }
        return has_mut, len(self.columns)

    def get_tmb(self):
        meta_tmb = self.meta.get_tmb()
        if meta_tmb is not None:
            return meta_tmb
        elif not self.snv.maf.empty:
            tmb = self.snv.get_tmb()
            return tmb.reindex(index=self.columns).fillna(0)
        else:
            return None

    def sort_genes(self):
        if self.idx_order is not None:
            self.genes = pd.Index([c for c in self.idx_order if c in self.genes], name=self.genes.name)
        else:
            sorted_features = (
                (self.snv.has_snv.astype(int) + self.cnv.has_high_cnv.astype(int) + self.cnv.has_mid_cnv.astype(int))
                .fillna(0)
                .sum(axis=1)
                .to_frame("mut_count")
                .join(self.snv.has_snv.any(axis=1).to_frame("has_snv"))
                .join(self.cnv.has_low_cnv.any(axis=1).to_frame("has_low_cnv"))
                .join(self.cnv.has_low_cnv.astype(int).sum(axis=1).to_frame("low_cnv_count"))
                .sort_values(by=["mut_count", "has_snv", "has_low_cnv", "low_cnv_count"], ascending=True)
            )
            genes = sorted_features.index.to_list()
            # bring the interesting gene to the front if the heatmap:
            if self.interesting_gene is not None:
                genes.pop(genes.index(self.interesting_gene))
                genes += [self.interesting_gene]

            if not len(genes):
                return None

            # group genes in the same cytoband together since they are co-amplified or co-deleted.
            cytobands = self.cnv.gistic.cytoband.reindex(index=genes).fillna("")
            cytoband_groups = cytobands[::-1].drop_duplicates().values
            cytoband_key = {cb: i for i, cb in enumerate(cytoband_groups)}

            gene_key = {}
            feature_count = 0
            previous_row_tuple = None
            for gene, row in sorted_features.loc[genes[::-1]].iterrows():
                row_tuple = tuple(row)
                if previous_row_tuple is None or row_tuple != previous_row_tuple:
                    feature_count += 1
                gene_key[gene] = feature_count
                previous_row_tuple = row_tuple

            genes = cytobands.reset_index().set_axis(genes, axis=0).apply(
                lambda row: (cytoband_key[row[Gistic._cytoband]], gene_key[row[MAF.gene_name]], row[MAF.gene_name]),
                axis=1
            ).sort_values(ascending=False).index.to_list()

            # bring the interesting gene to the front (again):
            # and remove genes in the same or neighboring cytobands
            if self.interesting_gene is not None:
                # genes.pop(genes.index(interesting_gene))
                sorted_cytoband_groups = sorted(
                    cytoband_groups,
                    key=lambda cytoband: [c if i % 2 else int(c) for i, c in enumerate(re.split(r'(\d+)', cytoband)[1:-1])]
                )
                sorted_cytoband_key = {cb: i for i, cb in enumerate(sorted_cytoband_groups)}
                cytoband_diff = cytobands.apply(lambda c: sorted_cytoband_key[c]) - sorted_cytoband_key[cytobands.loc[self.interesting_gene]]
                close_genes = cytoband_diff.abs() < 5
                for g in close_genes.loc[close_genes].index:
                    genes.pop(genes.index(g))
                genes += [self.interesting_gene]

            self.genes = pd.Index(genes, name=MAF.gene_name)

    def sort_columns(self):
        if self.col_order is not None:
            self.columns = pd.Index([c for c in self.col_order if c in self.columns], name=self.columns.name)
        else:
            # COMUT: ORDER BY:
            # 1. has high amplification
            # 2. has mid-level amplification
            # 3. has high deletion
            # 4. has mid-level deletion
            # 5. has high/mid CNV and no SNV
            # 6. has SNV
            # 7. todo: SNV type
            def get_score(criterion_list):
                log_weights = range(len(criterion_list), 0, -1)
                return reduce(lambda a, b: a + b, [10 ** w * c for w, c in zip(log_weights, criterion_list)])

            has_high_mut = get_score([
                self.cnv.has_high_amp.astype(int),
                self.cnv.has_mid_amp.astype(int),
                self.cnv.has_high_del.astype(int),
                self.cnv.has_mid_del.astype(int),
                (~self.snv.has_snv & (self.cnv.has_high_cnv | self.cnv.has_mid_cnv)).fillna(False).astype(int),
                self.snv.has_snv.astype(int),
            ])
            if self.tmb is not None:
                if SIF.tmb in self.tmb:
                    has_burden = self.tmb[SIF.tmb].gt(0).astype(int).to_frame("has_burden").T
                    burden = self.tmb[[SIF.tmb]].T
                else:
                    has_burden = self.tmb.sum(axis=1).gt(0).astype(int).to_frame("has_burden").T
                    burden = self.tmb.sum(axis=1).to_frame(SIF.tmb).T
            else:
                has_burden = pd.DataFrame()
                burden = pd.DataFrame()
            has_any_cnv = self.cnv.has_cnv.any(axis=0).astype(int).to_frame("has_cnv").T
            has_low_cnv = get_score([
                self.cnv.has_low_amp.astype(int),
                self.cnv.has_low_del.astype(int),
            ])
            comut_features = [df for df in [self.snv.deleteriousness_score, burden, has_low_cnv, has_any_cnv, has_burden, has_high_mut] if not df.empty]

            features = []
            for col in reversed(self.column_sort_by):
                if col == "COMUT":
                    features += comut_features
                elif col in self.meta_data_rows:
                    if col in self.meta_data_rows_per_sample:
                        features.append(self.meta.df[col].apply(lambda l: l[0] if len(l) else 0).to_frame(col).T)
                    else:
                        features.append(self.meta.df[[col]].T)
                else:
                    pass

            if len(features):
                columns = (
                    pd.concat(features)
                    .T
                    .apply(lambda x: tuple(reversed(tuple(x))), axis=1)
                    .sort_values(ascending=False)
                    .index
                )
                self.columns = pd.Index(columns, name=self.columns.name)

    def get_model_annotation(self):
        return pd.DataFrame(
            [[g in self.snv_interesting_genes, g in self.cnv_interesting_genes] for g in self.genes],
            index=self.genes,
            columns=["snv", "cnv"]
        )
