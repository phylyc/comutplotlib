import pandas as pd
import re

from src.gistic import Gistic, join_gistics
from src.maf import MAF, join_mafs
from src.seg import SEG, join_segs
from src.sif import SIF, join_sifs
from src.snv import SNV
from src.cnv import CNV
from src.meta import Meta


class ComutData(object):

    uninteresting_effects = [
        # MAF.utr5, MAF.utr3, MAF.flank5, MAF.flank3,
        MAF.synonymous, MAF.silent, MAF.igr, MAF.intron
    ]

    def __init__(
        self,

        maf_paths: list[str] = (),
        maf_pool_as: dict | None = None,

        seg_paths: list[str] = (),
        gistic_paths: list[str] = (),

        mutsig_paths: list[str] = (),

        model_significances: list[str] = (),
        model_names: list[str] = (),

        sif_paths: list[str] = (),
        meta_data_rows: list[str] = (),
        meta_data_rows_per_sample: list[str] = (),

        by: str = MAF.patient,
        column_order: tuple[str] = None,
        index_order: tuple[str] = None,

        interesting_gene: str = None,
        interesting_gene_comut_percent_threshold: float = None,
        interesting_genes: set = None,
        ground_truth_genes: dict[str, list[str]] = None,
        snv_interesting_genes: set = None,
        cnv_interesting_genes: set = None,
        total_prevalence_threshold: float = None,

        low_amp_threshold: int | float = 1,
        high_amp_threshold: int | float = 2,
        low_del_threshold: int | float = -1,
        high_del_threshold: int | float = -2,
    ):
        self.maf = join_mafs([MAF.from_file(path_to_file=maf) for maf in maf_paths])
        self.maf.pool_annotations(pool_as=maf_pool_as, inplace=True)
        self.maf.select(selection={MAF.effect: self.uninteresting_effects}, complement=True, inplace=True)
        self.snv = None

        self.seg = join_segs([SEG.from_file(path_to_file=seg) for seg in seg_paths])
        self.gistic = join_gistics([Gistic.from_file(path_to_file=gistic) for gistic in gistic_paths])
        self.cnv = None

        self.mutsig = pd.concat(
            [pd.read_csv(path_to_file, index_col=0, sep="\t") for path_to_file in mutsig_paths]
        ) if mutsig_paths is not None else None

        self.model_significance = pd.DataFrame.from_dict(
            {
                name: pd.read_csv(path_to_file, index_col=0, sep="\t")
                for name, path_to_file in zip(model_names, model_significances)
            }
        ) if model_significances is not None else None
        self.model_names = model_names

        self.sif = join_sifs([SIF.from_file(path_to_file=sif) for sif in sif_paths])
        self.sif.add_annotations(inplace=True)
        self.meta = None
        self.meta_data_rows = meta_data_rows
        self.meta_data_rows_per_sample = meta_data_rows_per_sample

        self.tmb = None

        self.interesting_gene = interesting_gene
        self.interesting_gene_comut_percent_threshold = interesting_gene_comut_percent_threshold
        self.snv_interesting_genes = set(snv_interesting_genes) if snv_interesting_genes is not None else set()
        self.cnv_interesting_genes = set(cnv_interesting_genes) if cnv_interesting_genes is not None else set()
        self.interesting_genes = set(interesting_genes) if interesting_genes is not None else set()
        self.ground_truth_genes = ground_truth_genes
        self.total_prevalence_threshold = total_prevalence_threshold
        self.snv_recurrence_threshold = 5

        self.low_amp_threshold = low_amp_threshold
        self.high_amp_threshold = high_amp_threshold
        self.low_del_threshold = low_del_threshold
        self.high_del_threshold = high_del_threshold

        self.by = by
        self.col_order = column_order
        self.idx_order = index_order
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
        self.cnv = CNV(seg=self.seg, gistic=self.gistic, low_amp_threshold=self.low_amp_threshold, high_amp_threshold=self.high_amp_threshold, low_del_threshold=self.low_del_threshold, high_del_threshold=self.high_del_threshold)
        self.meta = Meta(sif=self.sif, by=self.columns.name, rows=self.meta_data_rows, rows_per_sample=self.meta_data_rows_per_sample)
        self.genes = self.get_genes()
        self.tmb = self.get_tmb()
        self.reindex_data()
        self.sort_genes()
        self.reindex_data()
        self.sort_columns()
        self.reindex_data()

    def reindex_data(self):
        if self.genes is not None:
            self.snv.reindex(index=self.genes)
            self.cnv.reindex(index=self.genes)
        if self.columns is not None:
            self.snv.reindex(columns=self.columns)
            self.cnv.reindex(columns=self.columns)
            self.meta.reindex(columns=self.columns)
            if self.mutsig is not None:
                self.mutsig = self.mutsig.reindex(index=self.columns)
            if self.tmb is not None:
                self.tmb = self.tmb.reindex(index=self.columns)

    def get_columns(self):
        if self.col_order is not None:
            return self.col_order
        if not self.sif.empty:
            if self.by == MAF.sample:
                columns = pd.Index(self.sif.samples, name=SIF.sample)
            elif self.by == MAF.patient:
                columns = pd.Index(self.sif.patients, name=SIF.patient)
            else:
                raise ValueError("The argument 'by' needs to be one of {" + f"{MAF.sample}, {MAF.patient}" + "} but received " + f"{self.by}.")
        return columns

    def get_genes(self):
        if self.interesting_gene is not None:
            has_mut = self.snv.has_snv | self.cnv.has_high_cnv
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

        if self.total_prevalence_threshold is not None:
            prevalence = self.get_total_prevalence()
            prevalence = prevalence.loc[prevalence >= self.total_prevalence_threshold]
            self.interesting_genes = {g for g in self.interesting_genes if g in prevalence.index}

        return pd.Index(self.interesting_genes, name=MAF.gene_name)

    def get_total_prevalence(self):
        has_mut = (self.snv.has_snv | self.cnv.has_high_cnv).fillna(False)
        return has_mut.astype(int).sum(axis=1) / len(self.columns)

    def get_total_prevalence_overall(self):
        has_mut = self.snv.has_snv | self.cnv.has_high_cnv
        return has_mut.any(axis=0).astype(int).sum() / len(self.columns)

    def get_tmb(self):
        meta_tmb = self.meta.get_tmb()
        if meta_tmb is not None:
            return meta_tmb
        else:
            tmb = self.snv.get_tmb()
            return tmb.reindex(index=self.columns).fillna(0)

    def sort_genes(self):
        sorted_features = (
            (self.snv.has_snv.astype(int) + self.cnv.has_high_cnv.astype(int))
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

        # group genes in the same cytoband together since they are co-amplified or co-deleted.
        cytobands = self.cnv.gistic.cytoband.loc[genes]
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
                key=lambda cytoband: [c if i % 2 else int(c) for i, c in enumerate(re.split('(\d+)', cytoband)[1:-1])]
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
            has_high_mut = (
                4 * self.snv.has_snv.astype(int)
                + 3 * (self.snv.has_snv & ~self.cnv.has_high_cnv).fillna(False).astype(int)
                + 2 * self.cnv.has_high_amp.astype(int)
                + self.cnv.has_high_del.astype(int)
            )
            has_burden = self.tmb.gt(0).any(axis=1).astype(int).to_frame("burden").T
            has_any_cnv = self.cnv.has_cnv.any(axis=0).astype(int).to_frame("has_cnv").T
            has_low_cnv = (
                + 2 * self.cnv.has_low_amp.astype(int)
                + self.cnv.has_low_del.astype(int)
            )
            columns = (
                pd.concat([has_low_cnv, has_any_cnv, has_burden, has_high_mut])
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
