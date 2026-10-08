from collections import defaultdict
from collections.abc import Mapping, Sequence
from functools import reduce
import os
import numpy as np
import pandas as pd
import re
from scipy.cluster.hierarchy import linkage, leaves_list, optimal_leaf_ordering
from scipy.spatial.distance import squareform

from comutplotlib.gistic import Gistic, join_gistics
from comutplotlib.maf import MAF, join_mafs
from comutplotlib.mark import Mark, join_marks
from comutplotlib.seg import SEG, join_segs
from comutplotlib.sif import SIF, join_sifs
from comutplotlib.snv import SNV
from comutplotlib.cnv import CNV
from comutplotlib.epi import EPI
from comutplotlib.mutational_signature_set import MutationalSignatureSet
from comutplotlib.meta import Meta
from comutplotlib.grid import (
    GridPartition,
    build_column_groups,
    build_gene_groups,
    trivial_partition,
)


class ComutData(object):

    uninteresting_effects = [
        # MAF.utr5, MAF.utr3, MAF.flank5, MAF.flank3,
        MAF.synonymous, MAF.silent, MAF.igr, MAF.intron
    ]

    def __init__(
        self,
        cohort_name: str | None = None,

        maf_paths: list[str] | None = None,
        maf_pool_as: dict | None = None,

        seg_paths: list[str] | None = None,
        gistic_paths: list[str] | None = None,

        mark_paths: list[str] | None = None,

        signatures_paths: list[str] | None = None,
        group_signatures_by_etiology: bool = False,

        sif_paths: list[str] | None = None,
        meta_data_rows: Sequence[str] = (),
        meta_data_rows_per_sample: Sequence[str] = (),

        drop_empty_columns: bool = False,

        by: str = MAF.patient,
        column_order: tuple[str] | None = None,
        index_order: tuple[str] | None = None,
        column_sort_by: tuple[str] | None = None,
        sort_method: str | None = None,

        interesting_gene: str | None = None,
        interesting_gene_comut_percent_threshold: float | None = None,
        interesting_genes: set | None = None,
        ground_truth_genes: dict[str, list[str]] | None = None,
        snv_interesting_genes: set | None = None,
        cnv_interesting_genes: set | None = None,
        total_recurrence_threshold: float | None = None,
        snv_recurrence_threshold: int = 5,

        low_amp_threshold: int | float = 1,
        mid_amp_threshold: int | float = 1.5,
        high_amp_threshold: int | float = 2,
        baseline: int | float = 0,
        low_del_threshold: int | float = -1,
        mid_del_threshold: int | float = -1.5,
        high_del_threshold: int | float = -2,
        show_low_level_cnvs: bool = False,

        column_group_by: Sequence[str] | None = None,
        gene_groups: dict[str, list[str]] | None = None,
        group_order: Sequence[str] | None = None,
        column_group_labels: Sequence[str] | Mapping[str, str] | None = None,
        drop_ungrouped_genes: bool = False,
        na_group_label: str = "NA",
        other_gene_group_label: str = "Other",
    ):
        self.name = cohort_name
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

        self.mark = join_marks([Mark.from_file(path_to_file=mark) for mark in mark_paths]) if mark_paths is not None else Mark()

        self.signatures = pd.concat([
            pd.read_csv(path_to_file, index_col=0, sep="," if path_to_file.endswith(".csv") else "\t").fillna(0)
            for path_to_file in signatures_paths]
        ) if signatures_paths is not None else None
        self.group_signatures_by_etiology = group_signatures_by_etiology

        self.sif = join_sifs([SIF.from_file(path_to_file=sif) for sif in sif_paths]) if sif_paths is not None else SIF()
        self.sif.add_annotations(inplace=True)

        self.by = by

        self.snv = None
        self.cnv = None
        self.epi = None
        self.tmb = None
        self.meta = None

        self.meta_data_rows = meta_data_rows
        self.meta_data_rows_per_sample = meta_data_rows_per_sample

        self.drop_empty_columns = drop_empty_columns

        self.interesting_gene = interesting_gene
        self.interesting_gene_comut_percent_threshold = interesting_gene_comut_percent_threshold
        self.snv_interesting_genes = set(snv_interesting_genes) if snv_interesting_genes is not None else set()
        self.cnv_interesting_genes = set(cnv_interesting_genes) if cnv_interesting_genes is not None else set()
        self.interesting_genes = set(interesting_genes) if interesting_genes is not None else None
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

        self.col_order = column_order
        self.idx_order = index_order
        self.column_sort_by = column_sort_by if column_sort_by is not None else ()
        self.columns = None
        self.genes = None
        self.gene_sort_method = sort_method
        self.column_sort_method = sort_method
        self.cluster_cnv_weight = 1
        self.cluster_snv_weight = 0.75

        # Grid sub-panel configuration (see GRID_SUBPANELS_PLAN.md). ``self.grid``
        # is None until ``build_grid()`` is called; a trivial (1x1) partition
        # reproduces the legacy single-panel behaviour.
        self.column_group_by = list(column_group_by) if column_group_by is not None else []
        self.gene_groups_config = gene_groups
        self.group_order = list(group_order) if group_order is not None else None
        # Display labels for the grid columns: either a positional list (in
        # grid-column order) or a ``group key -> label`` mapping (used to share
        # the case cohort's labels with the control cohort).
        self.column_group_labels = (
            dict(column_group_labels) if isinstance(column_group_labels, Mapping)
            else list(column_group_labels) if column_group_labels is not None
            else None
        )
        self.drop_ungrouped_genes = drop_ungrouped_genes
        self.na_group_label = na_group_label
        self.other_gene_group_label = other_gene_group_label
        self.grid: GridPartition | None = None

    def get_dimensions(self) -> tuple[int, int, int]:
        n_genes = len(self.genes) if self.genes is not None else 0
        n_samples = len(self.columns) if self.columns is not None else 0
        n_meta = len(self.meta.rows) if self.meta.rows is not None else 0
        return n_genes, n_samples, n_meta

    @property
    def empty(self) -> bool:
        n_genes, n_samples, n_meta = self.get_dimensions()
        return n_samples == 0

    def preprocess(self):
        self.align_signatures()
        self.columns = self.get_columns()
        self.snv = SNV(maf=self.maf, by=self.by)
        self.cnv = CNV(seg=self.seg, gistic=self.gistic, baseline=self.baseline,
                       low_amp_threshold=self.low_amp_threshold, mid_amp_threshold=self.mid_amp_threshold, high_amp_threshold=self.high_amp_threshold,
                       low_del_threshold=self.low_del_threshold, mid_del_threshold=self.mid_del_threshold, high_del_threshold=self.high_del_threshold)
        self.epi = EPI(mark=self.mark)
        self.meta = Meta(sif=self.sif, by=self.columns.name, rows=self.meta_data_rows, rows_per_sample=self.meta_data_rows_per_sample)
        self.reindex_data()
        self.genes = self.get_genes()
        self.tmb = self.get_tmb()
        self.resort_data()

    def align_signatures(self):
        if self.signatures is None:
            return None
        self.signatures.columns = pd.Index([c.replace("Signature_", "SBS") for c in self.signatures.columns])
        if self.group_signatures_by_etiology:
            # Collapse individual signatures into their etiology (the
            # ``MutationalSignatureSet.signature_sets`` keys). The result is
            # already ordered by etiology, so no further column sorting is
            # needed. Grouping is idempotent, which matters because the joint
            # cohort re-runs ``preprocess`` on already-grouped exposures.
            self.signatures = MutationalSignatureSet.group_by_etiology(self.signatures)
        else:
            sorted_columns = MutationalSignatureSet.sort_signatures(signatures=self.signatures.columns)
            self.signatures = self.signatures.reindex(columns=sorted_columns)
        if self.sif is None or self.by != MAF.patient:
            return None
        if self.signatures.index.name == SIF.patient:
            return None
        sample_to_patient_map = self.sif.data[[SIF.sample, SIF.patient]].set_index(SIF.sample)[SIF.patient].to_dict()
        # Map sample-indexed signatures to patients and aggregate. Entries that
        # are not known sample ids are assumed to already be patient ids and keep
        # their own id (``.get(s, s)``); using ``.get(s)`` would map them to
        # ``None``, and the subsequent ``groupby`` would silently drop every row
        # (NaN group keys), wiping out signatures that were provided per patient.
        self.signatures[SIF.patient] = self.signatures.index.map(lambda s: sample_to_patient_map.get(s, s))
        self.signatures = self.signatures.groupby(SIF.patient).agg("sum")

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

    def resort_data(self):
        self.reindex_data()
        if self.drop_empty_columns:
            self._drop_empty_columns()
        self.sort_genes()
        self.reindex_data()
        self.sort_columns()
        self.reindex_data()

    # --- Grid sub-panels -----------------------------------------------------

    @property
    def grouping_enabled(self) -> bool:
        return bool(self.column_group_by) or bool(self.gene_groups_config)

    def _meta_key_frame(self):
        """DataFrame indexed by columns exposing the ``column_group_by`` keys.

        ``Meta.df`` is indexed by samples/patients with meta rows as columns, so a
        column subset yields exactly the per-sample stratification keys.
        """
        if self.meta is None or self.meta.df is None or self.meta.df.empty:
            return None
        keys = [k for k in self.column_group_by if k in self.meta.df.columns]
        if not keys:
            return None
        return self.meta.df[keys]

    def build_grid(self, gene_groups=None):
        """Populate ``self.grid`` from config (or adopt ``gene_groups``).

        The returned partition uses *config-first* member order (unsorted); call
        :meth:`apply_grid_sort` to sort within groups and reindex the carriers.
        A control cohort should pass the case's (already sorted) ``gene_groups``
        so both cohorts share an identical gene (row) axis.
        """
        if not self.grouping_enabled and gene_groups is None:
            self.grid = trivial_partition(self.genes, self.columns)
            return self.grid

        column_groups = build_column_groups(
            columns=self.columns,
            key_frame=self._meta_key_frame(),
            keys=self.column_group_by,
            order=self.group_order,
            na_label=self.na_group_label,
            labels=self.column_group_labels,
        )
        if gene_groups is not None:
            gene_groups = tuple(
                g.with_members([m for m in g.members if m in set(self.genes)])
                for g in gene_groups
            )
        else:
            gene_groups = build_gene_groups(
                genes=self.genes,
                gene_sets=self.gene_groups_config,
                order=self.group_order,
                other_label=self.other_gene_group_label,
                drop_ungrouped=self.drop_ungrouped_genes,
            )
        self.grid = GridPartition(gene_groups=gene_groups, column_groups=column_groups)
        return self.grid

    def column_group_label_map(self) -> dict[str, str]:
        """Resolved ``group key -> label`` map of the current grid columns.

        Used to hand the case cohort's (possibly user-supplied) column-group
        labels to the control cohort, which builds its own column groups and may
        resolve them in a different order.
        """
        if self.grid is None:
            return {}
        return {g.key: g.label for g in self.grid.column_groups}

    def apply_grid_sort(self, gene_groups=None):
        """Build/adopt the grid, sort within each group, and reindex.

        Evidence is *config-first* (Decision 1): gene order within each gene group
        is derived from the first (config-order) column group; column order within
        each column group is derived from the first (config-order) gene group.

        Returns the finalized (sorted) gene groups so a control cohort can adopt
        the exact same gene axis.
        """
        self.build_grid(gene_groups=gene_groups)
        if self.grid is None or self.grid.is_trivial():
            return self.grid.gene_groups if self.grid is not None else ()

        # config-first evidence (membership fixed; only within-group order changes)
        ref_columns = list(self.grid.reference_column_group.members)
        ref_genes = list(self.grid.reference_gene_group.members)

        if gene_groups is None:
            # case: sort genes within each group using the reference column group
            sorted_gene_groups = tuple(
                g.with_members(self._ordered_or_original(self._compute_gene_order(genes=g.members, columns=ref_columns), g.members))
                for g in self.grid.gene_groups
            )
        else:
            # control: adopt the case gene ordering verbatim
            sorted_gene_groups = self.grid.gene_groups

        sorted_column_groups = tuple(
            g.with_members(self._ordered_or_original(self._compute_column_order(columns=g.members, genes=ref_genes), g.members))
            for g in self.grid.column_groups
        )

        self.grid = GridPartition(gene_groups=sorted_gene_groups, column_groups=sorted_column_groups)
        self.apply_grid()
        return self.grid.gene_groups

    @staticmethod
    def _ordered_or_original(order, members):
        """Return ``order`` as a list, falling back to ``members`` when empty/None."""
        if order is None or len(order) == 0:
            return list(members)
        return list(order)

    def apply_grid(self):
        """Reassemble ``self.genes``/``self.columns`` from ``self.grid`` and reindex."""
        if self.grid is None:
            return
        gene_name = self.genes.name if self.genes is not None else MAF.gene_name
        column_name = self.columns.name if self.columns is not None else SIF.sample
        self.genes = self.grid.all_genes(name=gene_name)
        self.columns = self.grid.all_columns(name=column_name)
        self.reindex_data()

    def invert_grid_columns(self):
        """Reverse the grid column order to reproduce the inward-out ordering.

        Mirrors the legacy case/control behaviour (``self.columns[::-1]``) for grid
        layouts: the order of the grid column groups is reversed *and* the panel
        columns within each group are reversed, so a case cohort drawn to the left
        of a control cohort reads inside-out (most-mutated columns towards the
        centre). No-op when there is no grid.
        """
        if self.grid is None:
            return
        inverted_column_groups = tuple(
            g.with_members(tuple(reversed(g.members)))
            for g in reversed(self.grid.column_groups)
        )
        self.grid = self.grid.with_column_groups(inverted_column_groups)
        self.apply_grid()

    def reindex_data(self):
        if self.genes is not None:
            self.snv.reindex(index=self.genes)
            self.cnv.reindex(index=self.genes)
            self.epi.reindex(index=self.genes)
        else:
            genes = self.snv.df.index.union(self.cnv.df.index)
            self.snv.reindex(index=genes)
            self.cnv.reindex(index=genes)
            self.epi.reindex(index=genes)

        if self.columns is not None:
            self.snv.reindex(columns=self.columns)
            self.cnv.reindex(columns=self.columns)
            self.epi.reindex(columns=self.columns)
            self.meta.reindex(columns=self.columns)
            if self.signatures is not None:
                self.signatures = self.signatures.reindex(index=self.columns)
            if self.tmb is not None:
                self.tmb = self.tmb.reindex(index=self.columns)
        else:
            columns = self.snv.df.columns.union(self.cnv.df.columns).union(self.meta.df.columns)
            self.snv.reindex(columns=columns)
            self.cnv.reindex(columns=columns)
            self.epi.reindex(columns=columns)
            self.meta.reindex(columns=columns)
            if self.signatures is not None:
                self.signatures = self.signatures.reindex(index=columns)
            if self.tmb is not None:
                self.tmb = self.tmb.reindex(index=columns)

    def get_columns(self):
        if not self.sif.empty:
            entity = self.sif
        elif not self.maf.empty:
            entity = self.maf
        elif not self.gistic.data.empty:
            entity = self.gistic
        else:
            entity = self.mark

        if self.by == MAF.sample:
            columns = pd.Index(entity.samples, name=SIF.sample)
        elif self.by == MAF.patient:
            columns = pd.Index(entity.patients, name=SIF.patient)
        else:
            raise ValueError("The argument 'by' needs to be one of {" + f"{MAF.sample}, {MAF.patient}" + "} but received " + f"{self.by}.")
        return columns

    def get_genes(self, categories=None):
        if self.genes is not None:
            return self.genes

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

        if self.interesting_genes is None:
            self.interesting_genes = self.snv_interesting_genes | self.cnv_interesting_genes

        # if ground_truth_genes is not None:
        #     for _, ground_truth_gene_list in ground_truth_genes.items():
        #         self.interesting_genes |= set(ground_truth_gene_list)

        if self.total_recurrence_threshold is not None:
            recurrence = self.get_total_recurrence(categories=categories)["high"]
            recurrence = recurrence.loc[recurrence >= self.total_recurrence_threshold]
            self.interesting_genes = {g for g in self.interesting_genes if g in recurrence.index}

        return pd.Index(self.interesting_genes, name=MAF.gene_name)

    def get_mutation_status(self, categories) -> dict[str, pd.DataFrame]:
        """
            Arguments:
                categories dict[str, list[str]]: genes to list of categories {'snv', 'amp', 'del'}
            Output:
                dict[str, pd.DataFrame]: {'low'/'high': data frame of genes by samples/patients containing booleans whether }
        """
        mut_status = defaultdict(list)
        has_cn_of_at_least_level = self.cnv.has_cn_of_at_least_level
        for gene_name, cats in categories.items():
            low_mut_status = [pd.Series(False, index=self.columns)]
            high_mut_status = [pd.Series(False, index=self.columns)]
            if "snv" in cats:
                high_mut_status.append(self.snv.has_snv.loc[gene_name])
            if "amp" in cats:
                high_mut_status.append(has_cn_of_at_least_level[self.cnv.high_amp_threshold].loc[gene_name])
                # high_mut_status.append(has_cn_of_at_least_level[self.cnv.mid_amp_threshold].loc[gene_name])
                for t in self.cnv.amp_thresholds:  # count high amps/dels to the low amp/del numbers (at least low)
                    low_mut_status.append(has_cn_of_at_least_level[t].loc[gene_name])
            if "del" in cats:
                high_mut_status.append(has_cn_of_at_least_level[self.cnv.high_del_threshold].loc[gene_name])
                # high_mut_status.append(has_cn_of_at_least_level[self.cnv.mid_del_threshold].loc[gene_name])
                for t in self.cnv.del_thresholds:  # count high amps/dels to the low amp/del numbers (at least low)
                    low_mut_status.append(has_cn_of_at_least_level[t].loc[gene_name])

            mut_status["low"].append(pd.concat(low_mut_status, axis=1).any(axis=1))
            mut_status["high"].append(pd.concat(high_mut_status, axis=1).any(axis=1))
        return {
            "low": pd.DataFrame(mut_status["low"], index=categories.keys(), columns=self.columns),
            "high": pd.DataFrame(mut_status["high"], index=categories.keys(), columns=self.columns),
        }

    def get_total_recurrence(self, categories=None):
        if categories is not None:
            return pd.concat([
                v.sum(axis=1).to_frame(k)
                for k, v in self.get_mutation_status(categories=categories).items()
            ], axis=1)
        else:
            has_high_mut = (self.snv.has_snv | self.cnv.has_high_cnv | self.cnv.has_mid_cnv).fillna(False)
            has_low_mut = (self.snv.has_snv | self.cnv.has_high_cnv | self.cnv.has_mid_cnv | self.cnv.has_low_cnv).fillna(False)
            has_mut = pd.concat([
                has_high_mut.astype(int).sum(axis=1).to_frame("high"),
                has_low_mut.astype(int).sum(axis=1).to_frame("low"),
            ], axis=1) / len(self.columns)
            return has_mut

    def get_total_recurrence_overall(self, categories=None, ground_truth_only=False):
        """ Measures percent of all patients that have a mutation of type 'low/high' """
        if categories is not None and self.ground_truth_genes is not None:
            gt_genes = pd.Index(set([g for glist in self.ground_truth_genes.values() for g in glist]))
            if not ground_truth_only or not len(gt_genes):
                gt_genes = self.genes
            has_mut = pd.Series({
                k: v.loc[gt_genes].any(axis=0).sum()
                for k, v in self.get_mutation_status(categories=categories).items()
            })
            return has_mut, len(self.columns)
        else:
            has_high_mut = (self.snv.has_snv | self.cnv.has_high_cnv | self.cnv.has_mid_cnv).fillna(False)
            has_low_mut = (self.snv.has_snv | self.cnv.has_high_cnv | self.cnv.has_mid_cnv | self.cnv.has_low_cnv).fillna(False)
            if ground_truth_only and self.ground_truth_genes is not None:
                gt_genes = pd.Index({g for glist in self.ground_truth_genes.values() for g in glist})
                gt_genes = gt_genes.intersection(has_high_mut.index)
                has_high_mut = has_high_mut.loc[gt_genes]
                has_low_mut = has_low_mut.loc[gt_genes]
            has_mut = pd.Series({
                "high": has_high_mut.any(axis=0).astype(int).sum(),
                "low": has_low_mut.any(axis=0).astype(int).sum(),
            })
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

    def _cnv_cell_distance(self, a, b):
        if pd.isna(a) or pd.isna(b):
            return 0.0

        if a == b:
            return 0.0

        if a == self.baseline or b == self.baseline:
            return 1.0

        if np.sign(a) == np.sign(b):
            return 0.35 * abs(abs(a) - abs(b))

        return 2.0 + 0.25 * abs(abs(a) - abs(b))

    def _pairwise_gene_distance_for_hclust(
        self,
        cnv_weight: float = 1.0,
        snv_weight: float = 1.0,
        genes=None,
        columns=None,
    ):
        genes = self.genes if genes is None else genes
        columns = self.columns if columns is None else columns
        cnv = (
            self.cnv.df
            .reindex(index=genes, columns=columns)
            .fillna(self.baseline)
        )

        snv = (
            self.snv.has_snv
            .reindex(index=genes, columns=columns)
            .fillna(False)
            .astype(int)
        )

        genes = cnv.index
        n = len(genes)
        D = np.zeros((n, n), dtype=float)

        for i in range(n):
            cnv_i = cnv.iloc[i].to_numpy()
            snv_i = snv.iloc[i].to_numpy()

            for j in range(i + 1, n):
                cnv_j = cnv.iloc[j].to_numpy()
                snv_j = snv.iloc[j].to_numpy()

                cnv_dist = np.mean([
                    self._cnv_cell_distance(a, b)
                    for a, b in zip(cnv_i, cnv_j)
                ])

                snv_dist = np.mean(snv_i != snv_j)

                D[i, j] = D[j, i] = (
                        cnv_weight * cnv_dist +
                        snv_weight * snv_dist
                )

        return pd.DataFrame(D, index=genes, columns=genes)

    def _pairwise_column_distance_for_hclust(
        self,
        cnv_weight: float = 1.0,
        snv_weight: float = 1.0,
        genes=None,
        columns=None,
    ):
        genes = self.genes if genes is None else genes
        columns = self.columns if columns is None else columns
        cnv = (
            self.cnv.df
            .reindex(index=genes, columns=columns)
            .fillna(self.baseline)
            .T
        )

        snv = (
            self.snv.has_snv
            .reindex(index=genes, columns=columns)
            .fillna(False)
            .astype(int)
            .T
        )

        columns = cnv.index
        n = len(columns)
        D = np.zeros((n, n), dtype=float)

        for i in range(n):
            cnv_i = cnv.iloc[i].to_numpy()
            snv_i = snv.iloc[i].to_numpy()

            for j in range(i + 1, n):
                cnv_j = cnv.iloc[j].to_numpy()
                snv_j = snv.iloc[j].to_numpy()

                cnv_dist = np.mean([
                    self._cnv_cell_distance(a, b)
                    for a, b in zip(cnv_i, cnv_j)
                ])

                snv_dist = np.mean(snv_i != snv_j)

                D[i, j] = D[j, i] = (
                        cnv_weight * cnv_dist +
                        snv_weight * snv_dist
                )

        return pd.DataFrame(D, index=columns, columns=columns)

    def _order_from_distance(self, D, method="average", optimal_ordering=True):
        if len(D) <= 2:
            return D.index

        condensed = squareform(D.to_numpy(), checks=False)

        if np.all(condensed == 0):
            return D.index

        Z = linkage(condensed, method=method)

        if optimal_ordering:
            Z = optimal_leaf_ordering(Z, condensed)

        return pd.Index(D.index[leaves_list(Z)], name=D.index.name)

    def sort_genes(self):
        order = self._compute_gene_order(genes=self.genes, columns=self.columns)
        if order is not None:
            self.genes = order

    def _compute_gene_order(self, genes, columns):
        """Return an ordered gene ``Index`` for ``genes`` using ``columns`` as the
        only sorting evidence. Pure: does not mutate ``self`` or the carriers.

        Passing ``genes=self.genes`` and ``columns=self.columns`` reproduces the
        legacy :meth:`sort_genes` result exactly (the ``reindex`` calls are
        identity operations in that case).
        """
        name = self.genes.name if self.genes is not None else MAF.gene_name
        genes = pd.Index(genes, name=name)

        if self.idx_order is not None:
            return pd.Index([c for c in self.idx_order if c in genes], name=name)

        if self.gene_sort_method == "hierarchical":
            D = self._pairwise_gene_distance_for_hclust(
                cnv_weight=self.cluster_cnv_weight,
                snv_weight=self.cluster_snv_weight,
                genes=genes,
                columns=columns,
            )
            return self._order_from_distance(D)

        has_snv = self.snv.has_snv.reindex(index=genes, columns=columns)
        has_high_cnv = self.cnv.has_high_cnv.reindex(index=genes, columns=columns)
        has_mid_cnv = self.cnv.has_mid_cnv.reindex(index=genes, columns=columns)
        has_low_cnv = self.cnv.has_low_cnv.reindex(index=genes, columns=columns)

        sorted_features = (
            (has_snv.astype(int) + has_high_cnv.astype(int) + has_mid_cnv.astype(int))
            .fillna(0)
            .sum(axis=1)
            .to_frame("mut_count")
            .join(has_snv.any(axis=1).to_frame("has_snv"))
            .join(has_low_cnv.any(axis=1).to_frame("has_low_cnv"))
            .join(has_low_cnv.astype(int).sum(axis=1).to_frame("low_cnv_count"))
            .sort_values(by=["mut_count", "has_snv", "has_low_cnv", "low_cnv_count"], ascending=True)
        )
        genes = sorted_features.index.to_list()
        # bring the interesting gene to the front if the heatmap:
        if self.interesting_gene is not None and self.interesting_gene in genes:
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
        if self.interesting_gene is not None and self.interesting_gene in cytobands.index:
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

        return pd.Index(genes, name=name)

    def sort_columns(self):
        order = self._compute_column_order(columns=self.columns, genes=self.genes)
        if order is not None:
            self.columns = order

    def _compute_column_order(self, columns, genes):
        """Return an ordered sample ``Index`` for ``columns`` using ``genes`` as the
        only sorting evidence. Pure: does not mutate ``self`` or the carriers.

        Passing ``columns=self.columns`` and ``genes=self.genes`` reproduces the
        legacy :meth:`sort_columns` result exactly.
        """
        name = self.columns.name
        columns = pd.Index(columns, name=name)

        if self.col_order is not None:
            return pd.Index([c for c in self.col_order if c in columns], name=name)

        if self.column_sort_method == "hierarchical":
            D = self._pairwise_column_distance_for_hclust(
                cnv_weight=self.cluster_cnv_weight,
                snv_weight=self.cluster_snv_weight,
                genes=genes,
                columns=columns,
            )
            return self._order_from_distance(D)

        def reidx(frame):
            return frame.reindex(index=genes, columns=columns)

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
            reidx(self.cnv.has_high_amp).astype(int),
            reidx(self.cnv.has_mid_amp).astype(int),
            reidx(self.cnv.has_high_del).astype(int),
            reidx(self.cnv.has_mid_del).astype(int),
            (~reidx(self.snv.has_snv) & (reidx(self.cnv.has_high_cnv) | reidx(self.cnv.has_mid_cnv))).fillna(False).astype(int),
            reidx(self.snv.has_snv).astype(int),
        ])
        if self.tmb is not None:
            tmb = self.tmb.reindex(index=columns)
            if SIF.tmb in tmb:
                has_burden = tmb[SIF.tmb].gt(0).astype(int).to_frame("has_burden").T
                burden = tmb[[SIF.tmb]].T
            else:
                has_burden = tmb.sum(axis=1).gt(0).astype(int).to_frame("has_burden").T
                burden = tmb.sum(axis=1).to_frame(SIF.tmb).T
        else:
            has_burden = pd.DataFrame()
            burden = pd.DataFrame()
        has_any_cnv = reidx(self.cnv.has_cnv).any(axis=0).astype(int).to_frame("has_cnv").T
        has_low_cnv = get_score([
            reidx(self.cnv.has_low_amp).astype(int),
            reidx(self.cnv.has_low_del).astype(int),
        ])
        comut_features = [df for df in [reidx(self.snv.deleteriousness_score), burden, has_low_cnv, has_any_cnv, has_burden, has_high_mut] if not df.empty]

        features = []
        for col in reversed(self.column_sort_by):
            if col == "COMUT":
                features += comut_features
            elif col == "TMB":
                features += [burden]
            elif col in self.meta_data_rows:
                if col in self.meta_data_rows_per_sample:
                    features.append(self.meta.df[col].reindex(index=columns).apply(lambda l: l[0] if len(l) else 0).to_frame(col).T)
                else:
                    features.append(self.meta.df[[col]].reindex(index=columns).T)
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
        return pd.Index(columns, name=name)

    def get_model_annotation(self):
        return pd.DataFrame(
            [[g in self.cnv_interesting_genes, g in self.snv_interesting_genes] for g in self.genes],
            index=self.genes,
            columns=["cnv", "snv"]
        )
