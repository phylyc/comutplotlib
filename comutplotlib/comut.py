from __future__ import annotations

from collections.abc import Sequence
from copy import deepcopy
import warnings
import numpy as np
import pandas as pd

from comutplotlib.comut_data import ComutData
from comutplotlib.comut_layout import ComutLayout
from comutplotlib.comut_panels import ComutPanels, meta_legend, control, block, grid_row_panel, grid_col_panel
from comutplotlib.comut_plotter import ComutPlotter
from comutplotlib.functional_effect import sort_functional_effects
from comutplotlib.mathutils import fold_change
from comutplotlib.mutation_annotation import MutationAnnotation as MutA
from comutplotlib.palette import Palette
from comutplotlib.sample_annotation import SampleAnnotation as SA

from comutplotlib.gistic import join_gistics
from comutplotlib.mark import join_marks
from comutplotlib.maf import join_mafs
from comutplotlib.seg import join_segs
from comutplotlib.sif import join_sifs


class Comut(object):
    """ Caption:
    This comutation plot (central panel) visualizes the mutation landscape across a patient cohort. Rows represent genes, and columns correspond to patients. Each cell indicates a gene’s mutation status in a patient: rectangles denote copy-number variations (CNVs), and ellipses indicate short nucleotide variations (SNVs), with multiple SNVs shown as subdivided wedges. Colors encode mutation types and functional effects.

    The top panel displays tumor mutation burden (TMB) per patient, with high TMB (≥10/Mb) highlighted in red. The mutational signature panel shows the relative fraction of exposures to different mutational signatures for each patient or sample. The left panel summarizes mutation recurrence, showing SNV and CNV frequencies per gene, annotated with recurrent protein alterations. It also reports the percentage of patients with high-level CNVs, supplemented by low-level CNVs in brackets, along with the total percentage of patients carrying either an SNV or a high-level CNV in the gene.

    The bottom panel presents patient- and sample-level metadata. For patients with multiple samples, metadata cells are subdivided accordingly.
    """

    def __init__(
        self,
        output: str = "./comut.pdf",

        maf: list[str] | None = None,
        maf_pool_as: dict | None = None,

        seg: list[str] | None = None,
        gistic: list[str] | None = None,

        mark: list[str] | None = None,
        mark_label: str = "Epigenetic\nMarks",

        signatures: list[str] | None = None,
        group_signatures_by_etiology: bool = False,

        cohort_label: str | None = None,
        cohort_tag: str | None = None,

        model_significances: list[str] | None = None,
        model_names: Sequence[str] = (),

        sif: list[str] | None = None,
        meta_data_rows: Sequence[str] = (),
        meta_data_rows_per_sample: Sequence[str] = (),

        gene_meta_data: str | None = None,
        gene_meta_columns: list[str] | None = None,

        control_maf: list[str] | None = None,
        control_seg: list[str] | None = None,
        control_gistic: list[str] | None = None,
        control_mark: list[str] | None = None,
        control_sif: list[str] | None = None,
        control_signatures: list[str] | None = None,
        control_cohort_label: str | None = None,
        control_cohort_tag: str | None = None,
        control_position: str = "right",

        drop_empty_columns: bool = False,

        by: str = MutA.patient,

        column_order: tuple[str] | None = None,
        control_column_order: tuple[str] | None = None,
        index_order: tuple[str] | None = None,
        column_sort_by: tuple[str] = ("COMUT",),
        index_sort_by: str = "COMUT",
        index_sort_by_direction: str = "ascending",
        sort_method: str | None = None,

        interesting_gene: str | None = None,
        interesting_gene_comut_percent_threshold: float | None = None,
        interesting_genes: set | None = None,
        snv_interesting_genes: set | None = None,
        cnv_interesting_genes: set | None = None,
        total_recurrence_threshold: float | None = None,
        scale_recurrence: bool = False,
        recurrence_categories: dict[str, list[str] | dict[str, list[str]]] | None = None,
        snv_recurrence_threshold: int = 5,
        min_fold_change: float | None = None,
        max_fold_change: float | None = None,
        collapse_cytobands: bool = False,
        collapse_cytobands_range: int = 0,

        low_amp_threshold: int | float = 1,
        mid_amp_threshold: int | float = 1.5,
        high_amp_threshold: int | float = 2,
        baseline: int | float = 0,
        low_del_threshold: int | float = -1,
        mid_del_threshold: int | float = -1.5,
        high_del_threshold: int | float = -2,
        show_low_level_cnvs: bool = True,

        panels_to_plot: Sequence[str] = (),
        palette: dict[str, dict] | None = None,
        ground_truth_genes: dict[str, list[str]] | None = None,  # todo: refactor as a palette class
        max_xfigsize: int | None = None,
        max_xfigsize_scale: float = 1,
        label_columns: bool = False,

        column_group_by: Sequence[str] | None = None,
        hide_grouped_meta_data: bool = False,
        gene_groups: dict[str, list[str]] | None = None,
        group_order: Sequence[str] | None = None,
        column_group_labels: Sequence[str] | None = None,
        drop_ungrouped_genes: bool = False,
        na_group_label: str = "NA",
        other_gene_group_label: str = "Other",
        meta_data_position: str = "bottom",
        **kwargs
    ):
        # Configuration is stored on typed attributes so the pipeline can run as
        # discrete, independently callable stages (load -> preprocess ->
        # build_layout -> render) without a monolithic constructor or hidden disk
        # writes. Keeping one attribute per argument preserves the declared types
        # (a single config dict would collapse them into one big union).
        self._output = output
        self._palette = palette
        self._cohort_label = cohort_label
        # The tag is a short cohort identifier used for in-plot annotations (e.g.
        # the fold-change arrows), kept separate from the label used for titles.
        self._cohort_tag = cohort_tag if cohort_tag is not None else "case"
        self._maf = maf
        self._maf_pool_as = maf_pool_as
        self._seg = seg
        self._gistic = gistic
        self._mark = mark
        self._mark_label = mark_label
        self._signatures = signatures
        self._group_signatures_by_etiology = group_signatures_by_etiology
        self._sif = sif
        self._meta_data_rows = meta_data_rows
        self._meta_data_rows_per_sample = meta_data_rows_per_sample
        if meta_data_position not in {"top", "bottom"}:
            raise ValueError("meta_data_position must be either 'top' or 'bottom'.")
        self._meta_data_position = meta_data_position
        self._gene_meta_data = gene_meta_data
        self._gene_meta_columns = gene_meta_columns
        self._model_significances = model_significances
        self._model_names = model_names
        self._control_maf = control_maf
        self._control_seg = control_seg
        self._control_gistic = control_gistic
        self._control_mark = control_mark
        self._control_sif = control_sif
        self._control_signatures = control_signatures
        self._control_cohort_label = control_cohort_label
        self._control_cohort_tag = control_cohort_tag if control_cohort_tag is not None else "control"
        if control_position not in {"left", "right"}:
            raise ValueError("control_position must be either 'left' or 'right'.")
        self._control_position = control_position
        # When ``control_position == "left"`` the case and control cohorts are
        # mirrored around the central fold-change panel. This single flag drives
        # column reversal, recurrence axis inversion and layout placement.
        self._control_on_left = control_position == "left"
        self._drop_empty_columns = drop_empty_columns
        self._by = by
        self._column_order = column_order
        self._control_column_order = control_column_order
        self._index_order = index_order
        self._column_sort_by = column_sort_by
        self._index_sort_by = index_sort_by
        self._index_sort_by_direction = index_sort_by_direction
        self._sort_method = sort_method
        self._interesting_gene = interesting_gene
        self._interesting_gene_comut_percent_threshold = interesting_gene_comut_percent_threshold
        self._interesting_genes = interesting_genes
        self._snv_interesting_genes = snv_interesting_genes
        self._cnv_interesting_genes = cnv_interesting_genes
        self._total_recurrence_threshold = total_recurrence_threshold
        self._scale_recurrence = scale_recurrence
        self._recurrence_categories = (
            recurrence_categories if recurrence_categories is not None
            else {"global": ["snv", "amp", "del"]}
        )
        self._snv_recurrence_threshold = snv_recurrence_threshold
        self._min_fold_change = min_fold_change
        self._max_fold_change = max_fold_change
        self._collapse_cytobands = collapse_cytobands
        self._collapse_cytobands_range = collapse_cytobands_range
        self._low_amp_threshold = low_amp_threshold
        self._mid_amp_threshold = mid_amp_threshold
        self._high_amp_threshold = high_amp_threshold
        self._baseline = baseline
        self._low_del_threshold = low_del_threshold
        self._mid_del_threshold = mid_del_threshold
        self._high_del_threshold = high_del_threshold
        self._show_low_level_cnvs = show_low_level_cnvs
        self._panels_to_plot = list(panels_to_plot)  # mutated during panel pruning
        self._ground_truth_genes = ground_truth_genes
        self._max_xfigsize = max_xfigsize
        self._max_xfigsize_scale = max_xfigsize_scale
        self._label_columns = label_columns

        # Grid sub-panel configuration (see GRID_SUBPANELS_PLAN.md).
        self._column_group_by = list(column_group_by) if column_group_by is not None else []
        self._hide_grouped_meta_data = hide_grouped_meta_data
        self._gene_groups = gene_groups
        self._group_order = list(group_order) if group_order is not None else None
        # Optional display labels for the grid columns, applied positionally in
        # grid-column order (see ``grid.apply_column_group_labels``).
        self._column_group_labels = list(column_group_labels) if column_group_labels is not None else None
        self._drop_ungrouped_genes = drop_ungrouped_genes
        self._na_group_label = na_group_label
        self._other_gene_group_label = other_gene_group_label

    @property
    def left_cohort(self):
        """Cohort data drawn on the left of the central fold-change panel."""
        return self.control if self._control_on_left else self.case

    @property
    def right_cohort(self):
        """Cohort data drawn on the right of the central fold-change panel."""
        return self.case if self._control_on_left else self.control

    def load(self):
        """Read inputs and build the case / control / joint data models."""
        self.plotter = ComutPlotter(
            output=self._output,
            extra_palette=self._palette["global"] if self._palette is not None else None
        )

        self.case = ComutData(
            cohort_name=self._cohort_label,
            maf_paths=self._maf,
            maf_pool_as=self._maf_pool_as,
            seg_paths=self._seg,
            gistic_paths=self._gistic,
            mark_paths=self._mark,
            signatures_paths=self._signatures,
            group_signatures_by_etiology=self._group_signatures_by_etiology,
            sif_paths=self._sif,
            meta_data_rows=self._meta_data_rows,
            meta_data_rows_per_sample=self._meta_data_rows_per_sample,
            drop_empty_columns=self._drop_empty_columns,
            by=self._by,
            column_order=self._column_order,
            index_order=self._index_order,
            column_sort_by=self._column_sort_by,
            sort_method=self._sort_method,
            interesting_gene=self._interesting_gene,
            interesting_gene_comut_percent_threshold=self._interesting_gene_comut_percent_threshold,
            interesting_genes=self._interesting_genes,
            ground_truth_genes=self._ground_truth_genes,
            snv_interesting_genes=self._snv_interesting_genes,
            cnv_interesting_genes=self._cnv_interesting_genes,
            total_recurrence_threshold=self._total_recurrence_threshold,
            snv_recurrence_threshold=self._snv_recurrence_threshold,
            low_amp_threshold=self._low_amp_threshold,
            mid_amp_threshold=self._mid_amp_threshold,
            high_amp_threshold=self._high_amp_threshold,
            baseline=self._baseline,
            low_del_threshold=self._low_del_threshold,
            mid_del_threshold=self._mid_del_threshold,
            high_del_threshold=self._high_del_threshold,
            show_low_level_cnvs=self._show_low_level_cnvs,
            column_group_by=self._column_group_by,
            gene_groups=self._gene_groups,
            group_order=self._group_order,
            column_group_labels=self._column_group_labels,
            drop_ungrouped_genes=self._drop_ungrouped_genes,
            na_group_label=self._na_group_label,
            other_gene_group_label=self._other_gene_group_label,
        )
        self.case.preprocess()

        self.control = ComutData(
            cohort_name=self._control_cohort_label,
            maf_paths=self._control_maf,
            maf_pool_as=self._maf_pool_as,
            seg_paths=self._control_seg,
            gistic_paths=self._control_gistic,
            mark_paths=self._control_mark,
            signatures_paths=self._control_signatures,
            group_signatures_by_etiology=self._group_signatures_by_etiology,
            sif_paths=self._control_sif,
            meta_data_rows=self._meta_data_rows,
            meta_data_rows_per_sample=self._meta_data_rows_per_sample,
            drop_empty_columns=self._drop_empty_columns,
            by=self._by,
            column_order=self._control_column_order,
            index_order=self.case.genes,
            column_sort_by=self._column_sort_by,
            sort_method=self._sort_method,
            interesting_gene=self.case.interesting_gene,
            interesting_gene_comut_percent_threshold=self._interesting_gene_comut_percent_threshold,
            interesting_genes=self.case.interesting_genes,
            ground_truth_genes=self.case.ground_truth_genes,
            snv_interesting_genes=self.case.snv_interesting_genes,
            cnv_interesting_genes=self.case.cnv_interesting_genes,
            snv_recurrence_threshold=self._snv_recurrence_threshold,
            low_amp_threshold=self._low_amp_threshold,
            mid_amp_threshold=self._mid_amp_threshold,
            high_amp_threshold=self._high_amp_threshold,
            baseline=self._baseline,
            low_del_threshold=self._low_del_threshold,
            mid_del_threshold=self._mid_del_threshold,
            high_del_threshold=self._high_del_threshold,
            show_low_level_cnvs=self._show_low_level_cnvs,
            column_group_by=self._column_group_by,
            group_order=self._group_order,
            column_group_labels=self._column_group_labels,
            na_group_label=self._na_group_label,
        )
        self.control.preprocess()

        self.model_significance = pd.DataFrame.from_dict(
            {
                name: pd.read_csv(path_to_file, index_col=0, sep="\t")
                for name, path_to_file in zip(self._model_names, self._model_significances)
            }
        ) if self._model_significances is not None else None
        self.model_names = self._model_names

        self.joint = deepcopy(self.case)
        self.joint.gistic = join_gistics([self.case.gistic, self.control.gistic])
        self.joint.mark = join_marks([self.case.mark, self.control.mark])
        self.joint.seg = join_segs([self.case.seg, self.control.seg])
        self.joint.maf = join_mafs([self.case.maf, self.control.maf])
        self.joint.sif = join_sifs([self.case.sif, self.control.sif])
        signatures = [s for s in [self.case.signatures, self.control.signatures] if s is not None]
        self.joint.signatures = pd.concat(signatures) if len(signatures) else None
        self.joint.preprocess()

        self.case.meta.reindex(index=self.joint.meta.rows)
        self.control.meta.reindex(index=self.joint.meta.rows)

    def preprocess(self):
        """Filter/sort genes, build colour maps, and prune unavailable panels."""
        self.amp_thresholds = [self.joint.high_amp_threshold, self.joint.mid_amp_threshold, self.joint.low_amp_threshold]
        self.del_thresholds = [self.joint.high_del_threshold, self.joint.mid_del_threshold, self.joint.low_del_threshold]

        self.scale_recurrence = self._scale_recurrence
        self.recurrence_categories = self.get_recurrence_categories_by_gene(self._recurrence_categories)
        categories = {
            "snv": self.joint.snv.effects,
            "amp": self.amp_thresholds,
            "del": self.del_thresholds
        }
        recurrence_fold_change = None
        good_genes = list(self.joint.genes)
        if len(self.control.columns) > 0:
            recurrence_fold_change = self.get_recurrence_fold_change_by_gene(alpha_ci=0.1, base=2)
            good_genes = []
            for g, row in recurrence_fold_change["mean"].iterrows():
                relevant_columns = [c for cat in self.recurrence_categories.get(g, categories.keys()) for c in categories[cat]]
                is_good = True
                if self._min_fold_change is not None:
                    is_good &= (row[relevant_columns] > np.log(self._min_fold_change) / np.log(2)).any()
                if self._max_fold_change is not None:
                    is_good &= (row[relevant_columns] < np.log(self._max_fold_change) / np.log(2)).any()
                if is_good:
                    good_genes.append(g)

        for data in [self.case, self.control, self.joint]:
            data.genes = pd.Index(good_genes, name=MutA.gene_name)
            data.sort_columns()
            data.reindex_data()

        self.recurrence_categories = self.get_recurrence_categories_by_gene(self._recurrence_categories)

        # Optionally re-order the gene (row) index by case-vs-control fold change
        # (--index-sort-by). An explicit --index-order always wins; "COMUT" keeps
        # the standard mutation-based ordering established above.
        if self._index_order is None and self._index_sort_by != "COMUT":
            ordered_genes = self._order_genes_by_fold_change(
                genes=list(self.joint.genes),
                recurrence_fold_change=recurrence_fold_change,
                categories=categories,
            )
            for data in [self.case, self.control, self.joint]:
                # Set idx_order so the grid path honours the same ordering when it
                # re-sorts genes within each gene group.
                data.idx_order = pd.Index(ordered_genes, name=MutA.gene_name)
                data.genes = pd.Index(ordered_genes, name=MutA.gene_name)
                data.reindex_data()

        self.tmb_cmap = self.plotter.palette.get_tmb_cmap(self.joint.tmb)
        self.snv_cmap = self.plotter.palette.get_snv_cmap(self.joint.snv)
        self.cnv_cmap, self.cnv_names = self.plotter.palette.get_cnv_cmap(self.joint.cnv)
        self.epi_alphas = self.joint.epi.legend_alphas
        self.signatures_cmap = self.plotter.palette.get_signatures_cmap(self.joint.signatures)
        self.meta_cmaps = self.plotter.palette.get_meta_cmaps(self.joint.meta)
        if self._palette is not None:
            for col, pal in self._palette["local"].items():
                if col in self.meta_cmaps:
                    self.meta_cmaps[col] |= pal
                else:
                    self.meta_cmaps[col] = pal
        self.meta_cmaps_condensed = self.plotter.palette.condense(self.meta_cmaps)

        grid_enabled = self.case.grouping_enabled

        if len(self.control.columns) > 0 and not grid_enabled:
            # The cohort drawn on the LEFT has its columns reversed so its densest
            # columns meet the central recurrence panels (mirrors the right cohort,
            # whose natural order already points its densest columns inward).
            self.left_cohort.columns = self.left_cohort.columns[::-1]
            self.left_cohort.reindex_data()

        def remove(col):
            if col in self._panels_to_plot:
                self._panels_to_plot.remove(col)

        if "tmb" in self.tmb_cmap.keys():
            remove(ComutPanels.tmb_legend)
        if self.joint.tmb is None:
            remove(ComutPanels.tmb)
            remove(ComutPanels.tmb_legend)
        if self.joint.signatures is None:
            remove(ComutPanels.mutational_signatures)
            remove(ComutPanels.mutational_signatures_legend)
        if self.joint.epi.empty:
            remove(ComutPanels.epi_legend)

        if self._gene_meta_data is not None:
            self.gene_meta_data = (
                pd.read_csv(self._gene_meta_data, sep="\t")
                    .set_index("Gene")
                    .reindex(index=self.case.genes, columns=self._gene_meta_columns)
                    .replace(0, np.nan)
                    .dropna(axis=1, how="all")
                    .fillna(0)
            )
            self.gene_meta_data = self.gene_meta_data[self.gene_meta_data.sum().sort_values().index]
        else:
            self.gene_meta_data = None

        # Build the grid partition last, once genes/columns/meta are finalized.
        # The case defines the shared gene (row) axis; the control adopts it and
        # derives its own column groups.
        if grid_enabled:
            case_gene_groups = self.case.apply_grid_sort()
            if not self.control.empty:
                # Hand over the case's resolved labels as a key -> label map so
                # both cohorts label identical column groups identically, even if
                # the control resolves them in a different order (group sizes
                # feed into the ordering).
                if self._column_group_labels:
                    self.control.column_group_labels = self.case.column_group_label_map()
                self.control.apply_grid_sort(gene_groups=case_gene_groups)
            else:
                self.control.grid = None
            # Reproduce the inward-out case/control ordering for grid layouts: with
            # a control cohort present, invert the LEFT cohort's grid columns (both
            # the grid column order and the panel columns within each grid column)
            # so it reads inside-out towards the central fold-change panel (mirrors
            # the non-grid ``left_cohort.columns[::-1]`` above).
            if len(self.control.columns) > 0:
                self.left_cohort.invert_grid_columns()

        if self._hide_grouped_meta_data and self._column_group_by:
            grouped_rows = set(self._column_group_by)
            for data in [self.case, self.control, self.joint]:
                visible_rows = [row for row in data.meta.rows if row not in grouped_rows]
                data.meta.reindex(index=visible_rows)

            # Metadata palettes are built before grid partitioning because they do
            # not otherwise depend on the grid. Remove hidden rows from both the
            # table and its legends after grouping has consumed their values.
            visible_rows = set(self.joint.meta.rows)
            self.meta_cmaps = {
                row: cmap for row, cmap in self.meta_cmaps.items()
                if row in visible_rows
            }
            self.meta_cmaps_condensed = self.plotter.palette.condense(self.meta_cmaps)

    def build_layout(self):
        """Compute panel dimensions and instantiate the figure layout."""
        n_meta_genes = self.gene_meta_data.shape[1] if self.gene_meta_data is not None else 0

        n_genes, n_samples_case, n_meta_case = self.case.get_dimensions()
        _, n_samples_control, n_meta_control = self.control.get_dimensions()

        n_meta = max(n_meta_case, n_meta_control)

        grid = self.case.grid if (self.case.grid is not None and not self.case.grid.is_trivial()) else None
        control_column_group_sizes = (
            [len(g) for g in self.control.grid.column_groups]
            if (grid is not None and self.control.grid is not None and not self.control.grid.is_trivial())
            else None
        )
        control_column_group_labels = (
            [g.label for g in self.control.grid.column_groups]
            if (grid is not None and self.control.grid is not None and not self.control.grid.is_trivial())
            else None
        )

        self.layout = ComutLayout(
            panels_to_plot=self._panels_to_plot,
            max_xfigsize=self._max_xfigsize,
            max_xfigsize_scale=self._max_xfigsize_scale,
            n_genes=n_genes,
            n_samples=n_samples_case,
            n_samples_control=n_samples_control,
            n_meta=n_meta,
            n_meta_genes=n_meta_genes,
            label_columns=self._label_columns,
            meta_data_position=self._meta_data_position,
            tmb_cmap=self.tmb_cmap,
            snv_cmap=self.snv_cmap,
            cnv_cmap=self.cnv_cmap,
            epi_alphas=self.epi_alphas,
            mutsig_cmap=self.signatures_cmap,
            meta_cmaps=self.meta_cmaps_condensed,
            grid=grid,
            control_column_group_sizes=control_column_group_sizes,
            control_column_group_labels=control_column_group_labels,
            control_position=self._control_position,
        )

    def get_recurrence_categories_by_gene(self, categories, ref_cohort=None) -> dict[str, list[str]]:
        ref_cohort = (
            self.case if ref_cohort in [None, "case", self.case]
            else self.control if ref_cohort in ["control", self.control]
            else self.joint
        )
        cnv_recurrence = ref_cohort.cnv.get_num_patients_by_gene_by_cn_level().reindex(index=ref_cohort.genes).fillna(0)
        snv_recurrence = ref_cohort.snv.get_num_patients_by_gene_by_effect().reindex(index=ref_cohort.genes).fillna(0)

        recurrence_categories_by_gene = {}
        for i, (gene_name, row) in enumerate(cnv_recurrence.iterrows()):
            non_nan_cnv_counter = {k: v for k, v in row.dropna().items() if v > 0}
            non_nan_snv_counter = {k: v for k, v in snv_recurrence.loc[gene_name].items() if v > 0}

            non_nan_counts = non_nan_cnv_counter | non_nan_snv_counter

            if not len(non_nan_cnv_counter) and not len(non_nan_snv_counter):
                continue

            present_amp_thresholds = sorted([t for t in self.amp_thresholds if t in non_nan_cnv_counter], reverse=True)
            present_del_thresholds = sorted([t for t in self.del_thresholds if t in non_nan_cnv_counter])
            present_snv_effects = sort_functional_effects([e for e in non_nan_snv_counter.keys()], ascending=False)

            amp_max = max(non_nan_counts[t] for t in present_amp_thresholds) if len(present_amp_thresholds) else 0
            del_max = max(non_nan_counts[t] for t in present_del_thresholds) if len(present_del_thresholds) else 0
            snv_max = sum(non_nan_counts[t] for t in present_snv_effects) if len(present_snv_effects) else 0

            idx = np.argmax([amp_max, del_max, snv_max])
            max_cat = ["amp", "del", "snv"][idx]

            categories_for_gene = categories.get("local", {}).get(gene_name, categories.get("global"))
            recurrence_categories_by_gene[gene_name] = [max_cat if c == "max" else c for c in categories_for_gene]

        return recurrence_categories_by_gene

    def get_recurrence_fold_change_by_gene(self, alpha_ci=0.2, base=2):
        case = pd.concat([
            self.case.snv.get_num_patients_by_gene_by_effect().reindex(columns=self.joint.snv.effects).fillna(0),
            self.case.cnv.get_num_patients_by_gene_of_at_least_cn_level()
        ], axis=1)
        n_case = len(self.case.columns)
        control = pd.concat([
            self.control.snv.get_num_patients_by_gene_by_effect().reindex(columns=self.joint.snv.effects).fillna(0),
            self.control.cnv.get_num_patients_by_gene_of_at_least_cn_level()
        ], axis=1)
        n_control = len(self.control.columns)
        return fold_change(case, n_case, control, n_control, alpha_ci=alpha_ci, base=base)

    def get_total_recurrence_fold_change(self, alpha_ci=0.2, base=2):
        case, n_case = self.case.get_total_recurrence_overall(categories=self.recurrence_categories)
        control, n_control = self.control.get_total_recurrence_overall(categories=self.recurrence_categories)
        return fold_change(case, n_case, control, n_control, alpha_ci=alpha_ci, base=base)

    # --- Index (gene) ordering by fold change (``--index-sort-by``) ----------

    @staticmethod
    def _peak_fold_change(values):
        """Return the fold-change value with the largest magnitude (sign kept).

        NaN entries (event absent in either cohort) are ignored; an empty/all-NaN
        input yields NaN so the gene sinks to the bottom of the ordering.
        """
        values = pd.Series(values).dropna()
        if values.empty:
            return np.nan
        return values.loc[values.abs().idxmax()]

    def _total_fold_change_scores(self, genes):
        """Per-gene total fold change across *all* categories.

        Uses :meth:`ComutData.get_total_recurrence` (categories=None), which yields
        the per-gene fraction of patients carrying any 'high' or any 'low' mutation.
        These fractions are converted back to counts and fed through
        :func:`fold_change`; the larger-magnitude of the 'high'/'low' fold change is
        taken as the gene's score.
        """
        n_case = len(self.case.columns)
        n_control = len(self.control.columns)
        case_counts = (self.case.get_total_recurrence(categories=None) * n_case).round()
        control_counts = (self.control.get_total_recurrence(categories=None) * n_control).round()
        idx = case_counts.index.union(control_counts.index)
        case_counts = case_counts.reindex(idx).fillna(0)
        control_counts = control_counts.reindex(idx).fillna(0)
        mean_fc = fold_change(case_counts, n_case, control_counts, n_control, alpha_ci=0.1, base=2)["mean"]
        return {
            g: mean_fc.loc[g, "low"] if g in mean_fc.index else np.nan
            for g in genes
        }

    def _selected_fold_change_scores(self, genes, recurrence_fold_change, categories):
        """Per-gene fold change restricted to the gene's selected recurrence
        categories, further filtered to 'high' (SNV + high/mid CNV) or 'low'
        (low CNV) columns depending on ``--index-sort-by``.
        """
        mean_fc = recurrence_fold_change["mean"]
        snv_cols = set(categories["snv"])
        high_cnv = {
            self.joint.high_amp_threshold, self.joint.mid_amp_threshold,
            self.joint.high_del_threshold, self.joint.mid_del_threshold,
        }
        low_cnv = {self.joint.low_amp_threshold, self.joint.low_del_threshold}

        scores = {}
        for g in genes:
            cats = self.recurrence_categories.get(g, list(categories))
            cols = [c for cat in cats for c in categories[cat]]
            if self._index_sort_by == "high-fold-change":
                cols = [c for c in cols if c in snv_cols or c in high_cnv]
            elif self._index_sort_by == "low-fold-change":
                cols = [c for c in cols if c in low_cnv]
            cols = [c for c in cols if c in mean_fc.columns]
            scores[g] = (
                self._peak_fold_change(mean_fc.loc[g, cols])
                if g in mean_fc.index and cols else np.nan
            )
        return scores

    def _order_genes_by_fold_change(self, genes, recurrence_fold_change, categories):
        """Return ``genes`` re-ordered by the ``--index-sort-by`` fold-change score.

        The incoming ``genes`` order (COMUT) is preserved as a stable tie-break.
        Falls back to the original order (with a warning) when no control cohort is
        available, since fold change is undefined without one.
        """
        genes = list(genes)
        if len(self.control.columns) == 0:
            warnings.warn(
                f"--index-sort-by={self._index_sort_by!r} requires a control cohort; "
                "falling back to the standard COMUT gene ordering.",
                stacklevel=2,
            )
            return genes

        if self._index_sort_by == "fold-change":
            scores = self._total_fold_change_scores(genes)
        else:
            scores = self._selected_fold_change_scores(genes, recurrence_fold_change, categories)

        ascending = self._index_sort_by_direction == "ascending"
        na_position = "first" if ascending else "last"
        ordered = (
            pd.Series([scores[g] for g in genes], index=genes)
            .sort_values(ascending=ascending, na_position=na_position, kind="stable")
        )
        return ordered.index.tolist()

    def make_comut(self):
        self.load()
        self.preprocess()
        self.build_layout()
        return self.render()

    def render(self):
        # Persist the finalized column/gene ordering next to the figure output.
        self.case.save(out_dir=self.plotter.out_dir, name=self.plotter.file_name + ".case")
        self.control.save(out_dir=self.plotter.out_dir, name=self.plotter.file_name + ".control")

        if self.case.grid is not None and not self.case.grid.is_trivial():
            return self._render_grid()

        self.layout.add_panels()

        def plot_comut_gen(data):
            def plot_comut(ax):
                self.plotter.plot_cnv_heatmap(
                    ax=ax,
                    cnv=data.cnv.df,
                    cnv_cmap=self.cnv_cmap,
                    inter_heatmap_linewidth=self.layout.inter_heatmap_linewidth,
                    aspect_ratio=self.layout.aspect_ratio
                )
                self.plotter.plot_epi_heatmap(
                    ax=ax,
                    epi=data.epi.df,
                    inter_heatmap_linewidth=self.layout.inter_heatmap_linewidth,
                    aspect_ratio=self.layout.aspect_ratio
                )
                self.plotter.plot_snv_heatmap(
                    ax=ax,
                    snv=data.snv.df,
                    snv_cmap=self.snv_cmap
                )
                self.plotter.plot_heatmap_layout(
                    ax=ax,
                    cna=data.cnv.df,
                    labelbottom=(
                        self.layout.show_patient_names
                        and (data.meta is None or self._meta_data_position == "top")
                    )
                )
            return plot_comut

        def meta_data_color(column, value):
            cmap = self.meta_cmaps[column]
            if isinstance(cmap, Palette):
                return cmap[value] if value in cmap else self.plotter.palette.get(value, self.plotter.palette.white)
            else:
                _cmap, _norm = cmap
                return _cmap(_norm(value))

        if self.joint.tmb is None:
            tmb_ymin, tmb_ymax = 5 * 1e-1, 1.01 * 1e2
        elif SA.tmb in self.joint.tmb.columns:
            # tmb_ymin = min(10 ** np.floor(np.log10(self.joint.tmb[SA.tmb].quantile(0.15))), 5 * 1e-1)
            tmb_ymin = min(self.joint.tmb[SA.tmb].quantile(0.1), 5 * 1e-1)
            tmb_ymax = np.clip(self.joint.tmb[SA.tmb].max(), a_min=1.01 * 1e2, a_max=1e4)
        else:
            tmb_ymin = 0
            tmb_ymax = np.clip(self.joint.tmb.sum(axis=1).max(), a_min=1.01 * 1e2, a_max=1e4)

        # case_cnv_recurrence = self.case.cnv.get_prevalence()
        # case_snv_recurrence = self.case.snv.get_patient_recurrence()
        # control_cnv_recurrence = self.control.cnv.get_prevalence()
        # control_snv_recurrence = self.control.snv.get_patient_recurrence()
        #
        # recurrence_max_xlim = pd.concat([
        #     case_cnv_recurrence.apply(lambda d: sum([v for k, v in d.items() if k in amp_thresholds])),
        #     case_cnv_recurrence.apply(lambda d: sum([v for k, v in d.items() if k in del_thresholds])),
        #     case_snv_recurrence.apply(lambda d: sum([v for k, v in d.items() if k is not None])),
        #     control_cnv_recurrence.apply(lambda d: sum([v for k, v in d.items() if k in amp_thresholds])),
        #     control_cnv_recurrence.apply(lambda d: sum([v for k, v in d.items() if k in del_thresholds])),
        #     control_snv_recurrence.apply(lambda d: sum([v for k, v in d.items() if k is not None]))
        #  ], axis=1).max(axis=1).max() if self.scale_recurrence else max(len(self.case.columns), len(self.control.columns))

        has_control = len(self.control.columns) > 0
        # Effective mirror flag (only meaningful with a control cohort).
        control_on_left = self._control_on_left and has_control
        if has_control:
            self.layout.set_plot_func(
                ComutPanels.recurrence_fold_change,
                self.plotter.plot_recurrence_fold_change,
                fold_change=self.get_recurrence_fold_change_by_gene(alpha_ci=0.2, base=2),
                genes=self.case.genes,
                cnv_cmap=self.cnv_cmap,
                categories=self.recurrence_categories,
                effects=self.joint.snv.effects,
                amp_thresholds=self.amp_thresholds,
                del_thresholds=self.del_thresholds,
                label_x=ComutPanels.total_recurrence_fold_change not in self.layout.panels,
                pad=0.01,
                case_on_left=not control_on_left,
                case_label=self._cohort_tag,
                control_label=self._control_cohort_tag,
            )
            self.layout.set_plot_func(
                ComutPanels.total_recurrence_fold_change,
                self.plotter.plot_total_recurrence_fold_change,
                fold_change=self.get_total_recurrence_fold_change(alpha_ci=0.2, base=2),
                shared_x_ax=self.layout.panels.get(ComutPanels.recurrence_fold_change).ax if ComutPanels.recurrence_fold_change in self.layout.panels and self.layout.panels.get(ComutPanels.recurrence_fold_change).plot_func is not None else None,
                pad=0.01,
            )

        tmb_ref = {"ax": None}
        cohorts = [(self.case, "", True, has_control), (self.control, " control", False, False)]
        if control_on_left:
            cohorts = cohorts[::-1]
        for data, label, is_case, special in cohorts:
            if len(data.columns) == 0:
                continue

            # Geometric side of this cohort: the left cohort keeps the legacy
            # (non-inverted) recurrence axis, the right cohort mirrors it.
            on_left = is_case != control_on_left

            self.layout.set_plot_func(
                ComutPanels.comutation + label,
                plot_comut_gen(data)
            )

            self.layout.set_plot_func(
                ComutPanels.cohort_label + label,
                self.plotter.plot_cohort_label,
                label=data.name
            )

            # self.layout.set_plot_func("coverage", self.plotter.plot_coverage)
            self.layout.set_plot_func(
                ComutPanels.mutational_signatures + label,
                self.plotter.plot_signatures,
                signatures=data.signatures,
                signatures_cmap=self.signatures_cmap,
                add_ylabel=on_left,
            )
            tmb_name = ComutPanels.tmb + label
            self.layout.set_plot_func(
                tmb_name,
                self.plotter.plot_tmb,
                tmb=data.tmb,
                ytickpad=0,
                shared_y_ax=tmb_ref["ax"],
                aspect_ratio=self.layout.aspect_ratio,
                ymin=tmb_ymin,
                ymax=tmb_ymax,
                median_on_right=on_left,
            )
            if tmb_ref["ax"] is None:
                tmb_panel = self.layout.panels.get(tmb_name)
                tmb_ref["ax"] = tmb_panel.ax if tmb_panel is not None and tmb_panel.plot_func is not None else None
            # self.layout.set_plot_func("tmb legend", self.plotter.plot_legend, cmap=self.tmb_cmap)
            self.layout.set_plot_func(
                ComutPanels.recurrence + label,
                self.plotter.plot_recurrence,
                snv=data.snv,
                cnv=data.cnv,
                genes=data.genes,
                cnv_cmap=self.cnv_cmap,
                categories=self.recurrence_categories,
                # max_xlim=recurrence_max_xlim,
                pad=0.01,
                invert_x=not on_left,
                label_bottom=(ComutPanels.total_recurrence_overall + label) not in self.layout.panels,
                set_joint_title=special if has_control and ComutPanels.recurrence_fold_change not in self.layout.panels_to_plot else None,
            )
            # self.layout.set_plot_func(
            #     "total recurrence" + label,
            #     self.plotter.plot_total_recurrence,
            #     total_recurrence_per_gene=data.get_total_recurrence(),
            #     pad=0.01,
            #     invert_x=not is_case,
            #     set_joint_title=special if has_control and "recurrence fold change" not in self.layout.panels_to_plot else None,
            # )
            self.layout.set_plot_func(
                ComutPanels.total_recurrence_overall + label,
                self.plotter.plot_total_recurrence_overall,
                total_recurrence_overall=data.get_total_recurrence_overall(categories=self.recurrence_categories),
                shared_x_ax=self.layout.panels.get(ComutPanels.recurrence + label).ax if (ComutPanels.recurrence + label) in self.layout.panels and self.layout.panels.get(ComutPanels.recurrence + label).plot_func is not None else None,
                pad=0.01,
                invert_x=not on_left,
                # set_joint_title=special if has_control and "recurrence fold change" not in self.layout.panels_to_plot else None,
            )
            self.layout.set_plot_func(
                ComutPanels.meta_data + label,
                self.plotter.plot_meta_data,
                meta_data=data.meta.df,
                meta_data_color=meta_data_color,
                legend_titles=data.meta.legend_titles,
                inter_heatmap_linewidth=self.layout.inter_heatmap_linewidth,
                aspect_ratio=self.layout.aspect_ratio,
                labelbottom=self.layout.show_patient_names and self._meta_data_position == "bottom",
                add_ylabel=on_left,
            )

        self.layout.set_plot_func(
            ComutPanels.model_annotation,
            self.plotter.plot_model_annotation,
            model_annotation=self.case.get_model_annotation()
        )
        self.layout.set_plot_func(
            ComutPanels.gene_names,
            self.plotter.plot_gene_names,
            genes=self.case.genes,
            ground_truth_genes=self.case.ground_truth_genes
        )
        self.layout.set_plot_func(
            ComutPanels.cytoband,
            self.plotter.plot_cytoband,
            cytobands=self.case.cnv.gistic.cytoband
        )
        self.layout.set_plot_func(
            ComutPanels.gene_meta_data,
            self.plotter.plot_gene_meta_data,
            gene_meta_data=self.gene_meta_data,
        )

        self.layout.set_plot_func(
            ComutPanels.model_significance,
            self.plotter.plot_model_significance,
            model_significance=self.model_significance
        )
        self.layout.set_plot_func(
            ComutPanels.mutational_signatures_legend,
            self.plotter.plot_legend,
            cmap=self.signatures_cmap,
            title="Mutational Signatures"
        )
        self.layout.set_plot_func(
            ComutPanels.snv_legend,
            self.plotter.plot_legend,
            cmap=self.snv_cmap,
            title="Short Nucleotide\nVariations"
        )
        self.layout.set_plot_func(
            ComutPanels.cnv_legend,
            self.plotter.plot_legend,
            cmap=self.cnv_cmap,
            names=self.cnv_names,
            title="Copy Number\nVariations"
        )
        self.layout.set_plot_func(
            ComutPanels.epi_legend,
            self.plotter.plot_epi_legend,
            alphas=self.epi_alphas,
            title=self._mark_label
        )
        self.layout.set_plot_func(
            ComutPanels.model_annotation_legend,
            self.plotter.plot_model_annotation_legend,
            names=self.model_names,
            title="Significant by Model"
        )

        for title, cmap in self.meta_cmaps_condensed.items():
            self.layout.set_plot_func(
                meta_legend(title),
                self.plotter.plot_legend,
                cmap=cmap,
                title=title,
                title_loc="bottom"
            )

        for _, panel in self.layout.panels.items():
            # some panels may have zero size and are not placed on the gridspec:
            if panel.ax is not None and panel.plot_func is not None:
                panel.plot_func(panel.ax)

        self.plotter.save_figure(fig=self.layout.fig, bbox_inches="tight")
        self.plotter.close_figure(fig=self.layout.fig)

        return self.case.genes, self.case.columns

    # --- Grid rendering ------------------------------------------------------

    def _compute_tmb_ylims(self):
        if self.joint.tmb is None:
            return 5 * 1e-1, 1.01 * 1e2
        elif SA.tmb in self.joint.tmb.columns:
            tmb_ymin = min(self.joint.tmb[SA.tmb].quantile(0.1), 5 * 1e-1)
            tmb_ymax = np.clip(self.joint.tmb[SA.tmb].max(), a_min=1.01 * 1e2, a_max=1e4)
            return tmb_ymin, tmb_ymax
        else:
            return 0, np.clip(self.joint.tmb.sum(axis=1).max(), a_min=1.01 * 1e2, a_max=1e4)

    def _make_meta_data_color(self):
        def meta_data_color(column, value):
            cmap = self.meta_cmaps[column]
            if isinstance(cmap, Palette):
                return cmap[value] if value in cmap else self.plotter.palette.get(value, self.plotter.palette.white)
            else:
                _cmap, _norm = cmap
                return _cmap(_norm(value))
        return meta_data_color

    def _fold_change_xmax(self, fc):
        at = self.amp_thresholds + self.del_thresholds
        return max(2.05, 1.2 * max(np.abs(np.min(fc["mean"][at])), np.max(fc["mean"][at])))

    def _render_grid(self):
        """Render a grid of comut sub-panels. Marginal scales/legends are global
        (Decision 5); recurrence marginals use the full cohort (Decision 4)."""
        self.layout.add_panels()

        tmb_ymin, tmb_ymax = self._compute_tmb_ylims()
        meta_data_color = self._make_meta_data_color()
        has_control = len(self.control.columns) > 0
        # Effective mirror flag (only meaningful with a control cohort).
        control_on_left = self._control_on_left and has_control
        global_fc = self.get_recurrence_fold_change_by_gene(alpha_ci=0.2, base=2) if has_control else None
        fc_xmax = self._fold_change_xmax(global_fc) if global_fc is not None else None

        def slice_df(df, index=None, columns=None):
            if df is None:
                return None
            if index is not None:
                df = df.reindex(index=index)
            if columns is not None:
                df = df.reindex(columns=columns)
            return df

        def make_plot_comut(cnv_df, snv_df, epi_df, labelbottom):
            def _p(ax):
                self.plotter.plot_cnv_heatmap(
                    ax=ax, cnv=cnv_df, cnv_cmap=self.cnv_cmap,
                    inter_heatmap_linewidth=self.layout.inter_heatmap_linewidth,
                    aspect_ratio=self.layout.aspect_ratio,
                )
                self.plotter.plot_epi_heatmap(
                    ax=ax, epi=epi_df,
                    inter_heatmap_linewidth=self.layout.inter_heatmap_linewidth,
                    aspect_ratio=self.layout.aspect_ratio,
                )
                self.plotter.plot_snv_heatmap(ax=ax, snv=snv_df, snv_cmap=self.snv_cmap)
                self.plotter.plot_heatmap_layout(ax=ax, cna=cnv_df, labelbottom=labelbottom)
            return _p

        def panel_ax(name):
            p = self.layout.panels.get(name)
            return p.ax if p is not None else None

        tmb_ref = {"ax": None}

        def set_cohort(data, grid, is_case, on_left):
            label = "" if is_case else " control"
            # The comut blocks are always placed as scaffold; only draw the
            # heatmap into them when the comutation panel is actually requested,
            # otherwise blank them so they stay invisible (see add_panels_grid).
            draw_comut = (ComutPanels.comutation + label) in self.layout.panels_to_plot
            n_rows = grid.n_rows
            for i, gene_group in enumerate(grid.gene_groups):
                genes_i = list(gene_group.members)
                for j, col_group in enumerate(grid.column_groups):
                    block_name = block(ComutPanels.comutation + label, i, j)
                    if not draw_comut:
                        self.layout.set_plot_func(block_name, self.plotter.plot_blank)
                        continue
                    cols_j = list(col_group.members)
                    cnv_df = slice_df(data.cnv.df, index=genes_i, columns=cols_j)
                    snv_df = slice_df(data.snv.df, index=genes_i, columns=cols_j)
                    epi_df = slice_df(data.epi.df, index=genes_i, columns=cols_j)
                    labelbottom = (
                        self.layout.show_patient_names
                        and (
                            ComutPanels.meta_data not in self.layout.panels_to_plot
                            or self._meta_data_position == "top"
                        )
                        and i == n_rows - 1
                    )
                    self.layout.set_plot_func(
                        block_name,
                        make_plot_comut(cnv_df, snv_df, epi_df, labelbottom),
                    )
            for j, col_group in enumerate(grid.column_groups):
                cols_j = list(col_group.members)
                if col_group.label:
                    self.layout.set_plot_func(
                        grid_col_panel(ComutPanels.column_group_label + label, j),
                        self.plotter.plot_group_label,
                        label=col_group.label,
                        is_top=True,
                        offset=j,
                    )
                if data.signatures is not None:
                    self.layout.set_plot_func(
                        grid_col_panel(ComutPanels.mutational_signatures + label, j),
                        self.plotter.plot_signatures,
                        signatures=slice_df(data.signatures, index=cols_j),
                        signatures_cmap=self.signatures_cmap,
                        add_ylabel=on_left and j == 0,
                    )
                if data.tmb is not None:
                    tmb_name = grid_col_panel(ComutPanels.tmb + label, j)
                    self.layout.set_plot_func(
                        tmb_name, self.plotter.plot_tmb,
                        tmb=slice_df(data.tmb, index=cols_j),
                        ytickpad=0,
                        shared_y_ax=tmb_ref["ax"],
                        aspect_ratio=self.layout.aspect_ratio,
                        ymin=tmb_ymin, ymax=tmb_ymax,
                        median_on_right=on_left,
                    )
                    if tmb_ref["ax"] is None:
                        tmb_ref["ax"] = panel_ax(tmb_name)
                self.layout.set_plot_func(
                    grid_col_panel(ComutPanels.cohort_label + label, j),
                    self.plotter.plot_cohort_label, label=data.name,
                )
                if data.meta is not None and not data.meta.df.empty:
                    self.layout.set_plot_func(
                        grid_col_panel(ComutPanels.meta_data + label, j),
                        self.plotter.plot_meta_data,
                        meta_data=slice_df(data.meta.df, index=cols_j),
                        meta_data_color=meta_data_color,
                        legend_titles=data.meta.legend_titles,
                        inter_heatmap_linewidth=self.layout.inter_heatmap_linewidth,
                        aspect_ratio=self.layout.aspect_ratio,
                        labelbottom=self.layout.show_patient_names and self._meta_data_position == "bottom",
                        add_ylabel=on_left and j == 0,
                    )

        cohorts = [(self.case, self.case.grid, True)]
        if has_control and self.control.grid is not None and not self.control.grid.is_trivial():
            cohorts.insert(0 if control_on_left else 1, (self.control, self.control.grid, False))
        # The left-most cohort owns the shared y-axis labels/ticks.
        for _k, (_data, _grid, _is_case) in enumerate(cohorts):
            set_cohort(_data, _grid, is_case=_is_case, on_left=_k == 0)

        gene_name = self.case.genes.name
        n_gene_rows = self.case.grid.n_rows
        last_row = n_gene_rows - 1
        # Total-recurrence summary strips (if placed) live beneath the last row's
        # recurrence columns; when present, the last row's recurrence panels must
        # not draw their own bottom axis (the strip carries it), mirroring non-grid.
        has_total_overall = grid_row_panel(ComutPanels.total_recurrence_overall, last_row) in self.layout.panels
        has_total_fc = grid_row_panel(ComutPanels.total_recurrence_fold_change, last_row) in self.layout.panels
        has_total_overall_control = grid_row_panel(control(ComutPanels.total_recurrence_overall), last_row) in self.layout.panels
        for i, gene_group in enumerate(self.case.grid.gene_groups):
            is_top_row = i == 0
            is_bottom_row = i == n_gene_rows - 1
            genes_i = pd.Index(list(gene_group.members), name=gene_name)
            if gene_group.label:
                self.layout.set_plot_func(
                    grid_row_panel(ComutPanels.gene_group_label, i),
                    self.plotter.plot_group_label,
                    label=gene_group.label,
                    is_top=False,
                    offset=i,
                )
            self.layout.set_plot_func(
                grid_row_panel(ComutPanels.model_annotation, i),
                self.plotter.plot_model_annotation,
                model_annotation=slice_df(self.case.get_model_annotation(), index=genes_i),
            )
            self.layout.set_plot_func(
                grid_row_panel(ComutPanels.gene_names, i),
                self.plotter.plot_gene_names,
                genes=genes_i,
                ground_truth_genes=self.case.ground_truth_genes,
            )
            self.layout.set_plot_func(
                grid_row_panel(ComutPanels.cytoband, i),
                self.plotter.plot_cytoband,
                cytobands=self.case.cnv.gistic.cytoband.reindex(index=genes_i),
                show_title=is_top_row,
            )
            if self.gene_meta_data is not None:
                self.layout.set_plot_func(
                    grid_row_panel(ComutPanels.gene_meta_data, i),
                    self.plotter.plot_gene_meta_data,
                    gene_meta_data=self.gene_meta_data.reindex(index=genes_i),
                    show_title=is_top_row,
                )
            self.layout.set_plot_func(
                grid_row_panel(ComutPanels.recurrence, i),
                self.plotter.plot_recurrence,
                snv=self.case.snv,
                cnv=self.case.cnv,
                genes=genes_i,
                cnv_cmap=self.cnv_cmap,
                categories=self.recurrence_categories,
                pad=0.01,
                invert_x=control_on_left,
                label_top=is_top_row,
                label_bottom=is_bottom_row and not has_total_overall,
            )
            if has_control and global_fc is not None:
                fc_i = {k: (v.reindex(index=genes_i) if isinstance(v, pd.DataFrame) else v) for k, v in global_fc.items()}
                self.layout.set_plot_func(
                    grid_row_panel(ComutPanels.recurrence_fold_change, i),
                    self.plotter.plot_recurrence_fold_change,
                    fold_change=fc_i,
                    genes=genes_i,
                    cnv_cmap=self.cnv_cmap,
                    categories=self.recurrence_categories,
                    effects=self.joint.snv.effects,
                    amp_thresholds=self.amp_thresholds,
                    del_thresholds=self.del_thresholds,
                    label_top=is_top_row,
                    label_x=is_bottom_row and not has_total_fc,
                    pad=0.01,
                    xmax=fc_xmax,
                    case_on_left=not control_on_left,
                    case_label=self._cohort_tag,
                    control_label=self._control_cohort_tag,
                )
            # Control recurrence (full control cohort), mirroring the case
            # recurrence marginal, when requested and control data is present.
            if has_control:
                self.layout.set_plot_func(
                    grid_row_panel(control(ComutPanels.recurrence), i),
                    self.plotter.plot_recurrence,
                    snv=self.control.snv,
                    cnv=self.control.cnv,
                    genes=genes_i,
                    cnv_cmap=self.cnv_cmap,
                    categories=self.recurrence_categories,
                    pad=0.01,
                    invert_x=not control_on_left,
                    label_top=is_top_row,
                    label_bottom=is_bottom_row and not has_total_overall_control,
                )

        # --- total recurrence summary strips at the very bottom (same as non-grid) ---
        def _row_ax(base):
            p = self.layout.panels.get(grid_row_panel(base, last_row))
            return p.ax if p is not None else None

        rec_ax = _row_ax(ComutPanels.recurrence)
        if has_total_overall and rec_ax is not None:
            self.layout.set_plot_func(
                grid_row_panel(ComutPanels.total_recurrence_overall, last_row),
                self.plotter.plot_total_recurrence_overall,
                total_recurrence_overall=self.case.get_total_recurrence_overall(categories=self.recurrence_categories),
                shared_x_ax=rec_ax, pad=0.01, invert_x=control_on_left,
            )
        if has_control:
            fc_ax = _row_ax(ComutPanels.recurrence_fold_change)
            if has_total_fc and fc_ax is not None:
                self.layout.set_plot_func(
                    grid_row_panel(ComutPanels.total_recurrence_fold_change, last_row),
                    self.plotter.plot_total_recurrence_fold_change,
                    fold_change=self.get_total_recurrence_fold_change(alpha_ci=0.2, base=2),
                    shared_x_ax=fc_ax, pad=0.01,
                )
            rec_control_ax = _row_ax(control(ComutPanels.recurrence))
            if has_total_overall_control and rec_control_ax is not None:
                self.layout.set_plot_func(
                    grid_row_panel(control(ComutPanels.total_recurrence_overall), last_row),
                    self.plotter.plot_total_recurrence_overall,
                    total_recurrence_overall=self.control.get_total_recurrence_overall(categories=self.recurrence_categories),
                    shared_x_ax=rec_control_ax, pad=0.01, invert_x=not control_on_left,
                )

        self._set_grid_legends()

        for _, panel in self.layout.panels.items():
            if panel.ax is not None and panel.plot_func is not None:
                panel.plot_func(panel.ax)

        self.plotter.save_figure(fig=self.layout.fig, bbox_inches="tight")
        self.plotter.close_figure(fig=self.layout.fig)
        return self.case.genes, self.case.columns

    def _set_grid_legends(self):
        self.layout.set_plot_func(
            ComutPanels.mutational_signatures_legend, self.plotter.plot_legend,
            cmap=self.signatures_cmap, title="Mutational Signatures",
        )
        self.layout.set_plot_func(
            ComutPanels.snv_legend, self.plotter.plot_legend,
            cmap=self.snv_cmap, title="Short Nucleotide\nVariations",
        )
        self.layout.set_plot_func(
            ComutPanels.cnv_legend, self.plotter.plot_legend,
            cmap=self.cnv_cmap, names=self.cnv_names, title="Copy Number\nVariations",
        )
        self.layout.set_plot_func(
            ComutPanels.epi_legend, self.plotter.plot_epi_legend,
            alphas=self.epi_alphas, title=self._mark_label,
        )
        self.layout.set_plot_func(
            ComutPanels.model_annotation_legend, self.plotter.plot_model_annotation_legend,
            names=self.model_names, title="Significant by Model",
        )
        for title, cmap in self.meta_cmaps_condensed.items():
            self.layout.set_plot_func(
                meta_legend(title), self.plotter.plot_legend,
                cmap=cmap, title=title, title_loc="bottom",
            )
