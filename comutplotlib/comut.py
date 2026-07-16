from __future__ import annotations

from collections.abc import Sequence
from copy import deepcopy
import numpy as np
import pandas as pd

from comutplotlib.comut_data import ComutData
from comutplotlib.comut_layout import ComutLayout
from comutplotlib.comut_panels import ComutPanels, meta_legend
from comutplotlib.comut_plotter import ComutPlotter
from comutplotlib.functional_effect import sort_functional_effects
from comutplotlib.mathutils import fold_change
from comutplotlib.mutation_annotation import MutationAnnotation as MutA
from comutplotlib.palette import Palette
from comutplotlib.sample_annotation import SampleAnnotation as SA

from comutplotlib.gistic import join_gistics
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

        signatures: list[str] | None = None,

        cohort_label: str | None = None,

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
        control_sif: list[str] | None = None,
        control_signatures: list[str] | None = None,
        control_cohort_label: str | None = None,

        drop_empty_columns: bool = False,

        by: str = MutA.patient,

        column_order: tuple[str] | None = None,
        control_column_order: tuple[str] | None = None,
        index_order: tuple[str] | None = None,
        column_sort_by: tuple[str] = ("COMUT",),
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
        self._maf = maf
        self._maf_pool_as = maf_pool_as
        self._seg = seg
        self._gistic = gistic
        self._signatures = signatures
        self._sif = sif
        self._meta_data_rows = meta_data_rows
        self._meta_data_rows_per_sample = meta_data_rows_per_sample
        self._gene_meta_data = gene_meta_data
        self._gene_meta_columns = gene_meta_columns
        self._model_significances = model_significances
        self._model_names = model_names
        self._control_maf = control_maf
        self._control_seg = control_seg
        self._control_gistic = control_gistic
        self._control_sif = control_sif
        self._control_signatures = control_signatures
        self._control_cohort_label = control_cohort_label
        self._drop_empty_columns = drop_empty_columns
        self._by = by
        self._column_order = column_order
        self._control_column_order = control_column_order
        self._index_order = index_order
        self._column_sort_by = column_sort_by
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

        self.load()
        self.preprocess()
        self.build_layout()

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
            signatures_paths=self._signatures,
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
        )
        self.case.preprocess()

        self.control = ComutData(
            cohort_name=self._control_cohort_label,
            maf_paths=self._control_maf,
            maf_pool_as=self._maf_pool_as,
            seg_paths=self._control_seg,
            gistic_paths=self._control_gistic,
            signatures_paths=self._control_signatures,
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
        recurrence_fold_change = self.get_recurrence_fold_change_by_gene(alpha_ci=0.1, base=2)
        categories = {
            "snv": self.joint.snv.effects,
            "amp": self.amp_thresholds,
            "del": self.del_thresholds
        }
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

        self.tmb_cmap = self.plotter.palette.get_tmb_cmap(self.joint.tmb)
        self.snv_cmap = self.plotter.palette.get_snv_cmap(self.joint.snv)
        self.cnv_cmap, self.cnv_names = self.plotter.palette.get_cnv_cmap(self.joint.cnv)
        self.signatures_cmap = self.plotter.palette.get_signatures_cmap(self.joint.signatures)
        self.meta_cmaps = self.plotter.palette.get_meta_cmaps(self.joint.meta)
        if self._palette is not None:
            for col, pal in self._palette["local"].items():
                if col in self.meta_cmaps:
                    self.meta_cmaps[col] |= pal
                else:
                    self.meta_cmaps[col] = pal
        self.meta_cmaps_condensed = self.plotter.palette.condense(self.meta_cmaps)

        if len(self.control.columns) > 0:
            self.case.columns = self.case.columns[::-1]
            self.case.reindex_data()

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

    def build_layout(self):
        """Compute panel dimensions and instantiate the figure layout."""
        n_meta_genes = self.gene_meta_data.shape[1] if self.gene_meta_data is not None else 0

        n_genes, n_samples_case, n_meta_case = self.case.get_dimensions()
        _, n_samples_control, n_meta_control = self.control.get_dimensions()

        n_meta = max(n_meta_case, n_meta_control)

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
            tmb_cmap=self.tmb_cmap,
            snv_cmap=self.snv_cmap,
            cnv_cmap=self.cnv_cmap,
            mutsig_cmap=self.signatures_cmap,
            meta_cmaps=self.meta_cmaps_condensed,
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

    def make_comut(self):
        """Backwards-compatible alias for :meth:`render`."""
        return self.render()

    def render(self):
        # Persist the finalized column/gene ordering next to the figure output.
        self.case.save(out_dir=self.plotter.out_dir, name=self.plotter.file_name + ".case")
        self.control.save(out_dir=self.plotter.out_dir, name=self.plotter.file_name + ".control")

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
                self.plotter.plot_snv_heatmap(
                    ax=ax,
                    snv=data.snv.df,
                    snv_cmap=self.snv_cmap
                )
                self.plotter.plot_heatmap_layout(
                    ax=ax,
                    cna=data.cnv.df,
                    labelbottom=self.layout.show_patient_names and data.meta is None
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
            )
            self.layout.set_plot_func(
                ComutPanels.total_recurrence_fold_change,
                self.plotter.plot_total_recurrence_fold_change,
                fold_change=self.get_total_recurrence_fold_change(alpha_ci=0.2, base=2),
                shared_x_ax=self.layout.panels.get(ComutPanels.recurrence_fold_change).ax if ComutPanels.recurrence_fold_change in self.layout.panels and self.layout.panels.get(ComutPanels.recurrence_fold_change).plot_func is not None else None,
                pad=0.01,
            )

        for data, label, is_case, special in zip([self.case, self.control], ["", " control"], [True, False], [has_control, False]):
            if len(data.columns) == 0:
                continue

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
                add_ylabel=is_case,
            )
            self.layout.set_plot_func(
                ComutPanels.tmb + label,
                self.plotter.plot_tmb,
                tmb=data.tmb,
                ytickpad=0,
                fontsize=6,
                shared_y_ax=self.layout.panels.get(ComutPanels.tmb).ax if ComutPanels.tmb in self.layout.panels and self.layout.panels.get(ComutPanels.tmb).plot_func is not None else None,
                aspect_ratio=self.layout.aspect_ratio,
                ymin=tmb_ymin,
                ymax=tmb_ymax
            )
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
                invert_x=not is_case,
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
                invert_x=not is_case,
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
                labelbottom=self.layout.show_patient_names,
                add_ylabel=is_case,
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
            fontsize=6,
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
            title="Short Nucleotide Variations"
        )
        self.layout.set_plot_func(
            ComutPanels.cnv_legend,
            self.plotter.plot_legend,
            cmap=self.cnv_cmap,
            names=self.cnv_names,
            title="Copy Number Variations"
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
