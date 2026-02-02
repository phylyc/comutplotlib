from __future__ import annotations

from copy import deepcopy
import numpy as np
import pandas as pd
import scipy.stats as st

from comutplotlib.comut_data import ComutData
from comutplotlib.comut_layout import ComutLayout
from comutplotlib.comut_plotter import ComutPlotter
from comutplotlib.functional_effect import sort_functional_effects
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

        maf: list[str] = None,
        maf_pool_as: dict | None = None,

        seg: list[str] = None,
        gistic: list[str] = None,

        mutsig: list[str] = None,

        cohort_label: str = None,

        model_significances: list[str] = None,
        model_names: list[str] = (),

        sif: list[str] = None,
        meta_data_rows: list[str] = (),
        meta_data_rows_per_sample: list[str] = (),

        gene_meta_data: str = None,
        gene_meta_columns: list[str] = None,

        control_maf: list[str] = None,
        control_seg: list[str] = None,
        control_gistic: list[str] = None,
        control_sif: list[str] = None,
        control_mutsig: list[str] = None,
        control_cohort_label: str = None,

        drop_empty_columns: bool = False,

        by: str = MutA.patient,

        column_order: tuple[str] = None,
        control_column_order: tuple[str] = None,
        index_order: tuple[str] = None,
        column_sort_by: tuple[str] = ("COMUT",),

        interesting_gene: str = None,
        interesting_gene_comut_percent_threshold: float = None,
        interesting_genes: set = None,
        snv_interesting_genes: set = None,
        cnv_interesting_genes: set = None,
        total_recurrence_threshold: float = None,
        scale_recurrence: bool = False,
        recurrence_categories: dict[str, list[str] | dict[str, list[str]]] = {"global": ["snv", "amp", "del"]},
        snv_recurrence_threshold: int = 5,
        min_fold_change: float = None,
        max_fold_change: float = None,
        collapse_cytobands: bool = False,
        collapse_cytobands_range: int = 0,

        low_amp_threshold: int | float = 1,
        mid_amp_threshold: int | float = 1.5,
        high_amp_threshold: int | float = 2,
        baseline: int | float = 0,
        low_del_threshold: int | float = -1,
        mid_del_threshold: int | float = -1.5,
        high_del_threshold: int | float = -2,
        show_low_level_cnvs: bool = False,

        panels_to_plot: list = None,
        palette: dict[str, dict[str, tuple[float]], dict[str, dict[str, tuple[float]]]] = None,
        ground_truth_genes: dict[str, list[str]] = None,  # todo: refactor as a palette class
        max_xfigsize: int = None,
        max_xfigsize_scale: float = 1,
        label_columns: bool = False,
        **kwargs
    ):
        self.plotter = ComutPlotter(
            output=output,
            extra_palette=palette["global"] if palette is not None else None
        )

        self.case = ComutData(
            cohort_name=cohort_label,
            maf_paths=maf,
            maf_pool_as=maf_pool_as,
            seg_paths=seg,
            gistic_paths=gistic,
            mutsig_paths=mutsig,
            sif_paths=sif,
            meta_data_rows=meta_data_rows,
            meta_data_rows_per_sample=meta_data_rows_per_sample,
            drop_empty_columns=drop_empty_columns,
            by=by,
            column_order=column_order,
            index_order=index_order,
            column_sort_by=column_sort_by,
            interesting_gene=interesting_gene,
            interesting_gene_comut_percent_threshold=interesting_gene_comut_percent_threshold,
            interesting_genes=interesting_genes,
            ground_truth_genes=ground_truth_genes,
            snv_interesting_genes=snv_interesting_genes,
            cnv_interesting_genes=cnv_interesting_genes,
            total_recurrence_threshold=total_recurrence_threshold,
            snv_recurrence_threshold=snv_recurrence_threshold,
            low_amp_threshold=low_amp_threshold,
            mid_amp_threshold=mid_amp_threshold,
            high_amp_threshold=high_amp_threshold,
            baseline=baseline,
            low_del_threshold=low_del_threshold,
            mid_del_threshold=mid_del_threshold,
            high_del_threshold=high_del_threshold,
            show_low_level_cnvs=show_low_level_cnvs,
        )
        self.case.preprocess()

        self.control = ComutData(
            cohort_name=control_cohort_label,
            maf_paths=control_maf,
            maf_pool_as=maf_pool_as,
            seg_paths=control_seg,
            gistic_paths=control_gistic,
            mutsig_paths=control_mutsig,
            sif_paths=control_sif,
            meta_data_rows=meta_data_rows,
            meta_data_rows_per_sample=meta_data_rows_per_sample,
            drop_empty_columns=drop_empty_columns,
            by=by,
            column_order=control_column_order,
            index_order=self.case.genes,
            column_sort_by=column_sort_by,
            interesting_gene=self.case.interesting_gene,
            interesting_gene_comut_percent_threshold=interesting_gene_comut_percent_threshold,
            interesting_genes=self.case.interesting_genes,
            ground_truth_genes=self.case.ground_truth_genes,
            snv_interesting_genes=self.case.snv_interesting_genes,
            cnv_interesting_genes=self.case.cnv_interesting_genes,
            snv_recurrence_threshold=snv_recurrence_threshold,
            low_amp_threshold=low_amp_threshold,
            mid_amp_threshold=mid_amp_threshold,
            high_amp_threshold=high_amp_threshold,
            baseline=baseline,
            low_del_threshold=low_del_threshold,
            mid_del_threshold=mid_del_threshold,
            high_del_threshold=high_del_threshold,
            show_low_level_cnvs=show_low_level_cnvs,
        )
        self.control.preprocess()

        self.model_significance = pd.DataFrame.from_dict(
            {
                name: pd.read_csv(path_to_file, index_col=0, sep="\t")
                for name, path_to_file in zip(model_names, model_significances)
            }
        ) if model_significances is not None else None
        self.model_names = model_names

        self.joint = deepcopy(self.case)
        self.joint.gistic = join_gistics([self.case.gistic, self.control.gistic])
        self.joint.seg = join_segs([self.case.seg, self.control.seg])
        self.joint.maf = join_mafs([self.case.maf, self.control.maf])
        self.joint.sif = join_sifs([self.case.sif, self.control.sif])
        mutsigs = [s for s in [self.case.mutsig, self.control.mutsig] if s is not None]
        self.joint.mutsig = pd.concat(mutsigs) if len(mutsigs) else None
        self.joint.preprocess()

        self.case.meta.reindex(index=self.joint.meta.rows)
        self.control.meta.reindex(index=self.joint.meta.rows)

        self.amp_thresholds = [self.joint.high_amp_threshold, self.joint.mid_amp_threshold, self.joint.low_amp_threshold]
        self.del_thresholds = [self.joint.high_del_threshold, self.joint.mid_del_threshold, self.joint.low_del_threshold]

        self.scale_recurrence = scale_recurrence
        self.recurrence_categories = self.get_recurrence_categories_by_gene(recurrence_categories)
        recurrence_fold_change = self.get_recurrence_fold_change_by_gene(fdr=0.1, base=2)
        categories = {
            "snv": self.joint.snv.effects,
            "amp": self.amp_thresholds,
            "del": self.del_thresholds
        }
        good_genes = []
        for g, row in recurrence_fold_change["mean"].iterrows():
            relevant_columns = [c for cat in self.recurrence_categories[g] for c in categories[cat]]
            is_good = True
            if min_fold_change is not None:
                is_good &= (row[relevant_columns] > np.log(min_fold_change) / np.log(2)).any()
            if max_fold_change is not None:
                is_good &= (row[relevant_columns] < np.log(max_fold_change) / np.log(2)).any()
            if is_good:
                good_genes.append(g)

        for data in [self.case, self.control, self.joint]:
            data.genes = good_genes
            data.reindex_data()

        self.recurrence_categories = self.get_recurrence_categories_by_gene(recurrence_categories)

        self.tmb_cmap = self.plotter.palette.get_tmb_cmap(self.joint.tmb)
        self.snv_cmap = self.plotter.palette.get_snv_cmap(self.joint.snv)
        self.cnv_cmap, self.cnv_names = self.plotter.palette.get_cnv_cmap(self.joint.cnv)
        self.mutsig_cmap = self.plotter.palette.get_mutsig_cmap(self.joint.mutsig)
        self.meta_cmaps = self.plotter.palette.get_meta_cmaps(self.joint.meta)
        if palette is not None:
            for col, pal in palette["local"].items():
                if col in self.meta_cmaps:
                    self.meta_cmaps[col] |= pal
                else:
                    self.meta_cmaps[col] = pal
        self.meta_cmaps_condensed = self.plotter.palette.condense(self.meta_cmaps)

        if len(self.control.columns) > 0:
            self.case.columns = self.case.columns[::-1]
            self.case.reindex_data()

        def remove(col):
            if col in panels_to_plot:
                panels_to_plot.remove(col)

        if "tmb" in self.tmb_cmap.keys():
            remove("tmb legend")
        if self.joint.tmb is None:
            remove("tmb")
            remove("tmb legend")
        if self.joint.mutsig is None:
            remove("mutational signatures")
            remove("mutational signatures legend")

        if gene_meta_data is not None:
            self.gene_meta_data = (
                pd.read_csv(gene_meta_data, sep="\t")
                    .set_index("Gene")
                    .reindex(index=self.case.genes, columns=gene_meta_columns)
                    .replace(0, np.nan)
                    .dropna(axis=1, how="all")
                    .fillna(0)
            )
            self.gene_meta_data = self.gene_meta_data[self.gene_meta_data.sum().sort_values().index]
            n_meta_genes = self.gene_meta_data.shape[1]
        else:
            self.gene_meta_data = None
            n_meta_genes = 0

        n_genes, n_samples_case, n_meta_case = self.case.get_dimensions()
        _, n_samples_control, n_meta_control = self.control.get_dimensions()

        n_meta = max(n_meta_case, n_meta_control)

        self.layout = ComutLayout(
            panels_to_plot=panels_to_plot,
            max_xfigsize=max_xfigsize,
            max_xfigsize_scale=max_xfigsize_scale,
            n_genes=n_genes,
            n_samples=n_samples_case,
            n_samples_control=n_samples_control,
            n_meta=n_meta,
            n_meta_genes=n_meta_genes,
            label_columns=label_columns,
            tmb_cmap=self.tmb_cmap,
            snv_cmap=self.snv_cmap,
            cnv_cmap=self.cnv_cmap,
            mutsig_cmap=self.mutsig_cmap,
            meta_cmaps=self.meta_cmaps_condensed,
        )

        self.case.save(out_dir=self.plotter.out_dir, name=self.plotter.file_name + ".case")
        self.control.save(out_dir=self.plotter.out_dir, name=self.plotter.file_name + ".control")

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

    def get_recurrence_fold_change_by_gene(self, fdr=0.1, base=2):
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

        p_case = (case + 0.5) / (n_case + 0.5)
        p_control = (control + 0.5) / (n_control + 0.5)
        mean = np.log(p_case / p_control) / np.log(base)
        var = ((1 / p_case - 1) / n_case + (1 / p_control - 1) / n_control + 1e-12) / np.log(base) ** 2

        lo = pd.DataFrame(st.norm.ppf(fdr, loc=mean, scale=np.sqrt(var)), index=mean.index, columns=mean.columns)
        hi = pd.DataFrame(st.norm.ppf(1 - fdr, loc=mean, scale=np.sqrt(var)), index=mean.index, columns=mean.columns)

        is_present = (case > 0) & (control > 0) & ~case.isna() & ~control.isna()

        return {
            "mean": mean.mask(~is_present),
            "lo": lo.mask(~is_present),
            "hi": hi.mask(~is_present),
            "fdr": fdr,
            "base": base
        }

    def get_total_recurrence_fold_change(self, fdr=0.1, base=2):
        case, n_case = self.case.get_total_recurrence_overall(categories=self.recurrence_categories)
        control, n_control = self.control.get_total_recurrence_overall(categories=self.recurrence_categories)

        p_case = (case + 0.5) / (n_case + 0.5)
        p_control = (control + 0.5) / (n_control + 0.5)
        mean = np.log(p_case / p_control) / np.log(base)
        var = ((1 / p_case - 1) / n_case + (1 / p_control - 1) / n_control + 1e-12) / np.log(base) ** 2

        lo = pd.Series(st.norm.ppf(fdr, loc=mean, scale=np.sqrt(var)), index=mean.index)
        hi = pd.Series(st.norm.ppf(1 - fdr, loc=mean, scale=np.sqrt(var)), index=mean.index)

        is_present = (case > 0) & (control > 0) | case.isna() | control.isna()
        return {
            "mean": mean.mask(~is_present),
            "lo": lo.mask(~is_present),
            "hi": hi.mask(~is_present),
            "fdr": fdr,
            "base": base
        }

    def make_comut(self):
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

        if SA.tmb in self.joint.tmb.columns:
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
                "recurrence fold change",
                self.plotter.plot_recurrence_fold_change,
                fold_change=self.get_recurrence_fold_change_by_gene(fdr=0.1, base=2),
                genes=self.case.genes,
                cnv_cmap=self.cnv_cmap,
                categories=self.recurrence_categories,
                effects=self.joint.snv.effects,
                amp_thresholds=self.amp_thresholds,
                del_thresholds=self.del_thresholds,
                label_x="total recurrence fold change" not in self.layout.panels,
                pad=0.01,
            )
            self.layout.set_plot_func(
                "total recurrence fold change",
                self.plotter.plot_total_recurrence_fold_change,
                fold_change=self.get_total_recurrence_fold_change(fdr=0.1, base=2),
                shared_x_ax=self.layout.panels.get("recurrence fold change").ax if "recurrence fold change" in self.layout.panels and self.layout.panels.get("recurrence fold change").plot_func is not None else None,
                pad=0.01,
            )

        for data, label, is_case, special in zip([self.case, self.control], ["", " control"], [True, False], [has_control, False]):
            if len(data.columns) == 0:
                continue

            self.layout.set_plot_func(
                "comutation" + label,
                plot_comut_gen(data)
            )

            self.layout.set_plot_func(
                "cohort label" + label,
                self.plotter.plot_cohort_label,
                label=data.name
            )

            # self.layout.set_plot_func("coverage", self.plotter.plot_coverage)
            self.layout.set_plot_func(
                "mutational signatures" + label,
                self.plotter.plot_mutsig,
                mutsig=data.mutsig,
                mutsig_cmap=self.mutsig_cmap,
                add_ylabel=is_case,
            )
            self.layout.set_plot_func(
                "tmb" + label,
                self.plotter.plot_tmb,
                tmb=data.tmb,
                ytickpad=0,
                fontsize=6,
                shared_y_ax=self.layout.panels.get("tmb").ax if "tmb" in self.layout.panels and self.layout.panels.get("tmb").plot_func is not None else None,
                aspect_ratio=self.layout.aspect_ratio,
                ymin=tmb_ymin,
                ymax=tmb_ymax
            )
            # self.layout.set_plot_func("tmb legend", self.plotter.plot_legend, cmap=self.tmb_cmap)
            self.layout.set_plot_func(
                "recurrence" + label,
                self.plotter.plot_recurrence,
                snv=data.snv,
                cnv=data.cnv,
                genes=data.genes,
                cnv_cmap=self.cnv_cmap,
                categories=self.recurrence_categories,
                # max_xlim=recurrence_max_xlim,
                pad=0.01,
                invert_x=not is_case,
                label_bottom=("total recurrence overall" + label) not in self.layout.panels,
                set_joint_title=special if has_control and "recurrence fold change" not in self.layout.panels_to_plot else None,
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
                "total recurrence overall" + label,
                self.plotter.plot_total_recurrence_overall,
                total_recurrence_overall=data.get_total_recurrence_overall(categories=self.recurrence_categories),
                shared_x_ax=self.layout.panels.get("recurrence" + label).ax if ("recurrence" + label) in self.layout.panels and self.layout.panels.get("recurrence" + label).plot_func is not None else None,
                pad=0.01,
                invert_x=not is_case,
                # set_joint_title=special if has_control and "recurrence fold change" not in self.layout.panels_to_plot else None,
            )
            self.layout.set_plot_func(
                "meta data" + label,
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
            "model annotation",
            self.plotter.plot_model_annotation,
            model_annotation=self.case.get_model_annotation()
        )
        self.layout.set_plot_func(
            "gene names",
            self.plotter.plot_gene_names,
            genes=self.case.genes,
            ground_truth_genes=self.case.ground_truth_genes
        )
        self.layout.set_plot_func(
            "cytoband",
            self.plotter.plot_cytoband,
            cytobands=self.case.cnv.gistic.cytoband
        )
        self.layout.set_plot_func(
            "gene meta data",
            self.plotter.plot_gene_meta_data,
            gene_meta_data=self.gene_meta_data,
            fontsize=6,
        )

        self.layout.set_plot_func(
            "model significance",
            self.plotter.plot_model_significance,
            model_significance=self.model_significance
        )
        self.layout.set_plot_func(
            "mutational signatures legend",
            self.plotter.plot_legend,
            cmap=self.mutsig_cmap,
            title="Mutational Signatures"
        )
        self.layout.set_plot_func(
            "snv legend",
            self.plotter.plot_legend,
            cmap=self.snv_cmap,
            title="Short Nucleotide Variations"
        )
        self.layout.set_plot_func(
            "cnv legend",
            self.plotter.plot_legend,
            cmap=self.cnv_cmap,
            names=self.cnv_names,
            title="Copy Number Variations"
        )
        self.layout.set_plot_func(
            "model annotation legend",
            self.plotter.plot_model_annotation_legend,
            names=self.model_names,
            title="Significant by Model"
        )

        for title, cmap in self.meta_cmaps_condensed.items():
            self.layout.set_plot_func(
                f"meta data legend {title}",
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
