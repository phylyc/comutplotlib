from copy import deepcopy
import numpy as np
import pandas as pd

from comutplotlib.comut_data import ComutData
from comutplotlib.comut_layout import ComutLayout
from comutplotlib.comut_plotter import ComutPlotter
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

        model_significances: list[str] = None,
        model_names: list[str] = (),

        sif: list[str] = None,
        meta_data_rows: list[str] = (),
        meta_data_rows_per_sample: list[str] = (),

        control_maf: list[str] = None,
        control_seg: list[str] = None,
        control_gistic: list[str] = None,
        control_sif: list[str] = None,
        control_mutsig: list[str] = None,

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
        snv_recurrence_threshold: int = 5,

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
            total_recurrence_threshold=0,
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

        n_genes, n_samples_case, n_meta_case = self.case.get_dimensions()
        _, n_samples_control, n_meta_control = self.control.get_dimensions()

        n_meta = max(n_meta_case, n_meta_control)

        self.joint = deepcopy(self.case)
        self.joint.gistic = join_gistics([self.case.gistic, self.control.gistic])
        self.joint.seg = join_segs([self.case.seg, self.control.seg])
        self.joint.maf = join_mafs([self.case.maf, self.control.maf])
        self.joint.sif = join_sifs([self.case.sif, self.control.sif])
        mutsigs = [s for s in [self.case.mutsig, self.control.mutsig] if s is not None]
        self.joint.mutsig = pd.concat(mutsigs) if len(mutsigs) else None
        self.joint.preprocess()

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

        self.case.meta.reindex(index=self.joint.meta.rows)
        self.control.meta.reindex(index=self.joint.meta.rows)

        if len(self.control.columns) > 0:
            self.case.columns = self.case.columns[::-1]
            self.case.reindex_data()

        def remove(col):
            if col in panels_to_plot:
                panels_to_plot.remove(col)

        if self.case.tmb is None:
            remove("tmb")
            remove("tmb legend")
        if self.case.mutsig is None:
            remove("mutational signatures")
            remove("mutational signatures legend")

        self.layout = ComutLayout(
            panels_to_plot=panels_to_plot,
            max_xfigsize=max_xfigsize,
            max_xfigsize_scale=max_xfigsize_scale,
            n_genes=n_genes,
            n_samples=n_samples_case,
            n_samples_control=n_samples_control,
            n_meta=n_meta,
            label_columns=label_columns,
            tmb_cmap=self.tmb_cmap,
            snv_cmap=self.snv_cmap,
            cnv_cmap=self.cnv_cmap,
            mutsig_cmap=self.mutsig_cmap,
            meta_cmaps=self.meta_cmaps_condensed,
        )

        self.case.save(out_dir=self.plotter.out_dir, name=self.plotter.file_name + ".case")
        self.control.save(out_dir=self.plotter.out_dir, name=self.plotter.file_name + ".control")

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

        has_control = len(self.control.columns) > 0
        for data, label, is_case, special in zip([self.case, self.control], ["", " control"], [True, False], [has_control, False]):
            if len(data.columns) == 0:
                continue

            self.layout.set_plot_func(
                "comutation" + label,
                plot_comut_gen(data)
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
                shared_y_ax=self.layout.panels.get("tmb").ax if self.layout.panels.get("tmb").plot_func is not None else None,
                aspect_ratio=self.layout.aspect_ratio,
                ymin=tmb_ymin,
                ymax=tmb_ymax
            )
            # self.layout.set_plot_func("tmb legend", self.plotter.plot_legend, cmap=self.tmb_cmap)
            self.layout.set_plot_func(
                "recurrence" + label,
                self.plotter.plot_recurrence,
                mut_protein_data=data.snv.get_recurrence(index=data.genes),
                # mut_prevalence_counter=data.snv.get_effect_prevalence().reindex(index=data.genes),
                cna=data.cnv.df,
                cna_counter=data.cnv.get_prevalence(),
                # num_effects=sdata.snv.num_effects,
                genes=data.genes,
                columns=data.columns,
                cnv_cmap=self.cnv_cmap,
                amp_thresholds=[data.high_amp_threshold, data.mid_amp_threshold, data.low_amp_threshold],
                del_thresholds=[data.high_del_threshold, data.mid_del_threshold, data.low_del_threshold],
                snv_recurrence_threshold=data.snv_recurrence_threshold
            )
            self.layout.set_plot_func(
                "total recurrence" + label,
                self.plotter.plot_total_recurrence,
                total_recurrence_per_gene=data.get_total_recurrence(),
                pad=0.01,
                invert_x=special,
                set_joint_title=special if has_control else None,
            )
            self.layout.set_plot_func(
                "total recurrence overall" + label,
                self.plotter.plot_total_recurrence_overall,
                total_recurrence_overall=data.get_total_recurrence_overall(),
                pad=0.01,
                invert_x=special,
                case_control_spacing=special if has_control else None,
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

        return self.case.genes, self.case.columns
