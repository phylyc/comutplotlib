import warnings

from comutplotlib.comut_data import ComutData
from comutplotlib.comut_layout import ComutLayout
from comutplotlib.comut_plotter import ComutPlotter
from comutplotlib.mutation_annotation import MutationAnnotation as MutA
from comutplotlib.palette import Palette


class Comut(object):
    """ Caption:
    This comutation plot (central panel) visualizes the mutation landscape across a patient cohort. Rows represent genes, and columns correspond to patients. Each cell indicates a gene’s mutation status in a patient: rectangles denote copy-number variations (CNVs), and ellipses indicate short nucleotide variations (SNVs), with multiple SNVs shown as subdivided wedges. Colors encode mutation types and functional effects.

    The top panel displays tumor mutation burden (TMB) per patient, with high TMB (≥10/Mb) highlighted in red. The mutational signature panel shows the relative fraction of exposures to different mutational signatures for each patient or sample. The left panel summarizes mutation recurrence, showing SNV and CNV frequencies per gene, annotated with recurrent protein alterations. It also reports the percentage of patients with high-level CNVs, supplemented by low-level CNVs in brackets, along with the total percentage of patients carrying either an SNV or a high-level CNV in the gene.

    The bottom panel presents patient- and sample-level metadata. For patients with multiple samples, metadata cells are subdivided accordingly.
    """

    def __init__(
            self,
            output: str = "./comut.pdf",

            maf: list[str] = (),
            maf_pool_as: dict | None = None,

            seg: list[str] = (),
            gistic: list[str] = (),

            mutsig: list[str] = None,

            model_significances: list[str] = (),
            model_names: list[str] = (),

            sif: list[str] = (),
            meta_data_rows: list[str] = (),
            meta_data_rows_per_sample: list[str] = (),

            drop_empty_columns: bool = False,

            by: str = MutA.patient,

            column_order: tuple[str] = None,
            index_order: tuple[str] = None,

            interesting_gene: str = None,
            interesting_gene_comut_percent_threshold: float = None,
            interesting_genes: set = None,
            snv_interesting_genes: set = None,
            cnv_interesting_genes: set = None,
            total_recurrence_threshold: float = None,
            snv_recurrence_threshold: int = 5,

            low_amp_threshold: int = 1,
            high_amp_threshold: int = 2,
            low_del_threshold: int = -1,
            high_del_threshold: int = -2,
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
        self.data = ComutData(
            maf_paths=maf,
            maf_pool_as=maf_pool_as,
            seg_paths=seg,
            gistic_paths=gistic,
            mutsig_paths=mutsig,
            model_significances=model_significances,
            model_names=model_names,
            sif_paths=sif,
            meta_data_rows=meta_data_rows,
            meta_data_rows_per_sample=meta_data_rows_per_sample,
            drop_empty_columns=drop_empty_columns,
            by=by,
            column_order=column_order,
            index_order=index_order,
            interesting_gene=interesting_gene,
            interesting_gene_comut_percent_threshold=interesting_gene_comut_percent_threshold,
            interesting_genes=interesting_genes,
            ground_truth_genes=ground_truth_genes,
            snv_interesting_genes=snv_interesting_genes,
            cnv_interesting_genes=cnv_interesting_genes,
            total_recurrence_threshold=total_recurrence_threshold,
            snv_recurrence_threshold=snv_recurrence_threshold,
            low_amp_threshold=low_amp_threshold,
            high_amp_threshold=high_amp_threshold,
            low_del_threshold=low_del_threshold,
            high_del_threshold=high_del_threshold,
            show_low_level_cnvs=show_low_level_cnvs,
        )
        self.data.preprocess()
        self.data.save(out_dir=self.plotter.out_dir, name=self.plotter.file_name)

        n_genes, n_samples, n_meta = self.data.get_dimensions()

        self.tmb_cmap = self.plotter.palette.get_tmb_cmap(self.data)
        self.snv_cmap = self.plotter.palette.get_snv_cmap(self.data)
        self.cnv_cmap, self.cnv_names = self.plotter.palette.get_cnv_cmap(self.data)
        self.mutsig_cmap = self.plotter.palette.get_mutsig_cmap(self.data)
        self.meta_cmaps = self.plotter.palette.get_meta_cmaps(self.data)
        if palette is not None:
            for col, pal in palette["local"].items():
                if col in self.meta_cmaps:
                    self.meta_cmaps[col] |= pal
                else:
                    self.meta_cmaps[col] = pal
        self.meta_cmaps_condensed = self.plotter.palette.condense(self.meta_cmaps)

        def remove(col):
            if col in panels_to_plot:
                panels_to_plot.remove(col)

        if self.data.tmb is None:
            remove("tmb")
            remove("tmb legend")
        if self.data.mutsig is None:
            remove("mutational signatures")
            remove("mutsig legend")

        self.layout = ComutLayout(
            panels_to_plot=panels_to_plot,
            max_xfigsize=max_xfigsize,
            max_xfigsize_scale=max_xfigsize_scale,
            n_genes=n_genes,
            n_samples=n_samples,
            n_meta=n_meta,
            label_columns=label_columns,
            tmb_cmap=self.tmb_cmap,
            snv_cmap=self.snv_cmap,
            cnv_cmap=self.cnv_cmap,
            mutsig_cmap=self.mutsig_cmap,
            meta_cmaps=self.meta_cmaps_condensed,
        )

    def make_comut(self):
        self.layout.add_panels()

        def plot_comut(ax):
            self.plotter.plot_cnv_heatmap(
                ax=ax,
                cnv=self.data.cnv.df,
                cnv_cmap=self.cnv_cmap,
                inter_heatmap_linewidth=self.layout.inter_heatmap_linewidth,
                aspect_ratio=self.layout.aspect_ratio
            )
            self.plotter.plot_snv_heatmap(
                ax=ax,
                snv=self.data.snv.df,
                snv_cmap=self.snv_cmap
            )
            self.plotter.plot_heatmap_layout(
                ax=ax,
                cna=self.data.cnv.df,
                labelbottom=self.layout.show_patient_names and self.data.meta is None
            )

        def meta_data_color(column, value):
            cmap = self.meta_cmaps[column]
            if isinstance(cmap, Palette):
                if value not in cmap:
                    raise Exception(f"Something went wrong for column {column}: {value} not in {cmap}!")
                return cmap[value]
            else:
                _cmap, _norm = cmap
                return _cmap(_norm(value))

        self.layout.set_plot_func(
            "comutation",
            plot_comut
        )

        # self.layout.set_plot_func("coverage", self.plotter.plot_coverage)
        self.layout.set_plot_func(
            "mutational signatures",
            self.plotter.plot_mutsig,
            mutsig=self.data.mutsig,
            mutsig_cmap=self.mutsig_cmap
        )
        self.layout.set_plot_func(
            "tmb",
            self.plotter.plot_tmb,
            tmb=self.data.tmb,
            ytickpad=0,
            fontsize=6,
            aspect_ratio=self.layout.aspect_ratio
        )
        # self.layout.set_plot_func("tmb legend", self.plotter.plot_legend, cmap=self.tmb_cmap)

        self.layout.set_plot_func(
            "model annotation",
            self.plotter.plot_model_annotation,
            model_annotation=self.data.get_model_annotation()
        )
        self.layout.set_plot_func(
            "gene names",
            self.plotter.plot_gene_names,
            genes=self.data.genes,
            ground_truth_genes=self.data.ground_truth_genes
        )
        self.layout.set_plot_func(
            "cytoband",
            self.plotter.plot_cytoband,
            cytobands=self.data.cnv.gistic.cytoband
        )
        self.layout.set_plot_func(
            "total recurrence",
            self.plotter.plot_total_recurrence,
            total_recurrence_per_gene=self.data.get_total_recurrence(),
            pad=0.01
        )
        self.layout.set_plot_func(
            "total recurrence overall",
            self.plotter.plot_total_recurrence_overall,
            total_recurrence_overall=self.data.get_total_recurrence_overall(),
            pad=0.01
        )
        self.layout.set_plot_func(
            "recurrence",
            self.plotter.plot_recurrence,
            mut_protein_data=self.data.snv.get_recurrence(index=self.data.genes),
            # mut_prevalence_counter=self.data.snv.get_effect_prevalence().reindex(index=self.data.genes),
            cna=self.data.cnv.df,
            cna_counter=self.data.cnv.get_prevalence(),
            # num_effects=self.data.snv.num_effects,
            genes=self.data.genes,
            columns=self.data.columns,
            cnv_cmap=self.cnv_cmap,
            amp_thresholds=[self.data.high_amp_threshold, self.data.mid_amp_threshold, self.data.low_amp_threshold],
            del_thresholds=[self.data.high_del_threshold, self.data.mid_del_threshold, self.data.low_del_threshold],
            snv_recurrence_threshold=self.data.snv_recurrence_threshold
        )

        self.layout.set_plot_func(
            "model significance",
            self.plotter.plot_model_significance,
            model_significance=self.data.model_significance
        )
        self.layout.set_plot_func(
            "mutsig legend",
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
            names=self.data.model_names,
            title="Significant by Model"
        )

        self.layout.set_plot_func(
            "meta data",
            self.plotter.plot_meta_data,
            meta_data=self.data.meta.df,
            meta_data_color=meta_data_color,
            legend_titles=self.data.meta.legend_titles,
            inter_heatmap_linewidth=self.layout.inter_heatmap_linewidth,
            aspect_ratio=self.layout.aspect_ratio,
            labelbottom=self.layout.show_patient_names
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

        return self.data.genes, self.data.columns
