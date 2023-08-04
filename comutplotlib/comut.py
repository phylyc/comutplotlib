from src.comut_data import ComutData
from src.comut_layout import ComutLayout
from src.comut_plotter import ComutPlotter
from src.mutation_annotation import MutationAnnotation as MutA


class Comut(object):

    def __init__(
        self,
        output: str = "./comut.pdf",

        maf: list[str] = (),
        maf_pool_as: dict | None = None,

        seg: list[str] = (),
        gistic: list[str] = (),

        mutsig: list[str] = (),

        model_significances: list[str] = (),
        model_names: list[str] = (),

        sif: list[str] = (),
        meta_data_rows: list[str] = (),
        meta_data_rows_per_sample: list[str] = (),

        by: str = MutA.patient,

        column_order: tuple[str] = None,
        index_order: tuple[str] = None,

        interesting_gene: str = None,
        interesting_gene_comut_percent_threshold: float = None,
        interesting_genes: set = None,
        snv_interesting_genes: set = None,
        cnv_interesting_genes: set = None,
        total_prevalence_threshold: float = None,

        low_amp_threshold: int = 1,
        high_amp_threshold: int = 2,
        low_del_threshold: int = -1,
        high_del_threshold: int = -2,

        panels_to_plot: list = None,
        palette: dict[str, tuple[float]] = None,
        ground_truth_genes: dict[str, list[str]] = None,  # todo: refactor as a palette class
        max_xfigsize: int = None,
        label_columns: bool = False,
        **kwargs
    ):
        self.plotter = ComutPlotter(
            output=output,
            extra_palette=palette
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
            by=by,
            column_order=column_order,
            index_order=index_order,
            interesting_gene=interesting_gene,
            interesting_gene_comut_percent_threshold=interesting_gene_comut_percent_threshold,
            interesting_genes=interesting_genes,
            ground_truth_genes=ground_truth_genes,
            snv_interesting_genes=snv_interesting_genes,
            cnv_interesting_genes=cnv_interesting_genes,
            total_prevalence_threshold=total_prevalence_threshold,
            low_amp_threshold=low_amp_threshold,
            high_amp_threshold=high_amp_threshold,
            low_del_threshold=low_del_threshold,
            high_del_threshold=high_del_threshold,
        )
        self.data.preprocess()

        n_genes, n_samples, n_meta = self.data.get_dimensions()

        self.tmb_cmap = self.plotter.palette.get_tmb_cmap(self.data)
        self.snv_cmap = self.plotter.palette.get_snv_cmap(self.data)
        self.cnv_cmap, self.cnv_names = self.plotter.palette.get_cnv_cmap(self.data)
        self.mutsig_cmap = self.plotter.palette.get_mutsig_cmap(self.data)
        self.meta_cmaps = self.plotter.palette.get_meta_cmaps(self.data)

        self.layout = ComutLayout(
            panels_to_plot=panels_to_plot,
            max_xfigsize=max_xfigsize,
            n_genes=n_genes,
            n_samples=n_samples,
            n_meta=n_meta,
            label_columns=label_columns,
            tmb_cmap=self.tmb_cmap,
            snv_cmap=self.snv_cmap,
            cnv_cmap=self.cnv_cmap,
            mutsig_cmap=self.mutsig_cmap,
            meta_cmaps=self.meta_cmaps,
        )

    def make_comut(self):
        self.layout.add_panels()

        def plot_comut(ax):
            self.plotter.plot_cnv_heatmap(ax=ax, cnv=self.data.cnv.df, cnv_cmap=self.cnv_cmap, inter_heatmap_linewidth=self.layout.inter_heatmap_linewidth, aspect_ratio=self.layout.aspect_ratio)
            self.plotter.plot_snv_heatmap(ax=ax, snv=self.data.snv.df, snv_cmap=self.snv_cmap)
            self.plotter.plot_heatmap_layout(ax=ax, cna=self.data.cnv.df, labelbottom=self.layout.show_patient_names and self.data.meta is None)

        def meta_data_color(column, value):
            cmap = self.meta_cmaps[column]
            return cmap[value] if isinstance(cmap, dict) else cmap(value)

        self.layout.set_plot_func("comutation", plot_comut)

        self.layout.set_plot_func("mutational signatures", self.plotter.plot_mutsig, mutsig=self.data.mutsig, mutsig_cmap=self.mutsig_cmap)
        self.layout.set_plot_func("tmb", self.plotter.plot_tmb, tmb=self.data.tmb, ytickpad=0, fontsize=6)
        # self.layout.set_plot_func("tmb legend", self.plotter.plot_legend, cmap=self.tmb_cmap)

        self.layout.set_plot_func("model annotation", self.plotter.plot_model_annotation, model_annotation=self.data.get_model_annotation())
        self.layout.set_plot_func("gene names", self.plotter.plot_gene_names, genes=self.data.genes, ground_truth_genes=self.data.ground_truth_genes)
        self.layout.set_plot_func("cytoband", self.plotter.plot_cytoband, cytobands=self.data.cnv.gistic.cytoband)
        self.layout.set_plot_func("total prevalence", self.plotter.plot_total_prevalence, total_prevalence_per_gene=self.data.get_total_prevalence(), pad=0.01)
        self.layout.set_plot_func("total prevalence overall", self.plotter.plot_total_prevalence_overall, total_prevalence_overall=self.data.get_total_prevalence_overall(), pad=0.01)
        self.layout.set_plot_func("prevalence", self.plotter.plot_prevalence, self.data.snv.get_effect_prevalence().reindex(index=self.data.genes), self.data.cnv.df, self.data.cnv.get_prevalence(), self.data.snv.num_effects, self.data.genes, self.data.columns, self.cnv_cmap, [self.data.high_amp_threshold, self.data.low_amp_threshold], [self.data.high_del_threshold, self.data.low_del_threshold])
        self.layout.set_plot_func("recurrence", self.plotter.plot_recurrence, mut_protein_data=self.data.snv.get_recurrence(index=self.data.genes), genes=self.data.genes, snv_recurrence_threshold=self.data.snv_recurrence_threshold)

        self.layout.set_plot_func("model significance", self.plotter.plot_model_significance, model_significance=self.data.model_significance)
        self.layout.set_plot_func("mutsig legend", self.plotter.plot_legend, cmap=self.mutsig_cmap, title="Mutational Signatures")
        self.layout.set_plot_func("snv legend", self.plotter.plot_legend, cmap=self.snv_cmap, title="Short Nucleotide Alteration")
        self.layout.set_plot_func("cnv legend", self.plotter.plot_legend, cmap=self.cnv_cmap, names=self.cnv_names, title="Copy Number Alteration")
        self.layout.set_plot_func("model annotation legend", self.plotter.plot_model_annotation_legend, names=self.data.model_names, title="Significant by Model")

        self.layout.set_plot_func("meta data", self.plotter.plot_meta_data, meta_data=self.data.meta.df, meta_data_color=meta_data_color, legend_titles=self.data.meta.legend_titles, inter_heatmap_linewidth=self.layout.inter_heatmap_linewidth, aspect_ratio=self.layout.aspect_ratio, labelbottom=self.layout.show_patient_names)
        for title in self.data.meta.legend_titles:
            self.layout.set_plot_func(f"meta data legend {title}", self.plotter.plot_legend, cmap=self.meta_cmaps[title], title=title)

        for _, panel in self.layout.panels.items():
            if panel.ax is not None and panel.plot_func is not None:  # some panels may have zero size and are not placed on the gridspec
                panel.plot_func(panel.ax)

        self.plotter.save_figure(fig=self.layout.fig, bbox_inches="tight")
