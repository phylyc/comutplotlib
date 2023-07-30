import pandas as pd

from src.layout import Layout


class ComutLayout(Layout):

    def __init__(
        self,
        panels_to_plot: list[str],
        n_genes: int, n_samples: int, n_meta: int = 0, pad: int = 1,
        xfigsize: float = None, max_xfigsize: float = None, yfigsize: float = None,
        label_columns=False,
        tmb_cmap=(), snv_cmap=(), cnv_cmap=(), meta_cmaps={}
    ):
        self.panels_to_plot = panels_to_plot

        self.meta_data_legend_titles = list(meta_cmaps.keys())

        self.small_inter_legend_width = 2
        self.inter_legend_width = 7
        self.inter_legend_height = 3
        self.inter_heatmap_linewidth = 0.05

        # HEIGHTS
        tmb_height = 7
        mut_sig_height = 7
        comut_height = n_genes
        tmb_legend_height = len(tmb_cmap)
        snv_legend_height = len(snv_cmap)
        cnv_legend_height = len(cnv_cmap)
        model_annotation_legend_height = 2
        meta_height = n_meta
        meta_legend_continuous_height = 4
        meta_legend_heights = {
            legend_title: len(cmap) if isinstance(cmap, dict) else meta_legend_continuous_height
            for legend_title, cmap in meta_cmaps.items()
        }
        max_meta_legend_height = max(meta_legend_heights.values()) if len(meta_legend_heights) else 1

        # WIDTHS
        snv_recurrence_width = 4
        snv_prevalence_width = 3
        cnv_prevalence_width = snv_prevalence_width
        total_prevalence_width = 2
        cytoband_width = 2
        gene_label_width = 6
        model_annotation_width = 2
        model_significance_width = 3
        legend_width = 2
        num_legends = len(meta_cmaps)
        meta_legend_width = num_legends * (legend_width + self.inter_legend_width) - self.inter_legend_width

        left_of_comut_width = 0
        if "recurrence" in panels_to_plot:
            left_of_comut_width += snv_recurrence_width + pad
        if "prevalence" in panels_to_plot:
            left_of_comut_width += snv_prevalence_width + cnv_prevalence_width + pad
        if "total prevalence" in panels_to_plot:
            left_of_comut_width += total_prevalence_width + pad
        if "cytoband" in panels_to_plot:
            left_of_comut_width += cytoband_width
        if "gene names" in panels_to_plot:
            left_of_comut_width += gene_label_width
        if "model annotation" in panels_to_plot:
            left_of_comut_width += model_annotation_width

        non_heatmap_width = left_of_comut_width
        if any([part in panels_to_plot for part in ["tmb legend", "snv legend", "cnv legend", "model annotation legend"]]):
            non_heatmap_width += self.small_inter_legend_width + legend_width
        if "model significance" in panels_to_plot:
            non_heatmap_width += model_significance_width
        comut_width = (
            min(n_samples, int(10 * max_xfigsize - non_heatmap_width))
            if max_xfigsize is not None
            else n_samples
        )
        right_of_comut_width = 0
        if any([part in panels_to_plot for part in ["tmb legend", "snv legend", "cnv legend", "model annotation legend"]]):
            right_of_comut_width += self.small_inter_legend_width + legend_width
        if "model significance" in panels_to_plot:
            right_of_comut_width += model_significance_width
        if "meta data legend" in panels_to_plot:
            right_of_comut_width = max(right_of_comut_width, meta_legend_width - comut_width)

        self.aspect_ratio = comut_width / n_samples
        self.show_patient_names = (comut_width >= n_samples) and label_columns
        self.column_names_height = 5

        xsize = left_of_comut_width + comut_width + right_of_comut_width
        ysize = tmb_height + pad + comut_height
        if self.show_patient_names:
            ysize += self.column_names_height
        if "meta data" in panels_to_plot:
            ysize += pad + meta_height
        if "meta data legend" in panels_to_plot:
            ysize += 3 * pad + max_meta_legend_height

        xfigsize = xfigsize if xfigsize is not None else xsize / 10
        yfigsize = yfigsize if yfigsize is not None else ysize / 10

        self.dimensions = pd.DataFrame.from_dict(
            {
                "comutation": [comut_width, comut_height],

                "mutational signatures": [comut_width, mut_sig_height],
                "tmb": [comut_width, tmb_height],
                "tmb legend": [legend_width, tmb_legend_height],

                "model annotation": [model_annotation_width, comut_height],
                "gene names": [gene_label_width, comut_height],
                "cytoband": [cytoband_width, comut_height],
                "total prevalence": [total_prevalence_width, comut_height],
                "total prevalence overall": [total_prevalence_width, 1],
                "prevalence": [snv_prevalence_width + cnv_prevalence_width, comut_height],
                "recurrence": [snv_recurrence_width, comut_height],

                "model significance": [model_significance_width, comut_height],
                "snv legend": [legend_width, snv_legend_height],
                "cnv legend": [legend_width, cnv_legend_height],
                "model annotation legend": [legend_width, model_annotation_legend_height],

                "meta data": [comut_width, meta_height],
            } | {
                f"meta data legend {title}": [legend_width, meta_legend_heights[title]]
                for title in self.meta_data_legend_titles
            },
            orient="index",
            columns=["width", "height"],
        )

        super().__init__(xfigsize=xfigsize, yfigsize=yfigsize, pad=pad)

    def add_panel(self, name, **kwargs):
        return super().add_panel(name=name, **self.dimensions.loc[name].to_dict(), **kwargs)

    def add_panels(self):
        # CENTER PANEL - CORE
        assert "comutation" in self.panels_to_plot
        p_comut = self.add_panel(name="comutation")

        # LEFT PANELS
        p_ref = p_comut
        if "model annotation" in self.panels_to_plot:
            p_ref = self.add_panel(name="model annotation", left_of=p_ref, pad=0)
        if "gene names" in self.panels_to_plot:
            p_ref = self.add_panel(name="gene names", left_of=p_ref, pad=0)
        if "cytoband" in self.panels_to_plot:
            p_ref = self.add_panel(name="cytoband", left_of=p_ref, pad=0)
        if "total prevalence" in self.panels_to_plot:
            p_ref = self.add_panel(name="total prevalence", left_of=p_ref)
            if "total prevalence overall" in self.panels_to_plot:
                self.add_panel(name="total prevalence overall", below=p_ref)
        if "prevalence" in self.panels_to_plot:
            p_ref = self.add_panel(name="prevalence", left_of=p_ref)
        if "recurrence" in self.panels_to_plot:
            p_ref = self.add_panel(name="recurrence", left_of=p_ref)

        # TOP PANELS
        p_ref = p_comut
        if "mutational signatures" in self.panels_to_plot:
            p_ref = self.add_panel(name="mutational signatures", above=p_ref)
        if "tmb" in self.panels_to_plot:
            p_ref = self.add_panel(name="tmb", above=p_ref)
            if "tmb legend" in self.panels_to_plot:
                self.add_panel(name="tmb legend", right_of=p_ref, pad=self.small_inter_legend_width, align="top")

        # RIGHT PANELS
        p_ref = p_comut
        if "model significance" in self.panels_to_plot:
            p_ref = self.add_panel(name="model significance", right_of=p_ref)

        def first_legend_panel(name):
            return self.add_panel(name=name, right_of=p_ref, pad=self.small_inter_legend_width, align="top")

        def other_legend_panel(name, p_ref):
            return self.add_panel(name=name, below=p_ref, pad=self.inter_legend_height, align="left")

        if "snv legend" in self.panels_to_plot:
            p_ref = first_legend_panel(name="snv legend")
        if "cnv legend" in self.panels_to_plot:
            p_ref = (
                first_legend_panel(name="cnv legend")
                if p_ref.name == "comutation"
                else other_legend_panel(name="cnv legend", p_ref=p_ref)
            )
        if "model annotation legend" in self.panels_to_plot:
            p_ref = (
                first_legend_panel(name="model annotation legend")
                if p_ref.name == "comutation"
                else other_legend_panel(name="model annotation legend", p_ref=p_ref)
            )

        # BOTTOM PANELS
        p_ref = p_comut
        if "meta data" in self.panels_to_plot:
            p_ref = self.add_panel(name="meta data", below=p_ref)
            if "meta data legend" in self.panels_to_plot and len(self.meta_data_legend_titles):
                pad = self.inter_legend_height + (self.column_names_height if self.show_patient_names else 0)
                p_ref = self.add_panel(
                    name=f"meta data legend {self.meta_data_legend_titles[0]}", below=p_ref, pad=pad, align="left"
                )
                for title in self.meta_data_legend_titles[1:]:
                    p_ref = self.add_panel(name=f"meta data legend {title}", right_of=p_ref, pad=self.inter_legend_width, align="top")

        self.place_panels_on_gridspec(autoscale_figsize=True, scale=0.1)

    def set_plot_func(self, panel: str, plot_func, *args, **kwargs):
        if panel in self.panels:
            self.panels[panel].set_plot_func(plot_func, *args, **kwargs)
