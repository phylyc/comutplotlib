import pandas as pd

from comutplotlib.layout import Layout
from comutplotlib.palette import Palette
from comutplotlib.comut_panels import ComutPanels, control, meta_legend


class ComutLayout(Layout):

    def __init__(
        self,
        panels_to_plot: list[str],
        n_genes: int, n_samples: int, n_samples_control: int = 0, n_meta: int = 0, n_meta_genes: int = 0, pad: int = 1,
        xfigsize: float | None = None, max_xfigsize: float | None = None, max_xfigsize_scale: float = 1, yfigsize: float | None = None,
        label_columns=False,
        tmb_cmap=(), snv_cmap=(), cnv_cmap=(), mutsig_cmap=(), meta_cmaps=None,
    ):
        meta_cmaps = meta_cmaps if meta_cmaps is not None else {}
        self.panels_to_plot = panels_to_plot

        self.meta_data_legend_titles = list(meta_cmaps.keys())

        self.small_inter_legend_width = 2
        self.inter_legend_width = 5
        self.inter_legend_height = 3
        self.inter_heatmap_linewidth = 0.05

        # HEIGHTS
        title_height = 2
        tmb_height = 4
        mut_sig_height = 4
        coverage_height = 5
        comut_height = n_genes
        tmb_legend_height = len(tmb_cmap)
        snv_legend_height = len(snv_cmap)
        cnv_legend_height = len(cnv_cmap)
        mutsig_legend_height = len(mutsig_cmap)
        model_annotation_legend_height = 2
        meta_height = n_meta
        meta_legend_continuous_height = 4
        meta_legend_heights = {
            legend_title: len(cmap.drop_undefined()) if isinstance(cmap, Palette) else meta_legend_continuous_height
            for legend_title, cmap in meta_cmaps.items()
        }
        max_meta_legend_height = max(meta_legend_heights.values()) if len(meta_legend_heights) else 1

        # WIDTHS
        recurrence_width = 5
        recurrence_fold_change_width = 10
        total_recurrence_width = 2
        cytoband_width = 2
        gene_label_width = 6
        model_annotation_width = 2
        model_significance_width = 3
        legend_width = 2
        genes_meta_width = n_meta_genes
        num_legends = len(meta_cmaps)
        meta_legend_width = num_legends * (legend_width + self.inter_legend_width) - self.inter_legend_width

        left_of_comut_width = 0
        if ComutPanels.gene_meta_data in panels_to_plot:
            left_of_comut_width += genes_meta_width + (pad if ComutPanels.cytoband in panels_to_plot else 0)
        if ComutPanels.cytoband in panels_to_plot:
            left_of_comut_width += cytoband_width
        if ComutPanels.gene_names in panels_to_plot:
            left_of_comut_width += gene_label_width
        if ComutPanels.model_annotation in panels_to_plot:
            left_of_comut_width += model_annotation_width

        right_panel_legends = [ComutPanels.tmb_legend, ComutPanels.mutational_signatures_legend, ComutPanels.snv_legend, ComutPanels.cnv_legend, ComutPanels.model_annotation_legend]

        non_heatmap_width = left_of_comut_width
        if any([part in panels_to_plot for part in right_panel_legends]):
            non_heatmap_width += self.small_inter_legend_width + legend_width
        if ComutPanels.model_significance in panels_to_plot:
            non_heatmap_width += model_significance_width
        if ComutPanels.recurrence in panels_to_plot:
            non_heatmap_width += recurrence_width + pad
        if control(ComutPanels.recurrence) in panels_to_plot:
            non_heatmap_width += recurrence_width + (pad if ComutPanels.recurrence_fold_change in panels_to_plot else 0)
        if ComutPanels.recurrence_fold_change in panels_to_plot:
            non_heatmap_width += recurrence_fold_change_width + pad
        if ComutPanels.total_recurrence in panels_to_plot:
            non_heatmap_width += total_recurrence_width + pad
        comut_width = (
            max(min(n_samples, int(10 * max_xfigsize - non_heatmap_width)), 1)
            if max_xfigsize is not None
            else int(max_xfigsize_scale * n_samples)
        ) if ComutPanels.comutation in panels_to_plot else 0
        comut_width_control = (
            max(min(n_samples_control, int(10 * max_xfigsize - non_heatmap_width)), 1)
            if max_xfigsize is not None
            else int(max_xfigsize_scale * n_samples_control)
        ) if control(ComutPanels.comutation) in panels_to_plot else 0
        right_of_comut_width = 0
        if any([part in panels_to_plot for part in right_panel_legends]):
            right_of_comut_width += self.small_inter_legend_width + legend_width
        if ComutPanels.model_significance in panels_to_plot:
            right_of_comut_width += model_significance_width
        if ComutPanels.meta_data_legend in panels_to_plot:
            right_of_comut_width = max(right_of_comut_width, meta_legend_width - comut_width)

        self.aspect_ratio = comut_width / n_samples
        self.show_patient_names = (comut_width >= n_samples) and label_columns
        self.column_names_height = 5

        xsize = left_of_comut_width + comut_width + right_of_comut_width
        ysize = title_height + tmb_height + pad + comut_height
        if ComutPanels.total_recurrence_fold_change in panels_to_plot:
            ysize += 1
        if self.show_patient_names:
            ysize += self.column_names_height
        if ComutPanels.meta_data in panels_to_plot:
            ysize += pad + meta_height
        if ComutPanels.meta_data_legend in panels_to_plot:
            ysize += 3 * pad + max_meta_legend_height

        xfigsize = xfigsize if xfigsize is not None else xsize / 10
        yfigsize = yfigsize if yfigsize is not None else ysize / 10

        self.dimensions = pd.DataFrame.from_dict(
            {
                ComutPanels.anchor: [0, n_genes],

                ComutPanels.comutation: [comut_width, comut_height],
                control(ComutPanels.comutation): [comut_width_control, comut_height],
                ComutPanels.snv_legend: [legend_width, snv_legend_height],
                ComutPanels.cnv_legend: [legend_width, cnv_legend_height],

                ComutPanels.model_significance: [model_significance_width, comut_height],

                ComutPanels.model_annotation: [model_annotation_width, comut_height],
                ComutPanels.model_annotation_legend: [legend_width, model_annotation_legend_height],

                ComutPanels.cohort_label: [comut_width, title_height],
                control(ComutPanels.cohort_label): [comut_width, title_height],
                ComutPanels.coverage: [comut_width, coverage_height],
                control(ComutPanels.coverage): [comut_width_control, coverage_height],
                ComutPanels.tmb: [comut_width, tmb_height],
                control(ComutPanels.tmb): [comut_width_control, tmb_height],
                ComutPanels.tmb_legend: [legend_width, tmb_legend_height],

                ComutPanels.mutational_signatures: [comut_width, mut_sig_height],
                control(ComutPanels.mutational_signatures): [comut_width_control, mut_sig_height],
                ComutPanels.mutational_signatures_legend: [legend_width, mutsig_legend_height],

                ComutPanels.gene_names: [gene_label_width, comut_height],
                ComutPanels.cytoband: [cytoband_width, comut_height],
                ComutPanels.gene_meta_data: [genes_meta_width, comut_height],

                ComutPanels.total_recurrence: [total_recurrence_width, comut_height],
                control(ComutPanels.total_recurrence): [total_recurrence_width, comut_height],
                ComutPanels.recurrence: [recurrence_width, comut_height],
                control(ComutPanels.recurrence): [recurrence_width, comut_height],
                ComutPanels.total_recurrence_overall: [recurrence_width, 1],
                control(ComutPanels.total_recurrence_overall): [recurrence_width, 1],
                ComutPanels.recurrence_fold_change: [recurrence_fold_change_width, comut_height],
                ComutPanels.total_recurrence_fold_change: [recurrence_fold_change_width, 1],

                ComutPanels.meta_data: [comut_width, meta_height],
                control(ComutPanels.meta_data): [comut_width_control, meta_height],
            } | {
                meta_legend(title): [legend_width, meta_legend_heights[title]]
                for title in self.meta_data_legend_titles
            },
            orient="index",
            columns=["width", "height"],
        )

        super().__init__(xfigsize=xfigsize, yfigsize=yfigsize, pad=pad)

    def add_panel(self, name, ref=None, force_add=False, **kwargs):
        if name in self.panels_to_plot or force_add:
            return super().add_panel(name=name, **self.dimensions.loc[name].to_dict(), **kwargs)
        else:
            return ref

    def add_panels(self):
        p_anchor = self.add_panel(name=ComutPanels.anchor, force_add=True)

        # CENTER PANEL - CORE
        p_comut = self.add_panel(name=ComutPanels.comutation, ref=p_anchor, right_of=p_anchor, pad=0)

        # LEFT PANELS
        p_ref = p_anchor
        # for panel, pad in zip(["model significance", "model annotation", "gene names", "cytoband", "recurrence"], [1, 0, 0, 0, 1]):
        #     p_ref = self.add_panel(name=panel, ref=p_ref, left_of=p_ref, pad=pad)
        p_ref = self.add_panel(name=ComutPanels.model_annotation, ref=p_ref, left_of=p_ref, pad=0)
        p_ref = self.add_panel(name=ComutPanels.gene_names, ref=p_ref, left_of=p_ref, pad=0)
        p_ref = self.add_panel(name=ComutPanels.cytoband, ref=p_ref, left_of=p_ref, pad=0 if ComutPanels.gene_names in self.panels_to_plot else self.pad)
        p_ref = self.add_panel(name=ComutPanels.gene_meta_data, ref=p_ref, left_of=p_ref, pad=self.pad if ComutPanels.cytoband in self.panels_to_plot else 0)

        # TOP PANELS
        p_ref = p_anchor
        for panel in [ComutPanels.mutational_signatures, ComutPanels.coverage, ComutPanels.tmb, ComutPanels.cohort_label]:
            p_ref = self.add_panel(name=panel, ref=p_ref, above=p_ref, align="left")

        # BOTTOM PANELS
        p_ref = p_anchor
        p_ref = self.add_panel(name=ComutPanels.meta_data, ref=p_ref, below=p_ref, align="left")
        if ComutPanels.meta_data_legend in self.panels_to_plot and len(self.meta_data_legend_titles):
            pad = self.pad + (self.column_names_height if self.show_patient_names else 0)
            p_ref = self.add_panel(
                name=meta_legend(self.meta_data_legend_titles[0]), ref=p_ref, force_add=True, below=p_ref, pad=pad, align="left"
            )
            for title in self.meta_data_legend_titles[1:]:
                p_ref = self.add_panel(name=meta_legend(title), ref=p_ref, force_add=True, right_of=p_ref, pad=self.inter_legend_width, align="top")

        # RIGHT PANELS
        p_ref = p_comut
        p_ref = self.add_panel(name=ComutPanels.recurrence, ref=p_ref, right_of=p_ref)
        # p_ref = self.add_panel(name="total recurrence", ref=p_ref, right_of=p_ref)
        self.add_panel(name=ComutPanels.total_recurrence_overall, ref=p_ref, below=p_ref, pad=0)

        # LEFT PANELS
        p_ref = self.add_panel(name=ComutPanels.recurrence_fold_change, ref=p_ref, right_of=p_ref)
        self.add_panel(name=ComutPanels.total_recurrence_fold_change, ref=p_ref, below=p_ref, pad=0)
        p_ref = self.add_panel(name=control(ComutPanels.recurrence), ref=p_ref, right_of=p_ref, pad=self.pad if ComutPanels.recurrence_fold_change in self.panels_to_plot else 0)
        # p_ref = self.add_panel(name="total recurrence control", ref=p_ref, right_of=p_ref, pad=0)
        self.add_panel(name=control(ComutPanels.total_recurrence_overall), ref=p_ref, below=p_ref, pad=0)
        p_rec_control = p_ref

        p_ref = self.add_panel(name=control(ComutPanels.comutation), ref=p_ref, right_of=p_ref)
        p_comut_control = p_ref

        # TOP PANELS
        p_ref = p_comut_control
        for panel in [control(ComutPanels.mutational_signatures), control(ComutPanels.coverage), control(ComutPanels.tmb), control(ComutPanels.cohort_label)]:
            p_ref = self.add_panel(name=panel, ref=p_ref, above=p_ref)

        # BOTTOM PANELS
        p_ref = p_comut_control
        p_ref = self.add_panel(name=control(ComutPanels.meta_data), ref=p_ref, below=p_ref)

        # RIGHT PANELS
        # p_ref = p_comut_control
        # p_ref = self.add_panel(name="recurrence control", ref=p_ref, right_of=p_ref)

        p_ref = p_comut_control
        # p_ref = p_rec_control

        def first_legend_panel(name, p_ref):
            return self.add_panel(name=name, ref=p_ref, right_of=p_ref, pad=self.small_inter_legend_width, align="top")

        def above_legend_panel(name, p_ref):
            return self.add_panel(name=name, ref=p_ref, above=p_ref, pad=self.inter_legend_height, align="left")

        def below_legend_panel(name, p_ref):
            return self.add_panel(name=name, ref=p_ref, below=p_ref, pad=self.inter_legend_height, align="left")

        has_legend = False
        p_ref_top = None
        for panel in [ComutPanels.tmb_legend, ComutPanels.mutational_signatures_legend]:
            if panel in self.panels_to_plot:
                p_ref = (
                    first_legend_panel(name=panel, p_ref=p_ref)
                    if not has_legend
                    else below_legend_panel(name=panel, p_ref=p_ref_top)
                )
                p_ref_top = p_ref if p_ref_top is None else p_ref_top
                has_legend = True
        for panel in [ComutPanels.cnv_legend, ComutPanels.snv_legend, ComutPanels.model_annotation_legend]:
            if panel in self.panels_to_plot:
                p_ref = (
                    first_legend_panel(name=panel, p_ref=p_ref)
                    if not has_legend
                    else below_legend_panel(name=panel, p_ref=p_ref)
                )
                p_ref_top = p_ref if p_ref_top is None else p_ref_top
                has_legend = True
        self.place_panels_on_gridspec(autoscale_figsize=True, scale=0.1)

    def set_plot_func(self, panel: str, plot_func, *args, **kwargs):
        if panel in self.panels:
            self.panels[panel].set_plot_func(plot_func, *args, **kwargs)
