import pandas as pd

from comutplotlib.layout import Layout
from comutplotlib.palette import Palette
from comutplotlib.comut_panels import (
    ComutPanels,
    control,
    meta_legend,
    block,
    grid_row_panel,
    grid_col_panel,
    parse_grid_panel,
)


class ComutLayout(Layout):

    def __init__(
        self,
        panels_to_plot: list[str],
        n_genes: int, n_samples: int, n_samples_control: int = 0, n_meta: int = 0, n_meta_genes: int = 0, pad: int = 1,
        xfigsize: float | None = None, max_xfigsize: float | None = None, max_xfigsize_scale: float = 1, yfigsize: float | None = None,
        label_columns=False,
        tmb_cmap=(), snv_cmap=(), cnv_cmap=(), mutsig_cmap=(), meta_cmaps=None,
        grid=None, control_column_group_sizes=None, control_column_group_labels=None,
        meta_data_position: str = "bottom",
        control_position: str = "right",
    ):
        meta_cmaps = meta_cmaps if meta_cmaps is not None else {}
        if meta_data_position not in {"top", "bottom"}:
            raise ValueError("meta_data_position must be either 'top' or 'bottom'.")
        if control_position not in {"left", "right"}:
            raise ValueError("control_position must be either 'left' or 'right'.")
        self.panels_to_plot = panels_to_plot
        self.meta_data_position = meta_data_position
        # Only mirror the cohorts when a control cohort is actually present.
        self.control_on_left = (control_position == "left") and (n_samples_control > 0)
        self.grid = grid
        self._control_column_group_sizes = list(control_column_group_sizes) if control_column_group_sizes is not None else []
        self._control_column_group_labels = list(control_column_group_labels) if control_column_group_labels is not None else []

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
        recurrence_width = 4
        recurrence_fold_change_width = 5
        total_recurrence_width = 2
        cytoband_width = 2
        gene_label_width = 6
        model_annotation_width = 2
        gene_group_label_width = 2
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
                ComutPanels.column_group_label: [comut_width, title_height],
                control(ComutPanels.column_group_label): [comut_width_control, title_height],
                ComutPanels.gene_group_label: [gene_group_label_width, comut_height],
            } | {
                meta_legend(title): [legend_width, meta_legend_heights[title]]
                for title in self.meta_data_legend_titles
            },
            orient="index",
            columns=["width", "height"],
        )

        super().__init__(xfigsize=xfigsize, yfigsize=yfigsize, pad=pad)

        # Grid sub-panel sizing (see GRID_SUBPANELS_PLAN.md). Per-group widths are
        # distributed proportionally out of the total comut width so the overall
        # figure stays within the same width budget; per-group heights are the
        # gene-group sizes.
        self.grid_enabled = self.grid is not None and not self.grid.is_trivial()
        if self.grid_enabled:
            self.gene_group_heights = [len(g) for g in self.grid.gene_groups]
            col_sizes = [len(g) for g in self.grid.column_groups]
            self.col_group_widths = self._split_width(comut_width, col_sizes, n_samples)
            self.col_group_widths_control = (
                self._split_width(comut_width_control, self._control_column_group_sizes, n_samples_control)
                if n_samples_control and self._control_column_group_sizes else []
            )

    @staticmethod
    def _split_width(total: int, group_sizes: list[int], n_total: int) -> list[int]:
        """Distribute ``total`` grid columns across groups proportional to size."""
        if not group_sizes or not n_total or total <= 0:
            return [max(total, 1) for _ in group_sizes] if group_sizes else []
        return [max(int(round(total * s / n_total)), 1) for s in group_sizes]

    def _fixed_dim(self, base: str) -> tuple[int, int]:
        """Return the (width, height) a base panel declares in ``self.dimensions``."""
        row = self.dimensions.loc[base]
        return int(row["width"]), int(row["height"])

    def add_panel(self, name, ref=None, force_add=False, **kwargs):
        if name in self.panels_to_plot or force_add:
            return super().add_panel(name=name, **self.dimensions.loc[name].to_dict(), **kwargs)
        else:
            return ref

    def add_panels(self):
        if getattr(self, "grid_enabled", False):
            return self.add_panels_grid()

        p_anchor = self.add_panel(name=ComutPanels.anchor, force_add=True)

        # Cohort-side helpers: ``left_panel``/``right_panel`` map a base panel id to
        # the cohort drawn on that geometric side. With ``control_position="right"``
        # (default) the left cohort is the case and the right cohort is the control;
        # ``control_position="left"`` swaps them. The shared gene-axis marginals
        # (gene names, cytoband, ...) always stay on the far left.
        def left_panel(name):
            return control(name) if self.control_on_left else name

        def right_panel(name):
            return name if self.control_on_left else control(name)

        # CENTER PANEL - CORE (left cohort)
        p_comut = self.add_panel(name=left_panel(ComutPanels.comutation), ref=p_anchor, right_of=p_anchor, pad=0)
        draw_comut = left_panel(ComutPanels.comutation) in self.panels_to_plot

        # LEFT PANELS
        p_ref = p_anchor
        # for panel, pad in zip(["model significance", "model annotation", "gene names", "cytoband", "recurrence"], [1, 0, 0, 0, 1]):
        #     p_ref = self.add_panel(name=panel, ref=p_ref, left_of=p_ref, pad=pad)
        p_ref = self.add_panel(name=ComutPanels.model_annotation, ref=p_ref, left_of=p_ref, pad=0)
        p_ref = self.add_panel(name=ComutPanels.gene_names, ref=p_ref, left_of=p_ref, pad=0)
        p_ref = self.add_panel(name=ComutPanels.cytoband, ref=p_ref, left_of=p_ref, pad=0 if ComutPanels.gene_names in self.panels_to_plot else self.pad)
        p_ref = self.add_panel(name=ComutPanels.gene_meta_data, ref=p_ref, left_of=p_ref, pad=self.pad if ComutPanels.cytoband in self.panels_to_plot else 0)

        # TOP PANELS (left cohort)
        p_ref = p_anchor
        if self.meta_data_position == "top":
            p_ref = self.add_panel(name=left_panel(ComutPanels.meta_data), ref=p_ref, above=p_ref, align="left")
        for panel in [ComutPanels.mutational_signatures, ComutPanels.coverage, ComutPanels.tmb, ComutPanels.cohort_label]:
            p_ref = self.add_panel(name=left_panel(panel), ref=p_ref, above=p_ref, align="left")

        # BOTTOM PANELS (left cohort)
        p_ref = p_anchor
        if self.meta_data_position == "bottom":
            p_ref = self.add_panel(name=left_panel(ComutPanels.meta_data), ref=p_ref, below=p_ref, align="left")
        if ComutPanels.meta_data_legend in self.panels_to_plot and len(self.meta_data_legend_titles):
            pad = self.pad + (self.column_names_height if self.show_patient_names else 0)
            p_ref = self.add_panel(
                name=meta_legend(self.meta_data_legend_titles[0]), ref=p_ref, force_add=True, below=p_ref, pad=pad, align="left"
            )
            for title in self.meta_data_legend_titles[1:]:
                p_ref = self.add_panel(name=meta_legend(title), ref=p_ref, force_add=True, right_of=p_ref, pad=self.inter_legend_width, align="top")

        # RIGHT PANELS: left-cohort recurrence -> fold change -> right-cohort recurrence -> right-cohort comut
        p_ref = p_comut
        p_ref = self.add_panel(name=left_panel(ComutPanels.recurrence), ref=p_ref, right_of=p_ref, pad=self.pad if draw_comut else 0)
        # p_ref = self.add_panel(name="total recurrence", ref=p_ref, right_of=p_ref)
        self.add_panel(name=left_panel(ComutPanels.total_recurrence_overall), ref=p_ref, below=p_ref, pad=0)

        # central fold-change panel
        p_ref = self.add_panel(name=ComutPanels.recurrence_fold_change, ref=p_ref, right_of=p_ref)
        self.add_panel(name=ComutPanels.total_recurrence_fold_change, ref=p_ref, below=p_ref, pad=0)
        p_ref = self.add_panel(name=right_panel(ComutPanels.recurrence), ref=p_ref, right_of=p_ref, pad=self.pad if ComutPanels.recurrence_fold_change in self.panels_to_plot else 0)
        # p_ref = self.add_panel(name="total recurrence control", ref=p_ref, right_of=p_ref, pad=0)
        self.add_panel(name=right_panel(ComutPanels.total_recurrence_overall), ref=p_ref, below=p_ref, pad=0)
        p_rec_control = p_ref

        p_ref = self.add_panel(name=right_panel(ComutPanels.comutation), ref=p_ref, right_of=p_ref)
        p_comut_control = p_ref

        # TOP PANELS (right cohort)
        p_ref = p_comut_control
        if self.meta_data_position == "top":
            p_ref = self.add_panel(name=right_panel(ComutPanels.meta_data), ref=p_ref, above=p_ref)
        for panel in [ComutPanels.mutational_signatures, ComutPanels.coverage, ComutPanels.tmb, ComutPanels.cohort_label]:
            p_ref = self.add_panel(name=right_panel(panel), ref=p_ref, above=p_ref)

        # BOTTOM PANELS (right cohort)
        p_ref = p_comut_control
        if self.meta_data_position == "bottom":
            p_ref = self.add_panel(name=right_panel(ComutPanels.meta_data), ref=p_ref, below=p_ref)

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

    # --- Grid tiling ---------------------------------------------------------

    def _grid_add(self, name, width, height, ref=None, force=False, **kwargs):
        """Place one grid panel with explicit dims, honouring base-panel selection.

        Returns the newly placed panel, or the fallback ``ref`` when the base
        panel is not selected or the panel would be empty (so placement chains
        keep a valid reference).
        """
        base, _, _ = parse_grid_panel(name)
        if not (force or base in self.panels_to_plot):
            return ref
        if width <= 0 or height <= 0:
            return ref
        return Layout.add_panel(self, name=name, width=width, height=height, **kwargs)

    def add_panels_grid(self):
        N = self.grid.n_rows
        comut_height = sum(self.gene_group_heights) + (N - 1) * self.pad
        p_anchor = Layout.add_panel(self, name=ComutPanels.anchor, width=0, height=comut_height, x=0, y=0)

        # Cohort-side descriptors. The LEFT cohort provides the scaffold (defines
        # every row's height/y and carries the shared left gene-axis marginals);
        # the RIGHT cohort is the mirror block set. With ``control_position="right"``
        # (default) the left cohort is the case and the right cohort is the
        # control; ``control_position="left"`` swaps them. Each side keeps its own
        # column groups (widths + labels) and its own base panel names.
        def case_name(name):
            return name

        case_labels = [g.label for g in self.grid.column_groups]
        control_labels = list(self._control_column_group_labels)
        if self.control_on_left:
            left_name, left_widths, left_labels = control, self.col_group_widths_control, control_labels
            right_name, right_widths, right_labels = case_name, self.col_group_widths, case_labels
        else:
            left_name, left_widths, left_labels = case_name, self.col_group_widths, case_labels
            right_name, right_widths, right_labels = control, self.col_group_widths_control, control_labels
        M_left = len(left_widths)

        # --- central comut blocks (left cohort = scaffold) ---
        # The blocks are always placed as the grid *scaffold*: they define every
        # row's height/y and every column's width/x, so the marginal panels align
        # to them. When the comutation heatmap is NOT drawn, the blocks take zero
        # horizontal width (vertical scaffold only) and the inter-group
        # horizontal pads collapse to zero, so that the left/right marginals are
        # not separated by empty columns. The heatmap itself is drawn in
        # ``_render_grid``.
        draw_left_comut = left_name(ComutPanels.comutation) in self.panels_to_plot
        left_col_pad = self.pad if draw_left_comut else 0
        blocks = {}
        for i in range(N):
            hi = self.gene_group_heights[i]
            for j in range(M_left):
                wj = left_widths[j] if draw_left_comut else 0
                name = block(left_name(ComutPanels.comutation), i, j)
                if i == 0 and j == 0:
                    p = Layout.add_panel(self, name=name, width=wj, height=hi, right_of=p_anchor, pad=0, align="top")
                elif j == 0:
                    p = Layout.add_panel(self, name=name, width=wj, height=hi, below=blocks[(i - 1, 0)], pad=self.pad, align="left")
                else:
                    p = Layout.add_panel(self, name=name, width=wj, height=hi, right_of=blocks[(i, j - 1)], pad=left_col_pad, align="top")
                blocks[(i, j)] = p

        # --- top marginals (per left-cohort column group) ---
        top_panels = [ComutPanels.mutational_signatures, ComutPanels.coverage, ComutPanels.tmb, ComutPanels.cohort_label]
        _, meta_h = self._fixed_dim(ComutPanels.meta_data)
        for j in range(M_left):
            p_ref = blocks[(0, j)]
            wj = left_widths[j]
            if self.meta_data_position == "top":
                p_ref = self._grid_add(
                    grid_col_panel(left_name(ComutPanels.meta_data), j), wj, meta_h,
                    ref=p_ref, above=p_ref, align="left",
                )
            for base in top_panels:
                _, h = self._fixed_dim(base)
                p_ref = self._grid_add(grid_col_panel(left_name(base), j), wj, h, ref=p_ref, above=p_ref, align="left")
            # column-group header (drawn only for labelled groups, and only when
            # the ``column_group_label`` panel is requested in panels_to_plot)
            if j < len(left_labels) and left_labels[j] and ComutPanels.column_group_label in self.panels_to_plot:
                _, lh = self._fixed_dim(left_name(ComutPanels.column_group_label))
                self._grid_add(
                    grid_col_panel(left_name(ComutPanels.column_group_label), j), wj, lh,
                    ref=p_ref, above=p_ref, pad=self.pad, align="left", force=True,
                )

        # --- bottom marginal: meta data (per left-cohort column group) ---
        if self.meta_data_position == "bottom":
            for j in range(M_left):
                self._grid_add(
                    grid_col_panel(left_name(ComutPanels.meta_data), j), left_widths[j], meta_h,
                    ref=blocks[(N - 1, j)], below=blocks[(N - 1, j)], align="left",
                )

        # --- left marginals (per gene group) ---
        left_panels = [ComutPanels.model_annotation, ComutPanels.gene_names, ComutPanels.cytoband, ComutPanels.gene_meta_data]
        for i in range(N):
            p_ref = blocks[(i, 0)]
            hi = self.gene_group_heights[i]
            for base in left_panels:
                w, _ = self._fixed_dim(base)
                if base == ComutPanels.cytoband:
                    pad = 0 if ComutPanels.gene_names in self.panels_to_plot else self.pad
                elif base == ComutPanels.gene_meta_data:
                    pad = self.pad if ComutPanels.cytoband in self.panels_to_plot else 0
                else:
                    pad = 0
                p_ref = self._grid_add(grid_row_panel(base, i), w, hi, ref=p_ref, left_of=p_ref, pad=pad, align="top")
            # gene-group label in the far-left gutter (only for labelled groups,
            # and only when the ``gene_group_label`` panel is requested)
            if self.grid.gene_groups[i].label and ComutPanels.gene_group_label in self.panels_to_plot:
                gw, _ = self._fixed_dim(ComutPanels.gene_group_label)
                self._grid_add(
                    grid_row_panel(ComutPanels.gene_group_label, i), gw, hi,
                    ref=p_ref, left_of=p_ref, pad=self.pad, align="top", force=True,
                )

        # --- right marginals (per gene group): left-cohort recurrence + fold change
        # + right-cohort recurrence (all full-cohort) ---
        right_refs = {}
        for i in range(N):
            p_ref = blocks[(i, M_left - 1)]
            hi = self.gene_group_heights[i]
            for j, base in enumerate([left_name(ComutPanels.recurrence), ComutPanels.recurrence_fold_change, right_name(ComutPanels.recurrence)]):
                pad = self.pad if draw_left_comut or j > 0 else 0
                w, _ = self._fixed_dim(base)
                p_ref = self._grid_add(grid_row_panel(base, i), w, hi, ref=p_ref, right_of=p_ref, pad=pad, align="top")
            right_refs[i] = p_ref

        # --- total recurrence summary strips: appended once at the very bottom,
        # beneath the LAST gene-group row's recurrence columns (same 1-row strips
        # as the non-grid layout). ---
        last = N - 1
        for base, total_base in [
            (left_name(ComutPanels.recurrence), left_name(ComutPanels.total_recurrence_overall)),
            (ComutPanels.recurrence_fold_change, ComutPanels.total_recurrence_fold_change),
            (right_name(ComutPanels.recurrence), right_name(ComutPanels.total_recurrence_overall)),
        ]:
            parent = self.panels.get(grid_row_panel(base, last))
            if parent is None:
                continue
            tw, th = self._fixed_dim(total_base)
            self._grid_add(
                grid_row_panel(total_base, last), tw, th,
                ref=parent, below=parent, pad=0, align="left",
            )

        # --- right cohort mirror (shared gene rows, own column groups) ---
        rightmost = dict(right_refs)
        if right_widths:
            draw_right_comut = right_name(ComutPanels.comutation) in self.panels_to_plot
            right_col_pad = self.pad if draw_right_comut else 0
            Mc = len(right_widths)
            cblocks = {}
            for i in range(N):
                hi = self.gene_group_heights[i]
                p_row_ref = right_refs[i]
                for j in range(Mc):
                    wj = right_widths[j] if draw_right_comut else 0
                    p = Layout.add_panel(
                        self, name=block(right_name(ComutPanels.comutation), i, j), width=wj, height=hi,
                        right_of=p_row_ref, pad=right_col_pad, align="top",
                    )
                    cblocks[(i, j)] = p
                    p_row_ref = p
                rightmost[i] = p_row_ref
            for j in range(Mc):
                p_ref = cblocks[(0, j)]
                wj = right_widths[j]
                if self.meta_data_position == "top":
                    p_ref = self._grid_add(
                        grid_col_panel(right_name(ComutPanels.meta_data), j), wj, meta_h,
                        ref=p_ref, above=p_ref, align="left",
                    )
                for base in top_panels:
                    _, h = self._fixed_dim(base)
                    p_ref = self._grid_add(grid_col_panel(right_name(base), j), wj, h, ref=p_ref, above=p_ref, align="left")
                label = right_labels[j] if j < len(right_labels) else ""
                if label and ComutPanels.column_group_label in self.panels_to_plot:
                    _, lh = self._fixed_dim(right_name(ComutPanels.column_group_label))
                    self._grid_add(
                        grid_col_panel(right_name(ComutPanels.column_group_label), j), wj, lh,
                        ref=p_ref, above=p_ref, pad=self.pad, align="left", force=True,
                    )
            if self.meta_data_position == "bottom":
                for j in range(Mc):
                    self._grid_add(
                        grid_col_panel(right_name(ComutPanels.meta_data), j), right_widths[j], meta_h,
                        ref=cblocks[(N - 1, j)], below=cblocks[(N - 1, j)], align="left",
                    )

        # --- legends (global, stacked to the right of the top row) ---
        self._add_grid_legends(ref=rightmost.get(0, blocks[(0, M_left - 1)]))

        self.place_panels_on_gridspec(autoscale_figsize=True, scale=0.1)

    def _add_grid_legends(self, ref):
        legend_panels = [
            ComutPanels.tmb_legend, ComutPanels.mutational_signatures_legend,
            ComutPanels.cnv_legend, ComutPanels.snv_legend, ComutPanels.model_annotation_legend,
        ]
        names = [b for b in legend_panels if b in self.panels_to_plot]
        p_ref = None
        for name in names:
            w, h = self._fixed_dim(name)
            if p_ref is None:
                p_ref = Layout.add_panel(self, name=name, width=w, height=h, right_of=ref, pad=self.pad, align="top")
            else:
                p_ref = Layout.add_panel(self, name=name, width=w, height=h, below=p_ref, pad=self.inter_legend_height, align="left")

        # Meta-data legends sit below the bottom metadata marginal when metadata
        # is at the bottom. If metadata is moved to the top, keep its legends
        # below the comutation area so they do not separate metadata from comut.
        if ComutPanels.meta_data_legend in self.panels_to_plot and self.meta_data_legend_titles:
            anchor = None
            if self.meta_data_position == "bottom":
                anchor = self.panels.get(grid_col_panel(ComutPanels.meta_data, 0))
            if anchor is None:
                anchor = self.panels.get(block(ComutPanels.comutation, self.grid.n_rows - 1, 0))
            pad = self.pad + (self.column_names_height if self.show_patient_names else 0)
            m_ref = None
            for title in self.meta_data_legend_titles:
                name = meta_legend(title)
                w, h = self._fixed_dim(name)
                if m_ref is None:
                    m_ref = Layout.add_panel(self, name=name, width=w, height=h, below=anchor, pad=pad, align="left")
                else:
                    m_ref = Layout.add_panel(self, name=name, width=w, height=h, right_of=m_ref, pad=self.inter_legend_width, align="top")

    def set_plot_func(self, panel: str, plot_func, *args, **kwargs):
        if panel in self.panels:
            self.panels[panel].set_plot_func(plot_func, *args, **kwargs)
