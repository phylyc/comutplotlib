"""Unit tests for grid tiling in ``ComutLayout`` (gridspec placement only).

These construct the layout directly with a synthetic :class:`GridPartition`,
avoiding the palette/plotter/data layers. Matplotlib runs headless (Agg).
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from comutplotlib.comut_layout import ComutLayout  # noqa: E402
from comutplotlib.comut_panels import (  # noqa: E402
    ComutPanels,
    block,
    control,
    grid_col_panel,
    grid_row_panel,
)
from comutplotlib.grid import Group, GridPartition  # noqa: E402


def _group(key, n, offset=0):
    return Group(key=key, label=key, members=tuple(f"{key}_{offset + k}" for k in range(n)))


def _make_layout(gene_sizes, col_sizes, panels, control_sizes=None, **kwargs):
    gene_groups = tuple(_group(f"g{i}", n) for i, n in enumerate(gene_sizes))
    column_groups = tuple(_group(f"c{j}", n) for j, n in enumerate(col_sizes))
    grid = GridPartition(gene_groups=gene_groups, column_groups=column_groups)
    layout = ComutLayout(
        panels_to_plot=list(panels),
        n_genes=sum(gene_sizes),
        n_samples=sum(col_sizes),
        n_samples_control=sum(control_sizes) if control_sizes else 0,
        max_xfigsize_scale=1,
        grid=grid,
        control_column_group_sizes=control_sizes,
        **kwargs,
    )
    return layout


@pytest.fixture(autouse=True)
def _close_figs():
    yield
    plt.close("all")


# --- sizing ------------------------------------------------------------------

def test_grid_enabled_and_sizes():
    layout = _make_layout([3, 2], [4, 3, 3], [ComutPanels.comutation])
    assert layout.grid_enabled
    assert layout.gene_group_heights == [3, 2]
    assert layout.col_group_widths == [4, 3, 3]


def test_trivial_grid_uses_legacy_path():
    gene_groups = (Group("", "", ("A", "B")),)
    column_groups = (Group("", "", ("s1", "s2", "s3")),)
    grid = GridPartition(gene_groups=gene_groups, column_groups=column_groups)
    layout = ComutLayout(
        panels_to_plot=[ComutPanels.comutation],
        n_genes=2, n_samples=3, grid=grid, max_xfigsize_scale=1,
    )
    assert not layout.grid_enabled
    layout.add_panels()
    # legacy path places the bare "comutation" panel, not namespaced blocks
    assert ComutPanels.comutation in layout.panels


def test_legacy_meta_data_defaults_below_comutation():
    panels = [ComutPanels.comutation, ComutPanels.mutational_signatures, ComutPanels.meta_data]
    layout = ComutLayout(
        panels_to_plot=panels,
        n_genes=3,
        n_samples=4,
        n_meta=2,
        max_xfigsize_scale=1,
    )
    layout.add_panels()

    comut = layout.panels[ComutPanels.comutation]
    meta = layout.panels[ComutPanels.meta_data]
    signatures = layout.panels[ComutPanels.mutational_signatures]
    assert signatures.y < comut.y < meta.y


def test_legacy_meta_data_top_sits_between_signatures_and_comutation():
    panels = [ComutPanels.comutation, ComutPanels.mutational_signatures, ComutPanels.meta_data]
    layout = ComutLayout(
        panels_to_plot=panels,
        n_genes=3,
        n_samples=4,
        n_meta=2,
        max_xfigsize_scale=1,
        meta_data_position="top",
    )
    layout.add_panels()

    comut = layout.panels[ComutPanels.comutation]
    meta = layout.panels[ComutPanels.meta_data]
    signatures = layout.panels[ComutPanels.mutational_signatures]
    assert signatures.y < meta.y < comut.y
    assert meta.x == comut.x


# --- block tiling ------------------------------------------------------------

def test_all_blocks_placed():
    layout = _make_layout([3, 2], [4, 3, 3], [ComutPanels.comutation])
    layout.add_panels()
    for i in range(2):
        for j in range(3):
            assert block(ComutPanels.comutation, i, j) in layout.panels


def test_block_dimensions_match_group_sizes():
    layout = _make_layout([3, 2], [4, 3, 3], [ComutPanels.comutation])
    layout.add_panels()
    b00 = layout.panels[block(ComutPanels.comutation, 0, 0)]
    b12 = layout.panels[block(ComutPanels.comutation, 1, 2)]
    assert (b00.width, b00.height) == (4, 3)
    assert (b12.width, b12.height) == (3, 2)


def test_block_relative_offsets_include_padding():
    pad = 1
    layout = _make_layout([3, 2], [4, 3, 3], [ComutPanels.comutation])
    layout.add_panels()
    b00 = layout.panels[block(ComutPanels.comutation, 0, 0)]
    b01 = layout.panels[block(ComutPanels.comutation, 0, 1)]
    b10 = layout.panels[block(ComutPanels.comutation, 1, 0)]
    # next column starts after width(4) + pad
    assert b01.x - b00.x == 4 + pad
    # next gene-group row starts after height(3) + pad
    assert b10.y - b00.y == 3 + pad


# --- marginals ---------------------------------------------------------------

def test_top_marginals_are_per_column_group():
    panels = [ComutPanels.comutation, ComutPanels.tmb]
    layout = _make_layout([3, 2], [4, 3, 3], panels)
    layout.add_panels()
    for j in range(3):
        assert grid_col_panel(ComutPanels.tmb, j) in layout.panels
    # a top marginal aligns to its block width
    tmb0 = layout.panels[grid_col_panel(ComutPanels.tmb, 0)]
    assert tmb0.width == 4


def test_grid_meta_data_top_sits_between_signatures_and_comutation_blocks():
    panels = [ComutPanels.comutation, ComutPanels.mutational_signatures, ComutPanels.meta_data]
    layout = _make_layout([3, 2], [4, 3], panels, meta_data_position="top", n_meta=2)
    layout.add_panels()

    block00 = layout.panels[block(ComutPanels.comutation, 0, 0)]
    meta0 = layout.panels[grid_col_panel(ComutPanels.meta_data, 0)]
    sig0 = layout.panels[grid_col_panel(ComutPanels.mutational_signatures, 0)]
    assert sig0.y < meta0.y < block00.y
    assert meta0.x == block00.x


def test_grid_meta_data_default_remains_below_comutation_blocks():
    panels = [ComutPanels.comutation, ComutPanels.mutational_signatures, ComutPanels.meta_data]
    layout = _make_layout([3, 2], [4, 3], panels, n_meta=2)
    layout.add_panels()

    block_last = layout.panels[block(ComutPanels.comutation, 1, 0)]
    meta0 = layout.panels[grid_col_panel(ComutPanels.meta_data, 0)]
    sig0 = layout.panels[grid_col_panel(ComutPanels.mutational_signatures, 0)]
    assert sig0.y < block_last.y < meta0.y


def test_left_marginals_are_per_gene_group():
    panels = [ComutPanels.comutation, ComutPanels.gene_names]
    layout = _make_layout([3, 2], [4, 3], panels)
    layout.add_panels()
    for i in range(2):
        assert grid_row_panel(ComutPanels.gene_names, i) in layout.panels
    gn1 = layout.panels[grid_row_panel(ComutPanels.gene_names, 1)]
    assert gn1.height == 2  # matches gene group 1 size


def test_grid_cytoband_uses_padding_when_gene_names_absent():
    panels = [ComutPanels.comutation, ComutPanels.cytoband]
    layout = _make_layout([3, 2], [4, 3], panels)
    layout.add_panels()

    for i in range(2):
        cytoband = layout.panels[grid_row_panel(ComutPanels.cytoband, i)]
        comut = layout.panels[block(ComutPanels.comutation, i, 0)]
        assert comut.x - (cytoband.x + cytoband.width) == layout.pad


def test_grid_cytoband_has_no_padding_when_gene_names_present():
    panels = [ComutPanels.comutation, ComutPanels.gene_names, ComutPanels.cytoband]
    layout = _make_layout([3, 2], [4, 3], panels)
    layout.add_panels()

    for i in range(2):
        cytoband = layout.panels[grid_row_panel(ComutPanels.cytoband, i)]
        gene_names = layout.panels[grid_row_panel(ComutPanels.gene_names, i)]
        assert gene_names.x - (cytoband.x + cytoband.width) == 0


def test_grid_gene_meta_data_uses_padding_when_cytoband_present():
    panels = [ComutPanels.comutation, ComutPanels.cytoband, ComutPanels.gene_meta_data]
    layout = _make_layout([3, 2], [4, 3], panels, n_meta_genes=2)
    layout.add_panels()

    for i in range(2):
        gene_meta = layout.panels[grid_row_panel(ComutPanels.gene_meta_data, i)]
        cytoband = layout.panels[grid_row_panel(ComutPanels.cytoband, i)]
        assert cytoband.x - (gene_meta.x + gene_meta.width) == layout.pad


def test_grid_gene_meta_data_has_no_padding_when_cytoband_absent():
    panels = [ComutPanels.comutation, ComutPanels.gene_meta_data]
    layout = _make_layout([3, 2], [4, 3], panels, n_meta_genes=2)
    layout.add_panels()

    for i in range(2):
        gene_meta = layout.panels[grid_row_panel(ComutPanels.gene_meta_data, i)]
        comut = layout.panels[block(ComutPanels.comutation, i, 0)]
        assert comut.x - (gene_meta.x + gene_meta.width) == 0


def test_right_marginals_recurrence_per_gene_group():
    panels = [ComutPanels.comutation, ComutPanels.recurrence]
    layout = _make_layout([3, 2], [4, 3], panels)
    layout.add_panels()
    for i in range(2):
        assert grid_row_panel(ComutPanels.recurrence, i) in layout.panels


def test_unselected_panels_not_placed():
    panels = [ComutPanels.comutation]  # no tmb/recurrence/etc.
    layout = _make_layout([3, 2], [4, 3], panels)
    layout.add_panels()
    assert grid_col_panel(ComutPanels.tmb, 0) not in layout.panels
    assert grid_row_panel(ComutPanels.recurrence, 0) not in layout.panels


# --- control mirror ----------------------------------------------------------

def test_control_blocks_placed():
    panels = [ComutPanels.comutation, control(ComutPanels.comutation)]
    layout = _make_layout([3, 2], [4, 3], panels, control_sizes=[5, 5])
    layout.add_panels()
    assert layout.col_group_widths_control == [5, 5]
    for i in range(2):
        for j in range(2):
            assert block(control(ComutPanels.comutation), i, j) in layout.panels


def test_control_blocks_right_of_case():
    panels = [ComutPanels.comutation, control(ComutPanels.comutation)]
    layout = _make_layout([3], [4, 3], panels, control_sizes=[6])
    layout.add_panels()
    case = layout.panels[block(ComutPanels.comutation, 0, 0)]
    ctrl = layout.panels[block(control(ComutPanels.comutation), 0, 0)]
    assert ctrl.x > case.x


# --- control position (left/right mirror) ------------------------------------

def test_legacy_control_position_left_swaps_cohorts():
    panels = [ComutPanels.comutation, control(ComutPanels.comutation), ComutPanels.recurrence_fold_change]
    right = ComutLayout(
        panels_to_plot=panels, n_genes=3, n_samples=4, n_samples_control=4, max_xfigsize_scale=1,
    )
    right.add_panels()
    left = ComutLayout(
        panels_to_plot=panels, n_genes=3, n_samples=4, n_samples_control=4, max_xfigsize_scale=1,
        control_position="left",
    )
    left.add_panels()
    # default: control is right of case; left: control is left of case
    assert right.panels[control(ComutPanels.comutation)].x > right.panels[ComutPanels.comutation].x
    assert left.panels[control(ComutPanels.comutation)].x < left.panels[ComutPanels.comutation].x


def test_control_position_left_inert_without_control():
    layout = ComutLayout(
        panels_to_plot=[ComutPanels.comutation], n_genes=3, n_samples=4, max_xfigsize_scale=1,
        control_position="left",
    )
    assert layout.control_on_left is False


def test_grid_control_blocks_left_of_case_when_left():
    panels = [ComutPanels.comutation, control(ComutPanels.comutation)]
    layout = _make_layout([3], [4, 3], panels, control_sizes=[6], control_position="left")
    layout.add_panels()
    case = layout.panels[block(ComutPanels.comutation, 0, 0)]
    ctrl = layout.panels[block(control(ComutPanels.comutation), 0, 0)]
    assert ctrl.x < case.x


# --- gridspec integrity ------------------------------------------------------

def test_gridspec_builds_without_error():
    panels = [
        ComutPanels.comutation, ComutPanels.tmb, ComutPanels.gene_names,
        ComutPanels.recurrence, ComutPanels.meta_data, ComutPanels.snv_legend,
        ComutPanels.cnv_legend,
    ]
    layout = _make_layout([3, 2], [4, 3, 3], panels)
    layout.add_panels()
    assert layout.fig is not None
    assert layout.gs is not None
    # every placed panel got an axis
    placed = [p for p in layout.panels.values() if p.width and p.height]
    assert all(p.ax is not None for p in placed)


