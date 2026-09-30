"""Unit tests for panel-id helpers and the axis registry (``comut_panels``)."""

from comutplotlib.comut_panels import (
    AXIS_CENTER,
    AXIS_COLUMNS,
    AXIS_GENES,
    AXIS_NONE,
    ComutPanels,
    block,
    control,
    grid_col_panel,
    grid_row_panel,
    meta_legend,
    panel_axis,
    parse_grid_panel,
)


# --- legacy helpers unchanged ------------------------------------------------

def test_control_suffix_unchanged():
    assert control(ComutPanels.tmb) == "tmb control"


def test_meta_legend_unchanged():
    assert meta_legend("Histology") == "meta data legend Histology"


# --- grid id construction ----------------------------------------------------

def test_block_id_namespacing():
    assert block(ComutPanels.comutation, 0, 2) == "comutation @ r0 @ c2"


def test_grid_row_and_col_panel_ids():
    assert grid_row_panel(ComutPanels.recurrence, 1) == "recurrence @ r1"
    assert grid_col_panel(ComutPanels.tmb, 0) == "tmb @ c0"


# --- grid id parsing (round-trip) --------------------------------------------

def test_parse_block_roundtrip():
    assert parse_grid_panel(block(ComutPanels.comutation, 3, 4)) == ("comutation", 3, 4)


def test_parse_row_panel():
    assert parse_grid_panel(grid_row_panel(ComutPanels.gene_names, 2)) == ("gene names", 2, None)


def test_parse_col_panel():
    assert parse_grid_panel(grid_col_panel(ComutPanels.meta_data, 5)) == ("meta data", None, 5)


def test_parse_plain_base_id():
    assert parse_grid_panel(ComutPanels.comutation) == ("comutation", None, None)


def test_parse_retains_control_suffix():
    pid = block(control(ComutPanels.comutation), 0, 1)
    assert pid == "comutation control @ r0 @ c1"
    assert parse_grid_panel(pid) == ("comutation control", 0, 1)


# --- axis registry -----------------------------------------------------------

def test_panel_axis_center():
    assert panel_axis(ComutPanels.comutation) == AXIS_CENTER


def test_panel_axis_columns_family():
    for p in [ComutPanels.tmb, ComutPanels.coverage, ComutPanels.mutational_signatures, ComutPanels.cohort_label, ComutPanels.meta_data]:
        assert panel_axis(p) == AXIS_COLUMNS


def test_panel_axis_genes_family():
    for p in [ComutPanels.gene_names, ComutPanels.cytoband, ComutPanels.recurrence, ComutPanels.recurrence_fold_change, ComutPanels.model_annotation]:
        assert panel_axis(p) == AXIS_GENES


def test_panel_axis_legends_none():
    for p in [ComutPanels.tmb_legend, ComutPanels.snv_legend, ComutPanels.cnv_legend, ComutPanels.meta_data_legend]:
        assert panel_axis(p) == AXIS_NONE


def test_panel_axis_tolerates_control_suffix():
    assert panel_axis(control(ComutPanels.tmb)) == AXIS_COLUMNS
    assert panel_axis(control(ComutPanels.recurrence)) == AXIS_GENES

