"""End-to-end integration test for grid rendering via the ``Comut`` pipeline.

Uses the committed demo inputs and renders headless (Agg).
"""

import os

import matplotlib

matplotlib.use("Agg")
import pytest  # noqa: E402

DEMO = os.path.join(os.path.dirname(__file__), "..", "demo")


def _has_demo() -> bool:
    return os.path.exists(os.path.join(DEMO, "test.maf.tsv"))


pytestmark = pytest.mark.skipif(not _has_demo(), reason="demo inputs not available")



def _run(output, **extra):
    from comutplotlib.comut import Comut
    from comutplotlib.comut_panels import ComutPanels

    meta_data_rows = extra.pop("meta_data_rows", ["Sample Type"])
    panels = [
        ComutPanels.comutation, ComutPanels.tmb, ComutPanels.gene_names,
        ComutPanels.cytoband, ComutPanels.recurrence, ComutPanels.meta_data,
        ComutPanels.model_annotation, ComutPanels.snv_legend, ComutPanels.cnv_legend,
        ComutPanels.model_annotation_legend, ComutPanels.meta_data_legend,
        ComutPanels.column_group_label, ComutPanels.gene_group_label,
    ]
    return Comut(
        output=output,
        maf=[os.path.join(DEMO, "test.maf.tsv")],
        gistic=[os.path.join(DEMO, "test.all_thresholded.by_genes.txt")],
        sif=[os.path.join(DEMO, "test.sif.tsv")],
        meta_data_rows=meta_data_rows,
        snv_interesting_genes={f"Gene {i}" for i in range(1, 13)},
        cnv_interesting_genes={f"Gene {i}" for i in range(13, 21)},
        recurrence_categories={"global": ["max"]},
        panels_to_plot=panels,
        max_xfigsize=8,
        **extra,
    )


def test_hide_grouped_meta_data_preserves_grouping_and_other_rows(tmp_path):
    comut = _run(
        str(tmp_path / "hidden_group_meta.png"),
        meta_data_rows=["Sample Type", "Histology"],
        column_group_by=["Sample Type"],
        hide_grouped_meta_data=True,
    )
    comut.load()
    comut.preprocess()
    comut.build_layout()

    assert {group.key for group in comut.case.grid.column_groups} == {"BM", "EM"}
    for data in [comut.case, comut.control, comut.joint]:
        assert list(data.meta.rows) == ["Histology"]
        assert list(data.meta.df.columns) == ["Histology"]
    assert "Sample Type" not in comut.meta_cmaps
    assert list(comut.meta_cmaps) == ["Histology"]
    assert comut.layout.dimensions.loc["meta data", "height"] == 1


def test_grid_render_end_to_end(tmp_path):
    out = str(tmp_path / "grid.png")
    comut = _run(
        out,
        column_group_by=["Sample Type"],
        gene_groups={"SNV set": ["Gene 1", "Gene 2", "Gene 3", "Gene 4"],
                     "CNV set": ["Gene 13", "Gene 14", "Gene 15", "Gene 16"]},
    )
    genes, columns = comut.make_comut()
    assert os.path.exists(out)

    # grid is active and non-trivial
    assert comut.case.grid is not None and not comut.case.grid.is_trivial()
    keys = [g.key for g in comut.case.grid.gene_groups]
    assert keys[-1] == "Other"
    assert {"SNV set", "CNV set"}.issubset(set(keys))
    # column groups from Sample Type
    assert {g.key for g in comut.case.grid.column_groups} == {"BM", "EM"}

    # every gene-group's members are contiguous in the final gene axis
    gene_list = list(genes)
    for group in comut.case.grid.gene_groups:
        idx = [gene_list.index(m) for m in group.members]
        assert idx == list(range(min(idx), max(idx) + 1))

    # group-label panels are placed and rendered for every labelled group
    from comutplotlib.comut_panels import (
        ComutPanels, grid_col_panel, grid_row_panel,
    )
    for j in range(comut.case.grid.n_cols):
        name = grid_col_panel(ComutPanels.column_group_label, j)
        panel = comut.layout.panels.get(name)
        assert panel is not None and panel.ax is not None and panel.plot_func is not None
    for i in range(comut.case.grid.n_rows):
        name = grid_row_panel(ComutPanels.gene_group_label, i)
        panel = comut.layout.panels.get(name)
        assert panel is not None and panel.ax is not None and panel.plot_func is not None


def test_explicit_column_group_labels_are_rendered(tmp_path):
    """``column_group_labels`` replaces the auto-generated column-group titles."""
    from comutplotlib.comut_panels import ComutPanels, grid_col_panel

    out = str(tmp_path / "grid_labels.png")
    comut = _run(
        out,
        column_group_by=["Sample Type"],
        group_order=["BM", "EM"],
        column_group_labels=["Brain metastases", "Extracranial metastases"],
    )
    comut.make_comut()

    # keys keep the raw metadata values, labels are the user-supplied ones
    groups = comut.case.grid.column_groups
    assert {g.key for g in groups} == {"BM", "EM"}
    assert {g.key: g.label for g in groups} == {
        "BM": "Brain metastases", "EM": "Extracranial metastases",
    }
    # ... and they end up on the rendered header panels (text may be wrapped)
    for j, group in enumerate(groups):
        panel = comut.layout.panels.get(grid_col_panel(ComutPanels.column_group_label, j))
        assert panel is not None and panel.ax is not None
        drawn = [t.get_text().replace("\n", " ") for t in panel.ax.texts]
        assert drawn == [group.label]


def test_partial_column_group_labels_keep_defaults(tmp_path):
    out = str(tmp_path / "grid_partial_labels.png")
    comut = _run(
        out,
        column_group_by=["Sample Type"],
        group_order=["BM", "EM"],
        column_group_labels=["Brain metastases"],
    )
    comut.make_comut()
    assert {g.key: g.label for g in comut.case.grid.column_groups} == {
        "BM": "Brain metastases", "EM": "EM",
    }


def test_trivial_pipeline_still_renders(tmp_path):
    out = str(tmp_path / "plain.png")
    comut = _run(out)  # no grouping
    comut.make_comut()
    assert comut.case.grid is None or comut.case.grid.is_trivial()
    assert os.path.exists(out)


def test_single_cohort_ignores_fold_change_threshold(tmp_path):
    comut = _run(
        str(tmp_path / "single_cohort.png"),
        min_fold_change=1.05,
        gene_groups={"SNV set": ["Gene 1", "Gene 2"]},
    )
    comut.load()
    genes_before = list(comut.joint.genes)

    comut.preprocess()

    assert set(comut.case.genes) == set(genes_before)


def test_group_labels_only_when_requested(tmp_path):
    """The group-label panels are enumerated but not in DEFAULT_PANELS; they must
    only be placed when their strings are passed via panels_to_plot."""
    from comutplotlib.comut import Comut
    from comutplotlib.comut_panels import (
        ComutPanels, grid_col_panel, grid_row_panel,
    )

    # Same grid, but WITHOUT the label strings in panels_to_plot.
    panels = [
        ComutPanels.comutation, ComutPanels.tmb, ComutPanels.gene_names,
        ComutPanels.cytoband, ComutPanels.recurrence, ComutPanels.meta_data,
    ]
    comut = Comut(
        output=str(tmp_path / "no_labels.png"),
        maf=[os.path.join(DEMO, "test.maf.tsv")],
        gistic=[os.path.join(DEMO, "test.all_thresholded.by_genes.txt")],
        sif=[os.path.join(DEMO, "test.sif.tsv")],
        meta_data_rows=["Sample Type"],
        snv_interesting_genes={f"Gene {i}" for i in range(1, 13)},
        cnv_interesting_genes={f"Gene {i}" for i in range(13, 21)},
        recurrence_categories={"global": ["max"]},
        panels_to_plot=panels,
        max_xfigsize=8,
        column_group_by=["Sample Type"],
        gene_groups={"SNV set": ["Gene 1", "Gene 2"], "CNV set": ["Gene 13", "Gene 14"]},
    )
    comut.make_comut()
    assert comut.case.grid is not None and not comut.case.grid.is_trivial()
    # no group-label panel should have been placed
    for j in range(comut.case.grid.n_cols):
        assert comut.layout.panels.get(grid_col_panel(ComutPanels.column_group_label, j)) is None
    for i in range(comut.case.grid.n_rows):
        assert comut.layout.panels.get(grid_row_panel(ComutPanels.gene_group_label, i)) is None


def test_grid_rows_do_not_collapse_without_comutation(tmp_path):
    """Omitting the comutation panel must not collapse the gene-group rows: the
    central blocks are still placed as scaffold so every row keeps a distinct
    vertical position (regression for 'only the Other group is plotted')."""
    from comutplotlib.comut import Comut
    from comutplotlib.comut_panels import ComutPanels, grid_row_panel, block

    # Note: no "comutation" in panels_to_plot.
    panels = [ComutPanels.gene_names, ComutPanels.recurrence, ComutPanels.cytoband]
    comut = Comut(
        output=str(tmp_path / "no_comut.png"),
        maf=[os.path.join(DEMO, "test.maf.tsv")],
        gistic=[os.path.join(DEMO, "test.all_thresholded.by_genes.txt")],
        sif=[os.path.join(DEMO, "test.sif.tsv")],
        meta_data_rows=["Sample Type"],
        snv_interesting_genes={f"Gene {i}" for i in range(1, 13)},
        cnv_interesting_genes={f"Gene {i}" for i in range(13, 21)},
        recurrence_categories={"global": ["max"]},
        panels_to_plot=panels,
        max_xfigsize=8,
        gene_groups={"SNV set": ["Gene 1", "Gene 2", "Gene 3"],
                     "CNV set": ["Gene 13", "Gene 14", "Gene 15"]},
    )
    comut.make_comut()
    grid = comut.case.grid
    assert grid is not None and not grid.is_trivial()
    N = grid.n_rows
    assert N >= 3  # SNV set, CNV set, Other

    # central blocks are placed as scaffold even though comutation is off, but
    # with zero horizontal width so they don't add an extra empty column.
    for i in range(N):
        p = comut.layout.panels.get(block(ComutPanels.comutation, i, 0))
        assert p is not None, "scaffold block missing"
        assert p.width == 0, "scaffold block should take no horizontal space"

    # every gene-name row panel has a distinct top (no collapse to y=0)
    ys = []
    for i in range(N):
        p = comut.layout.panels.get(grid_row_panel(ComutPanels.gene_names, i))
        assert p is not None and p.ax is not None
        ys.append(p.y)
    assert len(set(ys)) == N, f"gene rows collapsed: y positions = {ys}"
    assert ys == sorted(ys) and ys[0] == 0

    # the gap between the gene-name column and the recurrence column is 0
    gn = comut.layout.panels.get(grid_row_panel(ComutPanels.gene_names, 0))
    rec = comut.layout.panels.get(grid_row_panel(ComutPanels.recurrence, 0))
    assert rec is not None and rec.ax is not None
    assert (rec.x - (gn.x + gn.width)) == 0


def test_grid_plots_control_recurrence_when_requested(tmp_path):
    """In grid mode with control data, requesting 'recurrence control' must place
    and plot a control-recurrence marginal per gene group (regression)."""
    from comutplotlib.comut import Comut
    from comutplotlib.comut_panels import ComutPanels, grid_row_panel, control

    panels = [
        ComutPanels.comutation, ComutPanels.gene_names,
        ComutPanels.recurrence, ComutPanels.recurrence_fold_change,
        control(ComutPanels.recurrence),
        ComutPanels.snv_legend, ComutPanels.cnv_legend,
    ]
    comut = Comut(
        output=str(tmp_path / "grid_rec_control.png"),
        maf=[os.path.join(DEMO, "test.maf.tsv")],
        gistic=[os.path.join(DEMO, "test.all_thresholded.by_genes.txt")],
        sif=[os.path.join(DEMO, "test.sif.tsv")],
        control_maf=[os.path.join(DEMO, "control.maf.tsv")],
        control_gistic=[os.path.join(DEMO, "control.all_thresholded.by_genes.txt")],
        control_sif=[os.path.join(DEMO, "control.sif.tsv")],
        snv_interesting_genes={f"Gene {i}" for i in range(1, 9)},
        cnv_interesting_genes={f"Gene {i}" for i in range(9, 17)},
        recurrence_categories={"global": ["max"]},
        panels_to_plot=panels,
        max_xfigsize=6,
        gene_groups={"SNV genes": ["Gene 1", "Gene 2", "Gene 3"],
                     "CNV genes": ["Gene 13", "Gene 14", "Gene 15"]},
    )
    comut.make_comut()
    grid = comut.case.grid
    assert grid is not None and not grid.is_trivial()
    # a control-recurrence marginal exists and is rendered for every gene group
    for i in range(grid.n_rows):
        p = comut.layout.panels.get(grid_row_panel(control(ComutPanels.recurrence), i))
        assert p is not None, "recurrence control marginal missing"
        assert p.ax is not None and p.plot_func is not None
        # it sits to the right of the case recurrence marginal
        rec = comut.layout.panels.get(grid_row_panel(ComutPanels.recurrence, i))
        assert p.x > rec.x


def test_grid_total_recurrence_strips_only_on_last_row(tmp_path):
    """The total-recurrence summary strips (overall, fold-change, overall-control)
    appear once, beneath the LAST gene-group row's recurrence columns."""
    from comutplotlib.comut import Comut
    from comutplotlib.comut_panels import ComutPanels, grid_row_panel, control

    panels = [
        ComutPanels.comutation, ComutPanels.gene_names,
        ComutPanels.recurrence, ComutPanels.recurrence_fold_change,
        control(ComutPanels.recurrence),
        ComutPanels.total_recurrence_overall,
        ComutPanels.total_recurrence_fold_change,
        control(ComutPanels.total_recurrence_overall),
        ComutPanels.snv_legend, ComutPanels.cnv_legend,
    ]
    comut = Comut(
        output=str(tmp_path / "grid_totals.png"),
        maf=[os.path.join(DEMO, "test.maf.tsv")],
        gistic=[os.path.join(DEMO, "test.all_thresholded.by_genes.txt")],
        sif=[os.path.join(DEMO, "test.sif.tsv")],
        control_maf=[os.path.join(DEMO, "control.maf.tsv")],
        control_gistic=[os.path.join(DEMO, "control.all_thresholded.by_genes.txt")],
        control_sif=[os.path.join(DEMO, "control.sif.tsv")],
        snv_interesting_genes={f"Gene {i}" for i in range(1, 9)},
        cnv_interesting_genes={f"Gene {i}" for i in range(9, 17)},
        recurrence_categories={"global": ["max"]},
        panels_to_plot=panels,
        max_xfigsize=6,
        gene_groups={"SNV genes": ["Gene 1", "Gene 2", "Gene 3"],
                     "CNV genes": ["Gene 13", "Gene 14", "Gene 15"]},
    )
    comut.make_comut()
    grid = comut.case.grid
    N = grid.n_rows
    assert N >= 3
    last = N - 1
    totals = [
        (ComutPanels.total_recurrence_overall, ComutPanels.recurrence),
        (ComutPanels.total_recurrence_fold_change, ComutPanels.recurrence_fold_change),
        (control(ComutPanels.total_recurrence_overall), control(ComutPanels.recurrence)),
    ]
    for total_base, parent_base in totals:
        # placed & rendered on the last row only
        for i in range(N):
            p = comut.layout.panels.get(grid_row_panel(total_base, i))
            if i == last:
                assert p is not None and p.ax is not None and p.plot_func is not None
                # a 1-tall strip directly beneath the parent recurrence column
                parent = comut.layout.panels.get(grid_row_panel(parent_base, last))
                assert p.height == 1
                assert p.x == parent.x
                assert p.y == parent.y + parent.height
            else:
                assert p is None, f"{total_base!r} must not appear on row {i}"
