"""Canonical panel identifiers for the comut figure — a single source of truth.

Panel ids are plain strings because they are used as layout/gridspec dict keys,
as CLI values (``--panels-to-plot``), and as f-string fragments (e.g.
``name + CONTROL_SUFFIX``). This module is the comut-specific counterpart to the
general-purpose :class:`comutplotlib.panel.Panel` class, mirroring the
Plotter/ComutPlotter and Layout/ComutLayout split.

Attribute names deliberately use the same lowercase style as the string
vocabularies in :class:`~comutplotlib.mutation_annotation.MutationAnnotation`
and :class:`~comutplotlib.sample_annotation.SampleAnnotation`.
"""


class ComutPanels:
    # internal anchor / reference panel (always placed)
    anchor = "_anchor"

    # center
    comutation = "comutation"

    # top panels
    cohort_label = "cohort label"
    coverage = "coverage"
    tmb = "tmb"
    mutational_signatures = "mutational signatures"

    # left panels (gene axis)
    gene_names = "gene names"
    cytoband = "cytoband"
    gene_meta_data = "gene meta data"
    model_annotation = "model annotation"
    model_significance = "model significance"

    # right panels (recurrence)
    recurrence = "recurrence"
    total_recurrence = "total recurrence"
    total_recurrence_overall = "total recurrence overall"
    recurrence_fold_change = "recurrence fold change"
    total_recurrence_fold_change = "total recurrence fold change"

    # bottom panels
    meta_data = "meta data"

    # grid group labels (only drawn when a non-trivial grid is active)
    gene_group_label = "gene group label"
    column_group_label = "column group label"

    # legends
    tmb_legend = "tmb legend"
    mutational_signatures_legend = "mutational signatures legend"
    snv_legend = "snv legend"
    cnv_legend = "cnv legend"
    model_annotation_legend = "model annotation legend"
    meta_data_legend = "meta data legend"


CONTROL_SUFFIX = " control"


def control(panel: str) -> str:
    """Return the control-cohort variant of a panel id (``"tmb" -> "tmb control"``)."""
    return panel + CONTROL_SUFFIX


def meta_legend(title: str) -> str:
    """Return the per-column metadata-legend panel id."""
    return f"{ComutPanels.meta_data_legend} {title}"


# --- Grid sub-panel ids ------------------------------------------------------
# When a figure is split into a grid (see ``GRID_SUBPANELS_PLAN.md``) every panel
# instance needs a unique id. These helpers namespace a base panel id with its
# grid row and/or column. They ALWAYS namespace; the legacy (single-panel) code
# path simply keeps using the bare base ids, so trivial-grid figures are
# unchanged. The separator is chosen to be unlikely to collide with the
# space-containing base ids in :class:`ComutPanels`.

GRID_SEP = " @ "


def block(panel: str, row: int, col: int) -> str:
    """Central comut block id, e.g. ``"comutation @ r0 @ c2"``."""
    return f"{panel}{GRID_SEP}r{row}{GRID_SEP}c{col}"


def grid_row_panel(panel: str, row: int) -> str:
    """Gene-axis (left/right) marginal id for one gene group, e.g. ``"recurrence @ r1"``."""
    return f"{panel}{GRID_SEP}r{row}"


def grid_col_panel(panel: str, col: int) -> str:
    """Column-axis (top/bottom) marginal id for one column group, e.g. ``"tmb @ c0"``."""
    return f"{panel}{GRID_SEP}c{col}"


def parse_grid_panel(panel_id: str) -> tuple[str, int | None, int | None]:
    """Decode a (possibly namespaced) panel id into ``(base, row, col)``.

    ``row``/``col`` are ``None`` when the corresponding suffix is absent. The
    ``base`` retains any ``control`` suffix (``"comutation control @ r0 @ c1"``
    -> ``("comutation control", 0, 1)``).
    """
    parts = panel_id.split(GRID_SEP)
    base = parts[0]
    row: int | None = None
    col: int | None = None
    for part in parts[1:]:
        if part.startswith("r") and part[1:].isdigit():
            row = int(part[1:])
        elif part.startswith("c") and part[1:].isdigit():
            col = int(part[1:])
    return base, row, col


# --- Panel axis registry -----------------------------------------------------
# Classifies each base panel by the axis it spans, so grid tiling can be derived
# generically instead of hand-wiring every panel:
#   "center"  -> tiled along BOTH axes (one instance per (row, col) block)
#   "columns" -> spans the sample axis; tiled along columns (one per column group)
#   "genes"   -> spans the gene axis; tiled along genes (one per gene group)
#   "none"    -> global (legends); never tiled
AXIS_CENTER = "center"
AXIS_COLUMNS = "columns"
AXIS_GENES = "genes"
AXIS_NONE = "none"

PANEL_AXIS = {
    ComutPanels.comutation: AXIS_CENTER,
    # top / bottom marginals (span the sample axis)
    ComutPanels.cohort_label: AXIS_COLUMNS,
    ComutPanels.coverage: AXIS_COLUMNS,
    ComutPanels.tmb: AXIS_COLUMNS,
    ComutPanels.mutational_signatures: AXIS_COLUMNS,
    ComutPanels.meta_data: AXIS_COLUMNS,
    # left / right marginals (span the gene axis)
    ComutPanels.gene_names: AXIS_GENES,
    ComutPanels.cytoband: AXIS_GENES,
    ComutPanels.gene_meta_data: AXIS_GENES,
    ComutPanels.model_annotation: AXIS_GENES,
    ComutPanels.model_significance: AXIS_GENES,
    ComutPanels.recurrence: AXIS_GENES,
    ComutPanels.total_recurrence: AXIS_GENES,
    ComutPanels.total_recurrence_overall: AXIS_GENES,
    ComutPanels.recurrence_fold_change: AXIS_GENES,
    ComutPanels.total_recurrence_fold_change: AXIS_GENES,
    # grid group labels tile like their axis
    ComutPanels.gene_group_label: AXIS_GENES,
    ComutPanels.column_group_label: AXIS_COLUMNS,
    # legends are global
    ComutPanels.tmb_legend: AXIS_NONE,
    ComutPanels.mutational_signatures_legend: AXIS_NONE,
    ComutPanels.snv_legend: AXIS_NONE,
    ComutPanels.cnv_legend: AXIS_NONE,
    ComutPanels.model_annotation_legend: AXIS_NONE,
    ComutPanels.meta_data_legend: AXIS_NONE,
}


def panel_axis(panel: str) -> str:
    """Return the tiling axis for a base panel id (``control`` suffix tolerated)."""
    base = panel[: -len(CONTROL_SUFFIX)] if panel.endswith(CONTROL_SUFFIX) else panel
    return PANEL_AXIS.get(base, AXIS_NONE)


# Default set of panels rendered by the CLI. Order is not significant here (the
# layout defines placement), but it is kept close to the visual layout for
# readability. Commented-out entries are intentionally retained as reminders of
# panels that are implemented but not enabled by default.
DEFAULT_PANELS = [
    # Left side
    ComutPanels.gene_meta_data,
    ComutPanels.cytoband,
    ComutPanels.gene_names,
    ComutPanels.model_annotation,
    # Top
    ComutPanels.cohort_label,
    control(ComutPanels.cohort_label),
    ComutPanels.tmb,
    control(ComutPanels.tmb),
    ComutPanels.tmb_legend,
    ComutPanels.mutational_signatures,
    control(ComutPanels.mutational_signatures),
    ComutPanels.mutational_signatures_legend,
    # Center
    ComutPanels.comutation,
    control(ComutPanels.comutation),
    ComutPanels.recurrence,
    control(ComutPanels.recurrence),
    ComutPanels.total_recurrence_overall,
    control(ComutPanels.total_recurrence_overall),
    ComutPanels.recurrence_fold_change,
    ComutPanels.total_recurrence_fold_change,
    # Bottom
    ComutPanels.meta_data,
    control(ComutPanels.meta_data),
    # Right
    ComutPanels.snv_legend,
    ComutPanels.cnv_legend,
    ComutPanels.model_annotation_legend,
    ComutPanels.meta_data_legend,
]

