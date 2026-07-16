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

