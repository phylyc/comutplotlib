"""Grid partition model for splitting a comut figure into sub-panels.

This module is intentionally dependency-light (only ``pandas`` for the ``Index``
convenience) and free of any matplotlib / I/O concerns so it can be unit-tested
in isolation. It is the data backbone for the "grid sub-panels" feature described
in ``GRID_SUBPANELS_PLAN.md``.

Two orthogonal stratifications are supported, each optional:

* **column groups** — an ordered partition of the samples/patients (the x axis),
  typically derived from one or more metadata keys.
* **gene groups** — an ordered partition of the genes (the y axis), typically
  derived from explicit gene sets.

When there is exactly one column group and one gene group the partition is
*trivial* and the rest of the pipeline collapses to the legacy single-panel
behaviour.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Iterable, Iterator, Sequence

import pandas as pd


@dataclass(frozen=True)
class Group:
    """An ordered, labelled subset of one axis (genes or samples).

    Attributes:
        key: Stable identifier used to build panel-id suffixes. Should be safe to
            embed in a string id (no reliance on human formatting).
        label: Human-readable title drawn on the figure.
        members: Ordered members (gene names or sample/patient ids).
    """

    key: str
    label: str
    members: tuple[str, ...]

    def __len__(self) -> int:
        return len(self.members)

    def with_members(self, members: Iterable[str]) -> "Group":
        """Return a copy of this group with a new (e.g. re-sorted) member order."""
        return Group(key=self.key, label=self.label, members=tuple(members))


@dataclass(frozen=True)
class GridPartition:
    """A 2-D partition of the comut matrix into ``n_rows x n_cols`` blocks."""

    gene_groups: tuple[Group, ...]
    column_groups: tuple[Group, ...]

    @property
    def n_rows(self) -> int:
        return len(self.gene_groups)

    @property
    def n_cols(self) -> int:
        return len(self.column_groups)

    def is_trivial(self) -> bool:
        """True when there is a single gene group and a single column group."""
        return self.n_rows == 1 and self.n_cols == 1

    @property
    def reference_gene_group(self) -> Group:
        """The first grid row — evidence for column sorting (config-first)."""
        return self.gene_groups[0]

    @property
    def reference_column_group(self) -> Group:
        """The first grid column — evidence for gene sorting (config-first)."""
        return self.column_groups[0]

    def blocks(self) -> Iterator[tuple[int, int, Group, Group]]:
        """Yield ``(row, col, gene_group, column_group)`` for every block."""
        for i, gene_group in enumerate(self.gene_groups):
            for j, column_group in enumerate(self.column_groups):
                yield i, j, gene_group, column_group

    def all_genes(self, name: str | None = None) -> pd.Index:
        """Concatenated gene order across all gene groups (group order preserved)."""
        genes = [g for group in self.gene_groups for g in group.members]
        return pd.Index(genes, name=name)

    def all_columns(self, name: str | None = None) -> pd.Index:
        """Concatenated sample order across all column groups (group order preserved)."""
        columns = [c for group in self.column_groups for c in group.members]
        return pd.Index(columns, name=name)

    def with_gene_groups(self, gene_groups: Sequence[Group]) -> "GridPartition":
        return GridPartition(gene_groups=tuple(gene_groups), column_groups=self.column_groups)

    def with_column_groups(self, column_groups: Sequence[Group]) -> "GridPartition":
        return GridPartition(gene_groups=self.gene_groups, column_groups=tuple(column_groups))


def trivial_partition(genes: Iterable[str], columns: Iterable[str]) -> GridPartition:
    """Build a 1x1 partition that reproduces the legacy single-panel figure."""
    genes = tuple(genes)
    columns = tuple(columns)
    return GridPartition(
        gene_groups=(Group(key="", label="", members=genes),),
        column_groups=(Group(key="", label="", members=columns),),
    )


def apply_column_group_labels(
    groups: Sequence[Group],
    labels: Sequence[str] | Mapping[str, str] | None,
) -> tuple[Group, ...]:
    """Override the *display* labels of ``groups`` (their ``key`` is untouched).

    Args:
        groups: Ordered column groups (grid-column order).
        labels: Either a mapping ``group key -> label`` (keys that are not
            present are ignored) or a flat sequence applied positionally in
            grid-column order. A shorter sequence leaves trailing groups with
            their default label; empty entries (``""``/``None``) also keep the
            default label of that group.

    Returns:
        The relabelled groups (a no-op copy when ``labels`` is empty/None).
    """
    groups = tuple(groups)
    if not labels:
        return groups
    if isinstance(labels, Mapping):
        return tuple(
            Group(key=g.key, label=str(labels[g.key]), members=g.members)
            if g.key in labels and labels[g.key] is not None
            else g
            for g in groups
        )
    labels = list(labels)
    return tuple(
        Group(key=g.key, label=str(labels[i]), members=g.members)
        if i < len(labels) and labels[i] not in (None, "")
        else g
        for i, g in enumerate(groups)
    )


def build_column_groups(
    columns: Iterable[str],
    key_frame: pd.DataFrame | None,
    keys: Sequence[str],
    order: Sequence[str] | None = None,
    na_label: str = "NA",
    separator: str = " | ",
    labels: Sequence[str] | Mapping[str, str] | None = None,
) -> tuple[Group, ...]:
    """Partition samples into column groups by one or more metadata keys.

    The grid columns are ordered **hierarchically** following the order of
    ``keys`` (i.e. of ``--column-group-by``): samples are grouped/sorted first by
    the first key's value, then by the second key's value, and so on. This keeps
    all samples sharing a first-key value contiguous.

    Within each hierarchy level, values are ordered by (decreasing priority):

        1. the explicit ``order`` — a flat list of individual metadata *values*
           (e.g. ``["neg", "pos"]``), applied at *every* level,
        2. the number of samples carrying that value at that level (descending),
        3. the value itself (lexicographic ascending).

    Args:
        columns: Ordered sample/patient ids (defines within-group member order).
        key_frame: DataFrame indexed by sample id whose columns include ``keys``.
            Cell values may be scalars or single-element lists (per-sample rows);
            list values use their first element. ``None`` or missing keys map
            every sample to ``na_label``.
        keys: Metadata column names to stratify by. Their order defines the
            nesting order of the grid columns (first key = outermost split).
        order: Optional explicit ordering of individual metadata *values* (not
            joined keys), applied at every hierarchy level (highest priority).
        na_label: Group label/key fragment used for missing/NaN/"unknown" values.
        separator: Joins multiple key values into a single group label/key.
        labels: Optional explicit display labels overriding the auto-generated
            ones, either positionally (grid-column order) or as a
            ``group key -> label`` mapping. See :func:`apply_column_group_labels`.

    Returns:
        Ordered tuple of column :class:`Group` s. When ``keys`` is empty a single
        group containing all ``columns`` is returned.
    """
    columns = list(columns)
    if not keys:
        return (Group(key="", label="", members=tuple(columns)),)

    present_keys = [k for k in keys if key_frame is not None and k in key_frame.columns]

    def value_for(sample: str, key: str) -> str:
        if key_frame is None or sample not in key_frame.index:
            return na_label
        value = key_frame.at[sample, key]
        if isinstance(value, (list, tuple)):
            value = value[0] if len(value) else None
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return na_label
        text = str(value).strip()
        return text if text and text.lower() not in {"nan", "unknown"} else na_label

    def parts_for(sample: str) -> tuple[str, ...]:
        if present_keys:
            return tuple(value_for(sample, key) for key in present_keys)
        return (na_label,)

    # Group samples by their per-key value tuple (order-preserving members).
    raw_groups: dict[tuple[str, ...], list[str]] = {}
    for sample in columns:
        raw_groups.setdefault(parts_for(sample), []).append(sample)

    if not raw_groups:
        # No samples to partition (e.g. an empty cohort) -> single empty group.
        return (Group(key="", label="", members=()),)

    # Per-level value frequencies for the size-based fallback (Decision 2:
    # explicit order > size > name), computed independently for each level.
    n_levels = len(present_keys) if present_keys else 1
    level_counts: list[dict[str, int]] = [defaultdict(int) for _ in range(n_levels)]
    for parts, members in raw_groups.items():
        for level, value in enumerate(parts):
            level_counts[level][value] += len(members)

    order_rank = {value: rank for rank, value in enumerate(order or [])}
    default_rank = len(order_rank)

    def value_key(level: int, value: str) -> tuple[int, int, str]:
        return (order_rank.get(value, default_rank), -level_counts[level][value], value)

    def group_sort_key(item: tuple[tuple[str, ...], list[str]]) -> tuple:
        parts, _ = item
        # Nesting: compare level 0 first, then level 1, ... (hierarchical).
        return tuple(value_key(level, value) for level, value in enumerate(parts))

    ordered = sorted(raw_groups.items(), key=group_sort_key)
    groups = tuple(
        Group(key=separator.join(parts), label=separator.join(parts), members=tuple(members))
        for parts, members in ordered
    )
    return apply_column_group_labels(groups, labels)


def build_gene_groups(
    genes: Iterable[str],
    gene_sets: Mapping[str, Sequence[str]] | None,
    order: Sequence[str] | None = None,
    other_label: str = "Other",
    drop_ungrouped: bool = False,
) -> tuple[Group, ...]:
    """Partition genes into gene groups by explicit gene sets.

    A gene may only belong to the first group (in resolved order) that lists it.
    Genes not covered by any set go into a trailing ``other_label`` group unless
    ``drop_ungrouped`` is set.

    Args:
        genes: Ordered gene names (defines within-group member order).
        gene_sets: Mapping of group label -> gene names. ``None``/empty yields a
            single group containing all ``genes``. The **insertion order of the
            mapping determines the row order** of the grid (top to bottom).
        order: Optional explicit group-name order (highest priority); names not
            listed keep their ``gene_sets`` insertion order after the listed ones.
        other_label: Label/key for the trailing ungrouped genes.
        drop_ungrouped: If True, ungrouped genes are dropped instead of collected.

    Returns:
        Ordered tuple of gene :class:`Group` s in ``order``-then-``gene_sets``
        insertion order. The ``other_label`` group, when present, is always placed
        last regardless of ordering.
    """
    genes = list(genes)
    if not gene_sets:
        return (Group(key="", label="", members=tuple(genes)),)

    gene_pos = {g: i for i, g in enumerate(genes)}
    assigned: set[str] = set()
    raw_groups: dict[str, list[str]] = {}
    labels: dict[str, str] = {}

    # Resolve which set claims each gene, honouring the requested group order so
    # that overlapping sets assign a gene to its highest-priority group. The
    # iteration order here becomes the final row order: explicit ``order`` first,
    # then the remaining groups in ``gene_sets`` insertion order.
    resolved_order = list(order) if order is not None else list(gene_sets.keys())
    resolved_order += [k for k in gene_sets.keys() if k not in resolved_order]

    for key in resolved_order:
        members = [g for g in gene_sets.get(key, []) if g in gene_pos and g not in assigned]
        members.sort(key=lambda g: gene_pos[g])
        if not members:
            continue
        raw_groups[key] = members
        labels[key] = key
        assigned.update(members)

    # Preserve the resolved (config/order) sequence for the row order rather than
    # re-sorting by group size.
    ordered = [
        Group(key=key, label=labels.get(key, key), members=tuple(members))
        for key, members in raw_groups.items()
    ]

    if not drop_ungrouped:
        leftover = [g for g in genes if g not in assigned]
        if leftover:
            ordered.append(Group(key=other_label, label=other_label, members=tuple(leftover)))

    return tuple(ordered)

