"""Unit tests for the grid partition model (``comutplotlib.grid``)."""

import pandas as pd
import pytest

from comutplotlib.grid import (
    Group,
    GridPartition,
    apply_column_group_labels,
    build_column_groups,
    build_gene_groups,
    trivial_partition,
)


# --- trivial partition -------------------------------------------------------

def test_trivial_partition_is_trivial():
    part = trivial_partition(genes=["A", "B"], columns=["s1", "s2", "s3"])
    assert part.is_trivial()
    assert part.n_rows == 1
    assert part.n_cols == 1
    assert part.reference_gene_group.members == ("A", "B")
    assert part.reference_column_group.members == ("s1", "s2", "s3")


def test_trivial_partition_all_genes_and_columns_roundtrip():
    part = trivial_partition(genes=["A", "B"], columns=["s1", "s2"])
    assert list(part.all_genes(name="gene")) == ["A", "B"]
    assert part.all_genes(name="gene").name == "gene"
    assert list(part.all_columns()) == ["s1", "s2"]


def test_blocks_iteration_shape():
    part = GridPartition(
        gene_groups=(Group("g0", "G0", ("A",)), Group("g1", "G1", ("B",))),
        column_groups=(Group("c0", "C0", ("s1",)), Group("c1", "C1", ("s2",)), Group("c2", "C2", ("s3",))),
    )
    assert not part.is_trivial()
    assert part.n_rows == 2 and part.n_cols == 3
    coords = [(i, j) for i, j, _, _ in part.blocks()]
    assert coords == [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]


def test_all_genes_preserves_group_order():
    part = GridPartition(
        gene_groups=(Group("g0", "G0", ("A", "B")), Group("g1", "G1", ("C",))),
        column_groups=(Group("c0", "C0", ("s1",)),),
    )
    assert list(part.all_genes()) == ["A", "B", "C"]


# --- column groups -----------------------------------------------------------

def test_build_column_groups_no_keys_returns_single_group():
    groups = build_column_groups(columns=["s1", "s2"], key_frame=None, keys=[])
    assert len(groups) == 1
    assert groups[0].members == ("s1", "s2")


def test_build_column_groups_single_key_partitions_and_orders_by_size():
    frame = pd.DataFrame(
        {"Histology": ["IDC", "ILC", "IDC", "IDC"]},
        index=["s1", "s2", "s3", "s4"],
    )
    groups = build_column_groups(columns=["s1", "s2", "s3", "s4"], key_frame=frame, keys=["Histology"])
    # IDC has 3 members, ILC has 1 -> IDC first (size desc)
    assert [g.key for g in groups] == ["IDC", "ILC"]
    assert groups[0].members == ("s1", "s3", "s4")
    assert groups[1].members == ("s2",)


def test_build_column_groups_explicit_order_wins_over_size():
    frame = pd.DataFrame(
        {"Histology": ["IDC", "ILC", "IDC", "IDC"]},
        index=["s1", "s2", "s3", "s4"],
    )
    groups = build_column_groups(
        columns=["s1", "s2", "s3", "s4"], key_frame=frame, keys=["Histology"], order=["ILC", "IDC"]
    )
    assert [g.key for g in groups] == ["ILC", "IDC"]


def test_build_column_groups_missing_values_go_to_na_label():
    frame = pd.DataFrame(
        {"Histology": ["IDC", None, "nan", "unknown"]},
        index=["s1", "s2", "s3", "s4"],
    )
    groups = build_column_groups(columns=["s1", "s2", "s3", "s4"], key_frame=frame, keys=["Histology"])
    na = {g.key: g for g in groups}["NA"]
    assert set(na.members) == {"s2", "s3", "s4"}


def test_build_column_groups_multi_key_tuple_labels_and_tiebreak():
    frame = pd.DataFrame(
        {"Histology": ["IDC", "IDC", "ILC"], "Platform": ["WES", "WGS", "WES"]},
        index=["s1", "s2", "s3"],
    )
    groups = build_column_groups(
        columns=["s1", "s2", "s3"], key_frame=frame, keys=["Histology", "Platform"]
    )
    # all groups size 1 -> ordered by key tuple lexicographically
    assert [g.key for g in groups] == ["IDC | WES", "IDC | WGS", "ILC | WES"]


def test_build_column_groups_multi_key_is_hierarchical():
    # HR (first key) should be the outer split: all HR=neg groups precede all
    # HR=pos groups, regardless of individual combination sizes.
    frame = pd.DataFrame(
        {
            "HR": ["pos", "pos", "pos", "neg", "neg"],
            "HER2": ["neg", "neg", "pos", "neg", "pos"],
        },
        index=["s1", "s2", "s3", "s4", "s5"],
    )
    cols = ["s1", "s2", "s3", "s4", "s5"]
    groups = build_column_groups(cols, frame, keys=["HR", "HER2"], order=["neg", "pos"])
    # group_order applies at every level: HR=neg block first, then HR=pos block,
    # and within each block HER2=neg before HER2=pos.
    assert [g.key for g in groups] == ["neg | neg", "neg | pos", "pos | neg", "pos | pos"]


def test_build_column_groups_key_order_defines_hierarchy():
    frame = pd.DataFrame(
        {
            "HR": ["pos", "neg", "pos", "neg"],
            "HER2": ["neg", "neg", "pos", "pos"],
        },
        index=["s1", "s2", "s3", "s4"],
    )
    cols = ["s1", "s2", "s3", "s4"]
    hr_first = build_column_groups(cols, frame, keys=["HR", "HER2"], order=["neg", "pos"])
    her2_first = build_column_groups(cols, frame, keys=["HER2", "HR"], order=["neg", "pos"])
    # the outer (first) key value is contiguous in both orderings
    assert [g.key.split(" | ")[0] for g in hr_first] == ["neg", "neg", "pos", "pos"]
    assert [g.key.split(" | ")[0] for g in her2_first] == ["neg", "neg", "pos", "pos"]
    # but the sample partition differs, because the nesting hierarchy differs
    assert [g.members for g in hr_first] != [g.members for g in her2_first]


def test_build_column_groups_order_is_individual_values_not_joined_keys():
    frame = pd.DataFrame(
        {"Status": ["pos", "neg", "unknown", "pos"]},
        index=["s1", "s2", "s3", "s4"],
    )
    cols = ["s1", "s2", "s3", "s4"]
    groups = build_column_groups(cols, frame, keys=["Status"], order=["neg", "pos"])
    # neg first, then pos, then the (unlisted) NA group from "unknown"
    assert [g.key for g in groups] == ["neg", "pos", "NA"]


def test_build_column_groups_list_cell_uses_first_element():
    frame = pd.DataFrame({"Histology": [["IDC"], ["ILC"]]}, index=["s1", "s2"])
    groups = build_column_groups(columns=["s1", "s2"], key_frame=frame, keys=["Histology"])
    assert {g.key for g in groups} == {"IDC", "ILC"}


def test_build_column_groups_missing_key_column_all_na():
    frame = pd.DataFrame({"Other": [1, 2]}, index=["s1", "s2"])
    groups = build_column_groups(columns=["s1", "s2"], key_frame=frame, keys=["Histology"])
    assert len(groups) == 1
    assert groups[0].key == "NA"


def test_build_column_groups_empty_columns_does_not_crash():
    # An empty cohort with a stratification key must not raise (regression for a
    # StopIteration on next(iter(raw_groups))).
    frame = pd.DataFrame({"Histology": []})
    groups = build_column_groups(columns=[], key_frame=frame, keys=["Histology"])
    assert len(groups) == 1
    assert groups[0].members == ()


# --- column group labels -----------------------------------------------------

def _status_frame():
    return pd.DataFrame(
        {"Status": ["pos", "neg", "pos", "neg"]},
        index=["s1", "s2", "s3", "s4"],
    )


def test_column_group_labels_applied_positionally():
    groups = build_column_groups(
        columns=["s1", "s2", "s3", "s4"], key_frame=_status_frame(), keys=["Status"],
        order=["neg", "pos"], labels=["Negative", "Positive"],
    )
    assert [g.key for g in groups] == ["neg", "pos"]
    assert [g.label for g in groups] == ["Negative", "Positive"]
    # members are untouched by the relabelling
    assert [g.members for g in groups] == [("s2", "s4"), ("s1", "s3")]


def test_column_group_labels_partial_and_blank_keep_defaults():
    groups = build_column_groups(
        columns=["s1", "s2", "s3", "s4"], key_frame=_status_frame(), keys=["Status"],
        order=["neg", "pos"], labels=[""],
    )
    assert [g.label for g in groups] == ["neg", "pos"]


def test_column_group_labels_accept_key_mapping():
    groups = build_column_groups(
        columns=["s1", "s2", "s3", "s4"], key_frame=_status_frame(), keys=["Status"],
        order=["neg", "pos"], labels={"pos": "Positive"},
    )
    assert [g.label for g in groups] == ["neg", "Positive"]


def test_apply_column_group_labels_is_noop_without_labels():
    groups = (Group("a", "a", ("s1",)), Group("b", "b", ("s2",)))
    assert apply_column_group_labels(groups, None) == groups
    assert apply_column_group_labels(groups, []) == groups


def test_apply_column_group_labels_ignores_extra_labels():
    groups = (Group("a", "a", ("s1",)),)
    relabelled = apply_column_group_labels(groups, ["A", "B", "C"])
    assert [g.label for g in relabelled] == ["A"]


# --- gene groups -------------------------------------------------------------

def test_build_gene_groups_no_sets_returns_single_group():
    groups = build_gene_groups(genes=["A", "B"], gene_sets=None)
    assert len(groups) == 1
    assert groups[0].members == ("A", "B")


def test_build_gene_groups_partitions_with_other_trailing():
    groups = build_gene_groups(
        genes=["EGFR", "KRAS", "TP53", "MDM2", "XYZ"],
        gene_sets={"RTK/RAS": ["EGFR", "KRAS"], "TP53": ["TP53", "MDM2"]},
    )
    assert [g.key for g in groups[:-1]]  # first groups from sets
    other = groups[-1]
    assert other.key == "Other"
    assert other.members == ("XYZ",)


def test_build_gene_groups_members_follow_input_gene_order():
    groups = build_gene_groups(
        genes=["KRAS", "EGFR", "BRAF"],
        gene_sets={"RTK/RAS": ["EGFR", "KRAS", "BRAF"]},
        drop_ungrouped=True,
    )
    assert groups[0].members == ("KRAS", "EGFR", "BRAF")


def test_build_gene_groups_overlap_assigned_to_priority_group():
    groups = build_gene_groups(
        genes=["A", "B"],
        gene_sets={"first": ["A", "B"], "second": ["B"]},
        order=["first", "second"],
        drop_ungrouped=True,
    )
    by_key = {g.key: g for g in groups}
    assert by_key["first"].members == ("A", "B")
    assert "second" not in by_key  # B already claimed, second becomes empty and dropped


def test_build_gene_groups_drop_ungrouped():
    groups = build_gene_groups(
        genes=["A", "B", "C"],
        gene_sets={"grp": ["A"]},
        drop_ungrouped=True,
    )
    assert [g.key for g in groups] == ["grp"]


def test_build_gene_groups_other_always_last_even_with_order():
    groups = build_gene_groups(
        genes=["A", "B", "C"],
        gene_sets={"g1": ["A"], "g2": ["B"]},
        order=["g2", "g1"],
    )
    assert [g.key for g in groups] == ["g2", "g1", "Other"]


def test_build_gene_groups_follow_dict_insertion_order():
    # Row order follows the gene_sets insertion order, NOT group size.
    genes = ["A", "B", "C", "D", "E"]
    small_first = build_gene_groups(genes, {"small": ["A"], "big": ["B", "C", "D"]})
    assert [g.key for g in small_first] == ["small", "big", "Other"]
    # reversing the dict flips the row order
    big_first = build_gene_groups(genes, {"big": ["B", "C", "D"], "small": ["A"]})
    assert [g.key for g in big_first] == ["big", "small", "Other"]


def test_build_gene_groups_dict_order_unaffected_by_unrelated_order():
    # A group_order meant for column values (e.g. ["neg", "pos"]) must not disturb
    # the gene row order, which stays in gene_sets insertion order.
    genes = ["A", "B", "C"]
    groups = build_gene_groups(
        genes, {"first": ["A"], "second": ["B"]}, order=["neg", "pos"]
    )
    assert [g.key for g in groups] == ["first", "second", "Other"]


