"""Integration tests for the grid-aware ``ComutData`` behaviour.

These use the committed demo inputs (``demo/*.tsv``) to exercise the real
SNV/CNV/meta pipeline. They deliberately avoid the palette/plotter layers so they
stay fast and independent of rendering.
"""

import copy
import os

import pytest

from comutplotlib.comut_data import ComutData
from comutplotlib.mutation_annotation import MutationAnnotation as MutA

DEMO = os.path.join(os.path.dirname(__file__), "..", "demo")


def _has_demo() -> bool:
    return os.path.exists(os.path.join(DEMO, "test.maf.tsv"))


pytestmark = pytest.mark.skipif(not _has_demo(), reason="demo inputs not available")


@pytest.fixture(scope="module")
def base_data():
    data = ComutData(
        maf_paths=[os.path.join(DEMO, "test.maf.tsv")],
        gistic_paths=[os.path.join(DEMO, "test.all_thresholded.by_genes.txt")],
        sif_paths=[os.path.join(DEMO, "test.sif.tsv")],
        meta_data_rows=["Sample Type", "Histology", "Platform", "Material"],
        by=MutA.patient,
        snv_interesting_genes={f"Gene {i}" for i in range(1, 13)},
        cnv_interesting_genes={f"Gene {i}" for i in range(9, 21)},
        column_sort_by=("COMUT", "TMB"),
    )
    data.preprocess()
    return data


@pytest.fixture
def data(base_data):
    # deepcopy so each test can mutate/reindex independently and cheaply
    return copy.deepcopy(base_data)


# --- refactor equivalence ----------------------------------------------------

def test_sort_refactor_is_idempotent(data):
    """The parameterized sort over the full slice reproduces the stored order."""
    genes = list(data.genes)
    columns = list(data.columns)
    assert list(data._compute_gene_order(data.genes, data.columns)) == genes
    assert list(data._compute_column_order(data.columns, data.genes)) == columns


def test_compute_gene_order_returns_permutation_for_subset(data):
    subset = list(data.genes[:8])
    order = list(data._compute_gene_order(subset, data.columns))
    assert sorted(order) == sorted(subset)


def test_compute_column_order_returns_permutation_for_subset(data):
    subset = list(data.columns[:30])
    order = list(data._compute_column_order(subset, data.genes))
    assert sorted(order) == sorted(subset)


# --- grid construction -------------------------------------------------------

def test_trivial_grid_when_no_grouping(data):
    genes_before = set(data.genes)
    columns_before = set(data.columns)
    gene_groups = data.apply_grid_sort()
    assert data.grid.is_trivial()
    assert len(gene_groups) == 1
    assert set(data.genes) == genes_before
    assert set(data.columns) == columns_before


def test_column_groups_from_metadata(data):
    data.column_group_by = ["Sample Type"]
    data.apply_grid_sort()
    keys = {g.key for g in data.grid.column_groups}
    assert keys == {"BM", "EM"}
    # partition covers all columns exactly once
    members = [m for g in data.grid.column_groups for m in g.members]
    assert sorted(members) == sorted(data.columns)
    assert len(members) == len(set(members))


def test_column_group_labels_override_default_labels(data):
    data.column_group_by = ["Sample Type"]
    data.column_group_labels = ["Brain metastases", "Extracranial metastases"]
    data.apply_grid_sort()
    # keys stay the raw metadata values; only the display labels change
    assert [g.key for g in data.grid.column_groups] != [g.label for g in data.grid.column_groups]
    assert [g.label for g in data.grid.column_groups] == [
        "Brain metastases", "Extracranial metastases"
    ]
    assert data.column_group_label_map() == {
        g.key: g.label for g in data.grid.column_groups
    }


def test_column_group_labels_as_key_mapping(data):
    data.column_group_by = ["Sample Type"]
    data.column_group_labels = {"BM": "Brain metastases"}
    data.apply_grid_sort()
    labels = {g.key: g.label for g in data.grid.column_groups}
    assert labels["BM"] == "Brain metastases"
    assert labels["EM"] == "EM"  # untouched


def test_gene_groups_with_other(data):
    data.gene_groups_config = {"setA": ["Gene 1", "Gene 2", "Gene 3"], "setB": ["Gene 10", "Gene 11"]}
    data.apply_grid_sort()
    keys = [g.key for g in data.grid.gene_groups]
    assert keys[-1] == "Other"
    assert set(keys) >= {"setA", "setB", "Other"}
    members = [m for g in data.grid.gene_groups for m in g.members]
    assert sorted(members) == sorted(data.genes)


def test_apply_grid_sort_reassembles_genes_and_columns(data):
    data.column_group_by = ["Sample Type"]
    data.gene_groups_config = {"setA": ["Gene 1", "Gene 2", "Gene 3"]}
    data.apply_grid_sort()
    assert list(data.genes) == list(data.grid.all_genes())
    assert list(data.columns) == list(data.grid.all_columns())
    # data carriers are aligned to the reassembled axes
    assert list(data.snv.df.columns) == list(data.columns)
    assert list(data.cnv.df.columns) == list(data.columns)


def test_within_group_columns_are_sorted_consistently(data):
    """Each column group is internally ordered; a group's members stay contiguous."""
    data.column_group_by = ["Sample Type"]
    data.apply_grid_sort()
    all_cols = list(data.columns)
    for group in data.grid.column_groups:
        idx = [all_cols.index(m) for m in group.members]
        assert idx == list(range(min(idx), max(idx) + 1)), "group members must be contiguous"


def test_control_adopts_case_gene_groups(base_data):
    case = copy.deepcopy(base_data)
    control = copy.deepcopy(base_data)
    case.gene_groups_config = {"setA": ["Gene 1", "Gene 2", "Gene 3"], "setB": ["Gene 10", "Gene 11"]}
    case.column_group_by = ["Sample Type"]
    case_gene_groups = case.apply_grid_sort()

    control.column_group_by = ["Sample Type"]
    control.apply_grid_sort(gene_groups=case_gene_groups)

    # gene axis (rows) identical between cohorts
    assert list(control.genes) == list(case.genes)
    assert [g.key for g in control.grid.gene_groups] == [g.key for g in case.grid.gene_groups]
    for cg, kg in zip(control.grid.gene_groups, case.grid.gene_groups):
        assert list(cg.members) == list(kg.members)

