"""Unit tests for CLI parsing/validation (``comutplotlib.comut_argparse``)."""

import sys

import pytest

from comutplotlib.comut_argparse import parse_args, validate_args
from comutplotlib.comut_panels import ComutPanels


def _parse(argv, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["ComutPlot"] + argv)
    return parse_args()


BASE = ["-o", "out.png", "--maf", "in.maf"]


def test_column_group_labels_parse_as_nargs(monkeypatch):
    args = _parse(
        BASE + [
            "--column-group-by", "Sample Type",
            "--column-group-labels", "Brain metastases", "Extracranial metastases",
        ],
        monkeypatch,
    )
    assert args.column_group_labels == ["Brain metastases", "Extracranial metastases"]


def test_column_group_labels_default_is_none(monkeypatch):
    args = _parse(BASE, monkeypatch)
    assert args.column_group_labels is None
    assert ComutPanels.column_group_label not in args.panels_to_plot


def test_hide_grouped_meta_data_flag(monkeypatch):
    args = _parse(
        BASE + ["--column-group-by", "Sample Type,Histology", "--hide-grouped-meta-data"],
        monkeypatch,
    )
    assert args.hide_grouped_meta_data is True
    assert args.column_group_by == ["Sample Type", "Histology"]
    assert {"Sample Type", "Histology"}.issubset(args.meta_data_rows)


def test_hide_grouped_meta_data_defaults_to_false(monkeypatch):
    args = _parse(BASE, monkeypatch)
    assert args.hide_grouped_meta_data is False


def test_column_group_labels_enable_the_label_panel(monkeypatch):
    args = _parse(
        BASE + ["--column-group-by", "Sample Type", "--column-group-labels", "A", "B"],
        monkeypatch,
    )
    assert ComutPanels.column_group_label in args.panels_to_plot
    # the stratification key is still auto-added to the metadata rows
    assert "Sample Type" in args.meta_data_rows


def test_column_group_labels_require_column_group_by(monkeypatch):
    with pytest.raises(ValueError, match="--column-group-labels requires --column-group-by"):
        _parse(BASE + ["--column-group-labels", "A"], monkeypatch)


def test_validate_args_does_not_duplicate_the_label_panel():
    class Args:
        maf = ["in.maf"]
        gistic = None
        mark = None
        control_maf = None
        control_gistic = None
        control_mark = None
        cohort_label = None
        control_cohort_label = None
        snv_interesting_genes = None
        cnv_interesting_genes = None
        signatures = None
        sif = None
        gene_meta_data = None
        meta_data_rows = ["Sample Type"]
        column_group_by = ["Sample Type"]
        column_group_labels = ["A"]
        panels_to_plot = [ComutPanels.column_group_label]

    args = Args()
    validate_args(args)
    assert args.panels_to_plot.count(ComutPanels.column_group_label) == 1
