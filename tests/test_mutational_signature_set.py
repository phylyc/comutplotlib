"""Unit tests for etiology grouping of mutational signatures."""

import sys

import pandas as pd
import pytest

from comutplotlib.comut_argparse import parse_args
from comutplotlib.comut_data import ComutData
from comutplotlib.mutational_signature_set import MutationalSignatureSet as MSS


def _exposures():
    return pd.DataFrame(
        {
            "SBS1": [1.0, 0.0],
            "SBS5": [2.0, 1.0],
            "SBS2": [3.0, 0.0],
            "SBS13": [4.0, 0.0],
            "SBS4": [0.0, 5.0],
        },
        index=pd.Index(["P1", "P2"], name="Patient"),
    )


def test_group_by_etiology_sums_members_of_a_set():
    grouped = MSS.group_by_etiology(_exposures())
    assert grouped["clock-like"].tolist() == [3.0, 1.0]
    assert grouped["APOBEC"].tolist() == [7.0, 0.0]
    assert grouped["Smoking"].tolist() == [0.0, 5.0]


def test_group_by_etiology_preserves_declaration_order():
    grouped = MSS.group_by_etiology(_exposures())
    assert list(grouped.columns) == ["clock-like", "APOBEC", "Smoking"]
    assert grouped.index.equals(_exposures().index)


def test_group_by_etiology_is_idempotent():
    grouped = MSS.group_by_etiology(_exposures())
    assert MSS.group_by_etiology(grouped).equals(grouped)


def test_group_by_etiology_keeps_unknown_signatures_as_own_column():
    df = _exposures().assign(SBS_made_up=[1.0, 2.0])
    grouped = MSS.group_by_etiology(df)
    assert grouped["SBS_made_up"].tolist() == [1.0, 2.0]
    # unknown categories are sorted to the end
    assert list(grouped.columns)[-1] == "SBS_made_up"


def test_group_by_etiology_passes_through_none():
    assert MSS.group_by_etiology(None) is None


def test_sort_signature_sets_appends_unknown_names():
    sets = pd.Index(["Smoking", "Mystery", "clock-like"])
    assert list(MSS.sort_signature_sets(sets)) == ["clock-like", "Smoking", "Mystery"]


@pytest.mark.parametrize("group", [True, False])
def test_comut_data_align_signatures_honours_the_flag(group):
    data = ComutData(group_signatures_by_etiology=group)
    data.signatures = _exposures()
    data.sif = None
    data.align_signatures()
    if group:
        assert list(data.signatures.columns) == ["clock-like", "APOBEC", "Smoking"]
    else:
        assert list(data.signatures.columns) == ["SBS1", "SBS5", "SBS2", "SBS13", "SBS4"]


def test_cli_flag_defaults_to_false_and_can_be_enabled(monkeypatch):
    base = ["ComutPlot", "-o", "out.png", "--maf", "in.maf"]
    monkeypatch.setattr(sys, "argv", base)
    assert parse_args().group_signatures_by_etiology is False
    monkeypatch.setattr(sys, "argv", base + ["--group-signatures-by-etiology"])
    assert parse_args().group_signatures_by_etiology is True

