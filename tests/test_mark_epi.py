"""Tests for the epigenetic mark input (``Mark``) and its ``EPI`` wrapper.

These exercise the file parsing, joining, and reindexing contract, plus the
``ComutData`` integration that aligns marks with the comutation gene/column axes.
"""

import os

import numpy as np
import pandas as pd
import pytest

from comutplotlib.comut_data import ComutData
from comutplotlib.epi import EPI
from comutplotlib.mark import Mark, join_marks
from comutplotlib.mutation_annotation import MutationAnnotation as MutA

DEMO = os.path.join(os.path.dirname(__file__), "..", "demo")


def _has_demo() -> bool:
    return os.path.exists(os.path.join(DEMO, "test.marks.by_genes.txt"))


def _mark(values, genes, samples):
    data = pd.DataFrame(values, index=pd.Index(genes, name=MutA.gene_name), columns=samples)
    return Mark(data=data)


def test_from_file_round_trip(tmp_path):
    mark = _mark([[0.0, 0.5], [1.0, np.nan]], ["Gene 1", "Gene 2"], ["S1", "S2"])
    path = str(tmp_path / "marks.txt")
    mark.to_csv(path)

    loaded = Mark.from_file(path_to_file=path)
    assert loaded.genes.name == MutA.gene_name
    assert list(loaded.genes) == ["Gene 1", "Gene 2"]
    assert list(loaded.samples) == ["S1", "S2"]
    assert list(loaded.patients) == ["S1", "S2"]
    assert loaded.num_loci == 2
    assert loaded.num_samples == 2
    pd.testing.assert_frame_equal(loaded.data, mark.data, check_names=False)


def test_from_file_missing_path_raises():
    with pytest.raises(FileNotFoundError):
        Mark.from_file(path_to_file="/does/not/exist.marks.txt")


def test_empty_mark_has_no_samples():
    mark = Mark()
    assert mark.data.empty
    assert mark.num_samples == 0
    assert EPI(mark=mark).empty


def test_sample_table_names_the_column_axis():
    mark = _mark([[0.5]], ["Gene 1"], ["S1"])
    assert mark.sample_table.columns.name == MutA.sample


def test_join_marks_unions_genes_and_keeps_distinct_samples():
    a = _mark([[0.1], [0.2]], ["Gene 1", "Gene 2"], ["S1"])
    b = _mark([[0.3], [0.4]], ["Gene 2", "Gene 3"], ["S2"])

    joined = join_marks([a, b])

    assert sorted(joined.genes) == ["Gene 1", "Gene 2", "Gene 3"]
    assert list(joined.samples) == ["S1", "S2"]
    assert joined.data.loc["Gene 2", "S2"] == 0.3
    assert np.isnan(joined.data.loc["Gene 1", "S2"])
    # joining must not mutate the inputs
    assert list(a.samples) == ["S1"]


def test_join_marks_with_zero_and_one_element():
    assert join_marks([]).data.empty
    mark = _mark([[0.5]], ["Gene 1"], ["S1"])
    assert join_marks([mark]).data.equals(mark.data)


def test_epi_has_mark_ignores_nan_and_zero():
    epi = EPI(mark=_mark([[0.0, 0.5], [np.nan, 1.0]], ["Gene 1", "Gene 2"], ["S1", "S2"]))
    expected = pd.DataFrame(
        [[False, True], [False, True]],
        index=pd.Index(["Gene 1", "Gene 2"], name=MutA.gene_name),
        columns=pd.Index(["S1", "S2"], name=MutA.sample),
    )
    pd.testing.assert_frame_equal(epi.has_mark, expected)
    pd.testing.assert_series_equal(
        epi.get_num_patients_by_gene(),
        pd.Series([1, 1], index=expected.index),
    )


def test_epi_reindex_aligns_both_axes():
    epi = EPI(mark=_mark([[0.5, 0.25]], ["Gene 1"], ["S1", "S2"]))
    epi.reindex(index=["Gene 1", "Gene 9"], columns=["S2", "S3"])

    assert list(epi.df.index) == ["Gene 1", "Gene 9"]
    assert list(epi.df.columns) == ["S2", "S3"]
    assert epi.df.loc["Gene 1", "S2"] == 0.25
    assert np.isnan(epi.df.loc["Gene 9", "S3"])


def test_epi_legend_alphas_are_graded_or_binary():
    binary = EPI(mark=_mark([[0.0, 1.0]], ["Gene 1"], ["S1", "S2"]))
    assert binary.is_binary
    assert binary.legend_alphas == EPI.legend_alphas_binary

    graded = EPI(mark=_mark([[0.0, 0.4]], ["Gene 1"], ["S1", "S2"]))
    assert not graded.is_binary
    assert graded.legend_alphas == EPI.legend_alphas_graded


def test_epi_is_binary_ignores_nan():
    epi = EPI(mark=_mark([[np.nan, 1.0]], ["Gene 1"], ["S1", "S2"]))
    assert epi.is_binary


@pytest.mark.skipif(not _has_demo(), reason="demo inputs not available")
def test_comut_data_aligns_marks_with_the_comut_axes():
    data = ComutData(
        maf_paths=[os.path.join(DEMO, "test.maf.tsv")],
        gistic_paths=[os.path.join(DEMO, "test.all_thresholded.by_genes.txt")],
        mark_paths=[os.path.join(DEMO, "test.marks.by_genes.txt")],
        sif_paths=[os.path.join(DEMO, "test.sif.tsv")],
        by=MutA.patient,
        snv_interesting_genes={f"Gene {i}" for i in range(1, 13)},
        cnv_interesting_genes={f"Gene {i}" for i in range(9, 21)},
    )
    data.preprocess()

    assert list(data.epi.df.index) == list(data.cnv.df.index)
    assert list(data.epi.df.columns) == list(data.cnv.df.columns)
    assert not data.epi.empty


@pytest.mark.skipif(not _has_demo(), reason="demo inputs not available")
def test_comut_data_without_marks_yields_an_empty_epi():
    data = ComutData(
        maf_paths=[os.path.join(DEMO, "test.maf.tsv")],
        gistic_paths=[os.path.join(DEMO, "test.all_thresholded.by_genes.txt")],
        sif_paths=[os.path.join(DEMO, "test.sif.tsv")],
        by=MutA.patient,
        snv_interesting_genes={f"Gene {i}" for i in range(1, 13)},
    )
    data.preprocess()

    assert data.epi.empty
    assert data.epi.legend_alphas == EPI.legend_alphas_binary
