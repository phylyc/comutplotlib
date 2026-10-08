from copy import deepcopy
import os
import pandas as pd
from tqdm import tqdm

from comutplotlib.config import load_config
from comutplotlib.mutation_annotation import MutationAnnotation as MutA


def join_marks(marks: list["Mark"], verbose=False):
    if len(marks) == 0:
        return Mark()
    elif len(marks) == 1:
        return marks[0].copy()
    else:
        mark = marks[0]
        for other in tqdm(marks[1:], desc="Joining Marks", total=len(marks), initial=1, disable=not verbose):
            mark = mark.join(other=other)
        return mark


class Mark(object):
    """ Epigenetic mark table: genes (rows) by samples/patients (columns).

    The file format mirrors the GISTIC sample table (same gene-symbol index, same
    tab-separated layout), but without the 'Gene ID', 'Locus ID', and 'Cytoband'
    annotation columns. Values are floats between 0 and 1.
    """

    # Column names of Mark input files
    _gene_symbol = "Gene Symbol"

    # HGNC-based deprecated -> approved gene-symbol aliases, loaded from
    # comutplotlib/config/gene_aliases.json (override via $COMUTPLOTLIB_CONFIG_DIR).
    _gene_name_map = load_config("gene_aliases.json", default={}).get("aliases", {})

    @classmethod
    def from_file(cls, path_to_file: str, encoding: str = "utf8"):
        if not os.path.exists(path_to_file):
            raise FileNotFoundError(
                f"No Mark instance created. No such file or directory: '{path_to_file}'"
            )
        try:
            with open(path_to_file, encoding=encoding) as f:
                data = pd.read_csv(filepath_or_buffer=f, sep="\t", engine="c", comment="#")
        except (OSError, UnicodeDecodeError):
            data = pd.read_csv(filepath_or_buffer=path_to_file, sep="\t", engine="c", comment="#")
        data = data.set_index(data.columns[0]).rename_axis(MutA.gene_name, axis=0)
        return Mark(data=data)

    def __init__(self, data=None):
        self.data = (
            data if data is not None
            else pd.DataFrame(None, index=pd.Index([], name=MutA.gene_name), columns=pd.Index([]))
        )
        self._standardize_gene_names()

    def _standardize_gene_names(self):
        for gene_name, replacement in self._gene_name_map.items():
            if gene_name in self.data.index:
                self.data.rename(index={gene_name: replacement}, inplace=True)

    @property
    def genes(self):
        return self.data.index

    @property
    def sample_table(self):
        return self.data.rename_axis(MutA.sample, axis=1)

    @property
    def gene_names(self):
        return self.data.index

    @property
    def num_loci(self):
        return self.data.shape[0]

    @property
    def samples(self):
        return self.sample_table.columns

    @property
    def patients(self):
        return self.samples

    @property
    def num_samples(self):
        return self.sample_table.shape[1]

    def select_genes(self, genes: list[str], inplace=False) -> "Mark":
        data = self.data.reindex(index=genes)
        if inplace:
            self.data = data
            return self
        return Mark(data=data)

    def select_samples(self, samples: list[str], inplace=False) -> "Mark":
        columns = list(samples)
        if inplace:
            self.data = self.data.reindex(columns=columns)
            return self
        return Mark(data=self.data.reindex(columns=columns))

    def copy(self) -> "Mark":
        mark = Mark()
        for name, value in self.__dict__.items():
            mark.__setattr__(name, deepcopy(value))
        return mark

    def join(self, other: "Mark") -> "Mark":
        other_samples = [s for s in other.sample_table.columns if s not in self.sample_table.columns]
        joined = self.data.join(other.sample_table[other_samples], how="outer")
        return Mark(data=joined)

    def to_csv(self, path_to_file: str, **kwargs) -> None:
        self.data.rename_axis(self._gene_symbol, axis="index").to_csv(
            path_or_buf=path_to_file,
            header=True,
            index=True,
            sep="\t",
            na_rep="nan",
            mode="w+",
            **kwargs,
        )
