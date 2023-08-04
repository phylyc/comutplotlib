from copy import deepcopy
import os
import pandas as pd
from tqdm import tqdm

from comutplotlib.mutation_annotation import MutationAnnotation as MutA


def join_gistics(gistics: list["Gistic"]):
    if len(gistics) == 0:
        return Gistic()
    elif len(gistics) == 1:
        return gistics[0].copy()
    else:
        gistic = gistics[0]
        for other in tqdm(gistics[1:], desc="Joining GISTICs", total=len(gistics), initial=1):
            gistic = gistic.join(other=other)
        return gistic


class Gistic(object):

    # Column names of Gistic output files
    _gene_symbol = "Gene Symbol"
    _locus_id = "Locus ID"
    _cytoband = "Cytoband"

    @classmethod
    def from_file(cls, path_to_file: str, encoding: str = "utf8"):
        if not os.path.exists(path_to_file):
            raise FileNotFoundError(
                f"No Gistic instance created. No such file or directory: '{path_to_file}'"
            )
        with open(path_to_file, encoding=encoding) as f:
            data = (
                pd.read_csv(filepath_or_buffer=f, sep="\t", engine="c", comment="#")
                .set_index(cls._gene_symbol)
                .rename_axis(MutA.gene_name, axis=0)
            )
        return Gistic(data=data)

    def __init__(self, data=None):
        self.data = (
            data if data is not None
            else pd.DataFrame(None, index=pd.Index([], name=MutA.gene_name), columns=pd.Index([self._locus_id, self._cytoband]))
        )
        self._standardize_gene_names()

    def _standardize_gene_names(self):
        gene_name_map = {
            "MLL2": "KMT2D",
            "HIST1H3A": "H3C1",
            "HIST1H3B": "H3C2",
            "HIST1H3C": "H3C3",
            "HIST1H3D": "H3C4",
            "HIST1H3E": "H3C6",
            "HIST1H3F": "H3C7",
            "HIST1H3G": "H3C8",
            "HIST1H3H": "H3C10",
            "HIST1H3I": "H3C11",
            "HIST1H3J": "H3C12",
            "HIST2H3D": "H3C13",
            "HIST2H3C": "H3C14",
            "HIST2H3A": "H3C15",
        }

        for gene_name, replacement in gene_name_map.items():
            if gene_name in self.data.index:
                self.data.rename(index={gene_name: replacement}, inplace=True)

    @property
    def genes(self):
        return self.data.index

    @property
    def locus_id(self):
        return self.data[self._locus_id]

    @property
    def cytoband(self):
        return self.data[self._cytoband]

    @property
    def sample_table(self):
        return self.data.drop([self._locus_id, self._cytoband], axis=1).rename_axis(MutA.sample, axis=1)

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
    def num_samples(self):
        return self.sample_table.shape[1]

    def select_genes(self, genes: list[str], inplace=True) -> "Gistic":
        data = self.data.reindex(index=genes)
        if inplace:
            self.data = data
            return self
        return Gistic(data=data)

    def select_samples(self, samples: list[str], inplace=True) -> "Gistic":
        columns = [self._locus_id, self._cytoband] + list(samples)
        if inplace:
            self.data = self.data.reindex(columns=columns)
            return self
        return Gistic(data=self.data.reindex(columns=columns))

    def copy(self) -> "Gistic":
        gistic = Gistic()
        for name, value in self.__dict__.items():
            gistic.__setattr__(name, deepcopy(value))
        return gistic

    def join(self, other: "Gistic") -> "Gistic":
        other_samples = [s for s in other.sample_table.columns if s not in self.sample_table.columns]
        joined = self.data.join(other.sample_table[other_samples], how="outer")
        return Gistic(data=joined)
