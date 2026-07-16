from copy import deepcopy
import os
import pandas as pd
from tqdm import tqdm

from comutplotlib.config import load_config
from comutplotlib.mutation_annotation import MutationAnnotation as MutA


def join_gistics(gistics: list["Gistic"], verbose=False):
    if len(gistics) == 0:
        return Gistic()
    elif len(gistics) == 1:
        return gistics[0].copy()
    else:
        gistic = gistics[0]
        for other in tqdm(gistics[1:], desc="Joining GISTICs", total=len(gistics), initial=1, disable=not verbose):
            gistic = gistic.join(other=other)
        return gistic


class Gistic(object):
    """ GISTIC2 output table: all_thresholded.by_genes.txt
    """

    # Column names of Gistic output files
    _gene_symbol = "Gene Symbol"
    _gene_id = "Gene ID"
    _locus_id = "Locus ID"
    _cytoband = "Cytoband"

    # HGNC-based deprecated -> approved gene-symbol aliases, loaded from
    # comutplotlib/config/gene_aliases.json (override via $COMUTPLOTLIB_CONFIG_DIR).
    _gene_name_map = load_config("gene_aliases.json", default={}).get("aliases", {})

    @classmethod
    def from_file(cls, path_to_file: str, encoding: str = "utf8"):
        if not os.path.exists(path_to_file):
            raise FileNotFoundError(
                f"No Gistic instance created. No such file or directory: '{path_to_file}'"
            )
        try:
            with open(path_to_file, encoding=encoding) as f:
                data = pd.read_csv(filepath_or_buffer=f, sep="\t", engine="c", comment="#")
        except:
            data = pd.read_csv(filepath_or_buffer=path_to_file, sep="\t", engine="c", comment="#")
        data = data.set_index(cls._gene_symbol).rename_axis(MutA.gene_name, axis=0)
        return Gistic(data=data)

    def __init__(self, data=None):
        self.data = (
            data if data is not None
            else pd.DataFrame(None, index=pd.Index([], name=MutA.gene_name), columns=pd.Index([self._locus_id, self._cytoband]))
        )
        self._standardize_gene_names()

    def _standardize_gene_names(self):
        # self.data.index = self.data.index.map(lambda g: g.split("|")[0])
        for gene_name, replacement in self._gene_name_map.items():
            if gene_name in self.data.index:
                self.data.rename(index={gene_name: replacement}, inplace=True)

    @property
    def genes(self):
        return self.data.index

    @property
    def gene_id(self):
        return self.data[self._gene_id]

    @property
    def locus_id(self):
        return self.data[self._locus_id]

    @property
    def cytoband(self):
        return self.data[self._cytoband]

    @property
    def sample_table(self):
        drop_cols = [c for c in [self._gene_id, self._locus_id, self._cytoband] if c in self.data.columns]
        return self.data.drop(drop_cols, axis=1).rename_axis(MutA.sample, axis=1)

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

    def select_genes(self, genes: list[str], inplace=False) -> "Gistic":
        data = self.data.reindex(index=genes)
        if inplace:
            self.data = data
            return self
        return Gistic(data=data)

    def select_samples(self, samples: list[str], inplace=False) -> "Gistic":
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
