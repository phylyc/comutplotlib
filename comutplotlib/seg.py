from copy import deepcopy
import numpy as np
import os
import pandas as pd
from tqdm import tqdm
import pyensembl


def join_segs(segs: list["SEG"], verbose=False):
    if len(segs) == 0:
        return SEG()
    elif len(segs) == 1:
        return segs[0].copy()
    else:
        seg = segs[0]
        for other in tqdm(segs[1:], desc="Joining SEGs", total=len(segs), initial=1, disable=not verbose):
            seg = seg.join(other=other)
        return seg


class SEG(object):
    """ ABSOLUTE segtab output
    """

    # Column names of SEG output files
    _sample = "Sample"
    _contig = "Chromosome"
    _start = "Start Position"
    _end = "End Position"
    _num_markers = "Num Markers"
    _log2_tCR = "Seg.CN"
    _gene_name = "gene_name"

    @classmethod
    def from_file(cls, path_to_file: str, encoding: str = "utf8"):
        if not os.path.exists(path_to_file):
            raise FileNotFoundError(
                f"No SEG instance created. No such file or directory: '{path_to_file}'"
            )
        with open(path_to_file, encoding=encoding) as f:
            data = pd.read_csv(filepath_or_buffer=f, sep="\t", engine="c", comment="#")
        return SEG(data=data)

    def __init__(self, data=None):
        self.data = (
            data if data is not None
            else pd.DataFrame(None, columns=pd.Index([self._sample, self._contig, self._start, self._end, self._log2_tCR, self._gene_name]))
        )
        self.add_gene_name()

    def add_gene_name(self, verbose=False):
        data_is_incomplete = (
            self._gene_name not in self.data.columns
            or self.data[self._gene_name].isna().any()
        )
        if data_is_incomplete:
            genome = pyensembl.ensembl_grch37

            def get_gene_names(df):
                _gene_names = {}
                for index, segment in tqdm(
                    iterable=df.iterrows(),
                    total=df.shape[0],
                    desc="Annotate segs with genes",
                    disable=not verbose,
                ):
                    _gene_names[index] = genome.gene_names_at_locus(
                        contig=segment.get(self._contig),
                        position=segment.get(self._end),
                        end=segment.get(self._end),
                    )
                return pd.Series(_gene_names)

            self.data[self._gene_name] = get_gene_names(df=self.data)

    def get_df(self) -> pd.DataFrame:
        return (
            self.data[[self._sample, self._log2_tCR, self._gene_name]]
            .explode(self._gene_name)
            .groupby([self._sample, self._gene_name])
            .agg({self._log2_tCR: lambda x: x.values[np.argmax(np.abs(x.values))]})
            [self._log2_tCR]
            .unstack(self._sample)
        )

    def select_genes(self, genes: list[str], inplace=True) -> "SEG":
        gene_col = self.data[self._gene_name].apply(lambda _genes: [g in genes for g in _genes])
        mask = gene_col.apply(lambda l: len(l) > 0)
        data = self.data.loc[mask]
        data[self._gene_name] = gene_col
        if inplace:
            self.data = data
            return self
        else:
            return SEG(data=data)

    def select_samples(self, samples: list[str], inplace=True) -> "SEG":
        mask = self.data[self._sample].isin(samples)
        if inplace:
            self.data = self.data.loc[mask]
            return self
        return SEG(data=self.data.loc[mask])

    def copy(self) -> "SEG":
        seg = SEG()
        for name, value in self.__dict__.items():
            seg.__setattr__(name, deepcopy(value))
        return seg

    def join(self, other: "SEG") -> "SEG":
        joined = pd.concat([self.data, other.data], ignore_index=True)
        return SEG(data=joined)
