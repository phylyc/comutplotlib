import numpy as np
import pandas as pd

from comutplotlib.mark import Mark


class EPI(object):

    legend_alphas_graded = (0.25, 0.5, 0.75, 1.0)
    legend_alphas_binary = (1.0,)

    def __init__(
        self,
        mark = None,
        baseline: int | float = 0,
    ):
        self.mark = mark if mark is not None else Mark()
        self.baseline = baseline
        self.df = self.get_df()

    def get_df(self):
        if not self.mark.data.empty:
            return self.mark.sample_table
        else:
            return pd.DataFrame(None)

    def reindex(self, index=None, columns=None):
        if index is not None:
            self.mark.select_genes(genes=index, inplace=True)
        if columns is not None:
            self.mark.select_samples(samples=columns, inplace=True)
        self.df = self.get_df()

    @property
    def empty(self):
        return self.df.empty or self.isna.all().all()

    @property
    def isna(self):
        return self.df.isna()

    @property
    def values(self):
        return self.df.to_numpy(dtype=float).flatten() if not self.df.empty else np.array([])

    @property
    def has_mark(self) -> pd.DataFrame:
        return ~self.isna & self.df.fillna(self.baseline).gt(self.baseline)

    @property
    def is_binary(self) -> bool:
        values = self.values
        values = values[~np.isnan(values)]
        return bool(np.isin(values, [0, 1]).all())

    @property
    def legend_alphas(self) -> tuple[float, ...]:
        return self.legend_alphas_binary if self.is_binary else self.legend_alphas_graded

    def get_num_patients_by_gene(self):
        return self.has_mark.sum(axis=1)
