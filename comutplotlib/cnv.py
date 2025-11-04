from collections import Counter

import pandas as pd


class CNV(object):

    def __init__(
        self, seg, gistic,
        low_amp_threshold: int | float = 1,
        mid_amp_threshold: int | float = 1.5,
        high_amp_threshold: int | float = 2,
        baseline: int | float = 0,
        low_del_threshold: int | float = -1,
        mid_del_threshold: int | float = -1.5,
        high_del_threshold: int | float = -2,
        higher_amp_thresholds: list[int | float] = tuple(),
    ):
        self.seg = seg
        self.gistic = gistic
        self.low_amp_threshold = low_amp_threshold
        self.mid_amp_threshold = mid_amp_threshold
        self.high_amp_threshold = high_amp_threshold
        self.baseline = baseline
        self.low_del_threshold = low_del_threshold
        self.mid_del_threshold = mid_del_threshold
        self.high_del_threshold = high_del_threshold
        self.higher_amp_thresholds = higher_amp_thresholds
        self.df = self.get_df()

    def get_df(self):
        if not self.gistic.data.empty:
            return self.gistic.sample_table
        elif not self.seg.data.empty:
            return self.seg.get_df()
        else:
            return pd.DataFrame(None)

    def reindex(self, index=None, columns=None):
        if index is not None:
            self.seg.select_genes(genes=index, inplace=True)
            self.gistic.select_genes(genes=index, inplace=True)
        if columns is not None:
            self.seg.select_samples(samples=columns, inplace=True)
            self.gistic.select_samples(samples=columns, inplace=True)
        self.df = self.get_df()

    def get_prevalence(self):
        return self.df.apply(Counter, axis=1)

    @property
    def isna(self):
        return self.df.isna()

    @property
    def has_high_amp(self):
        return ~self.isna & self.df.fillna(self.baseline).ge(self.high_amp_threshold)

    @property
    def has_high_del(self):
        return ~self.isna & self.df.fillna(self.baseline).le(self.high_del_threshold)

    @property
    def has_high_cnv(self):
        return self.has_high_amp | self.has_high_del

    @property
    def has_mid_amp(self):
        return ~self.isna & self.df.fillna(self.baseline).ge(self.mid_amp_threshold) & ~self.has_high_amp

    @property
    def has_mid_del(self):
        return ~self.isna & self.df.fillna(self.baseline).le(self.mid_del_threshold) & ~self.has_high_del

    @property
    def has_mid_cnv(self):
        return self.has_mid_amp | self.has_mid_del

    @property
    def has_low_amp(self):
        return ~self.isna & self.df.fillna(self.baseline).ge(self.low_amp_threshold) & ~self.has_high_amp & ~self.has_mid_amp

    @property
    def has_low_del(self):
        return ~self.isna & self.df.fillna(self.baseline).le(self.low_del_threshold) & ~self.has_high_del & ~self.has_mid_del

    @property
    def has_low_cnv(self):
        return self.has_low_amp | self.has_low_del

    @property
    def has_amp(self):
        return self.has_high_amp | self.has_mid_amp | self.has_low_amp

    @property
    def has_del(self):
        return self.has_high_del | self.has_mid_del | self.has_low_del

    @property
    def has_cnv(self):
        return self.has_amp | self.has_del
