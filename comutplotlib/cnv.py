from collections import Counter


class CNV(object):

    def __init__(
        self, seg, gistic,
        low_amp_threshold: int | float = 1,
        mid_amp_threshold: int | float = 1.5,
        high_amp_threshold: int | float = 2,
        low_del_threshold: int | float = -1,
        mid_del_threshold: int | float = -1.5,
        high_del_threshold: int | float = -2,
    ):
        self.seg = seg
        self.gistic = gistic
        self.low_amp_threshold = low_amp_threshold
        self.mid_amp_threshold = mid_amp_threshold
        self.high_amp_threshold = high_amp_threshold
        self.low_del_threshold = low_del_threshold
        self.mid_del_threshold = mid_del_threshold
        self.high_del_threshold = high_del_threshold
        self.df = self.get_df()

    def get_df(self):
        return self.gistic.sample_table

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
        return ~self.isna & self.df.fillna(0).ge(self.high_amp_threshold)

    @property
    def has_high_del(self):
        return ~self.isna & self.df.fillna(0).le(self.high_del_threshold)

    @property
    def has_high_cnv(self):
        return self.has_high_amp | self.has_high_del

    @property
    def has_mid_amp(self):
        return ~self.isna & self.df.fillna(0).ge(self.mid_amp_threshold) & ~self.has_high_amp

    @property
    def has_mid_del(self):
        return ~self.isna & self.df.fillna(0).le(self.mid_del_threshold) & ~self.has_high_del

    @property
    def has_mid_cnv(self):
        return self.has_mid_amp | self.has_mid_del

    @property
    def has_low_amp(self):
        return ~self.isna & self.df.fillna(0).ge(self.low_amp_threshold) & ~self.has_high_amp & ~self.has_mid_amp

    @property
    def has_low_del(self):
        return ~self.isna & self.df.fillna(0).le(self.low_del_threshold) & ~self.has_high_del & ~self.has_mid_del

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
