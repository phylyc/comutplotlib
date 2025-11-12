import pandas as pd


class CNV(object):

    def __init__(
        self, seg, gistic,
        high_amp_threshold: int | float = 2,
        mid_amp_threshold: int | float = 1.5,
        low_amp_threshold: int | float = 1,
        baseline: int | float = 0,
        low_del_threshold: int | float = -1,
        mid_del_threshold: int | float = -1.5,
        high_del_threshold: int | float = -2,
        higher_amp_thresholds: list[int | float] = tuple(),
    ):
        self.seg = seg
        self.gistic = gistic
        self.high_amp_threshold = high_amp_threshold
        self.mid_amp_threshold = mid_amp_threshold
        self.low_amp_threshold = low_amp_threshold
        self.baseline = baseline
        self.low_del_threshold = low_del_threshold
        self.mid_del_threshold = mid_del_threshold
        self.high_del_threshold = high_del_threshold
        self.higher_amp_thresholds = higher_amp_thresholds
        self.amp_thresholds = sorted([high_amp_threshold, mid_amp_threshold, low_amp_threshold], reverse=True)
        self.del_thresholds = sorted([high_del_threshold, mid_del_threshold, low_del_threshold])
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

    @property
    def has_cn_by_level(self):
        return {
            self.high_amp_threshold: self.has_high_amp,
            self.mid_amp_threshold: self.has_mid_amp,
            self.low_amp_threshold: self.has_low_amp,
            self.baseline: self.is_baseline,
            self.low_del_threshold: self.has_low_del,
            self.mid_del_threshold: self.has_mid_del,
            self.high_del_threshold: self.has_high_del,
        }

    @property
    def has_cn_of_at_least_level(self):
        return {
            self.high_amp_threshold: self.has_high_amp,
            self.mid_amp_threshold: self.has_mid_amp,
            self.low_amp_threshold: self.has_high_amp | self.has_mid_amp | self.has_low_amp,
            self.baseline: self.is_baseline,
            self.low_del_threshold: self.has_high_del | self.has_mid_del | self.has_low_del,
            self.mid_del_threshold: self.has_mid_del,
            self.high_del_threshold: self.has_high_del,
        }

    def get_num_patients_by_gene_by_cn_level(self):
        return pd.DataFrame({k: v.sum(axis=1) for k, v in self.has_cn_by_level.items()})

    def get_num_patients_by_gene_of_at_least_cn_level(self):
        return pd.DataFrame({k: v.sum(axis=1) for k, v in self.has_cn_of_at_least_level.items()})

    @property
    def isna(self):
        return self.df.isna()

    def eq(self, threshold):
        return ~self.isna & self.df.fillna(self.baseline).eq(threshold)

    def ge(self, threshold):
        return ~self.isna & self.df.fillna(self.baseline).ge(threshold)

    def le(self, threshold):
        return ~self.isna & self.df.fillna(self.baseline).le(threshold)

    @property
    def is_baseline(self):
        return self.eq(self.baseline)

    @property
    def has_high_amp(self):
        return self.ge(self.high_amp_threshold)

    @property
    def has_high_del(self):
        return self.le(self.high_del_threshold)

    @property
    def has_high_cnv(self):
        return self.has_high_amp | self.has_high_del

    @property
    def has_mid_amp(self):
        return self.ge(self.mid_amp_threshold) & ~self.has_high_amp

    @property
    def has_mid_del(self):
        return self.le(self.mid_del_threshold) & ~self.has_high_del

    @property
    def has_mid_cnv(self):
        return self.has_mid_amp | self.has_mid_del

    @property
    def has_low_amp(self):
        return self.ge(self.low_amp_threshold) & ~self.has_high_amp & ~self.has_mid_amp

    @property
    def has_low_del(self):
        return self.le(self.low_del_threshold) & ~self.has_high_del & ~self.has_mid_del

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
