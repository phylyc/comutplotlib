import numpy as np
import pandas as pd


class Meta(object):

    def __init__(self, sif, by, rows=(), rows_per_sample=()):
        self.sif = sif
        self.by = by
        self.rows = rows
        self.rows_per_sample = rows_per_sample
        self.df = self.get_df()

    def get_df(self):
        sif_meta_data = self.sif.data

        if self.sif.sample_type in sif_meta_data.columns:
            sif_meta_data["Sample Type"] = sif_meta_data[self.sif.sample_type].astype(str)
        if self.sif.material in sif_meta_data.columns:
            sif_meta_data["Material"] = sif_meta_data[self.sif.material].astype(str)
        if self.sif.contamination in sif_meta_data.columns:
            sif_meta_data["Contamination"] = sif_meta_data[self.sif.contamination].astype(float)
        if self.sif.tumor_purity in sif_meta_data.columns:
            sif_meta_data["Tumor Purity"] = sif_meta_data[self.sif.tumor_purity].astype(float)
        if self.sif.ploidy in sif_meta_data.columns:
            sif_meta_data["Ploidy"] = sif_meta_data[self.sif.ploidy].astype(float)
        if self.sif.genome_doublings in sif_meta_data.columns:
            sif_meta_data["WGD"] = sif_meta_data[self.sif.genome_doublings].fillna("nan")
        if self.sif.subclonal_genome_fraction in sif_meta_data.columns:
            sif_meta_data["Subclonal Fraction"] = sif_meta_data[self.sif.subclonal_genome_fraction].astype(float)
        if self.sif.platform_abv in sif_meta_data.columns:
            sif_meta_data["Platform"] = sif_meta_data[self.sif.platform_abv].astype(str)
        if self.sif.paired in sif_meta_data.columns:
            sif_meta_data["has matched N"] = sif_meta_data[self.sif.paired].replace(True, "yes").replace(False, "no").astype(str)
        if self.sif.sex in sif_meta_data.columns:
            sif_meta_data["Sex"] = sif_meta_data[self.sif.sex].astype(str)
        if self.sif.histology in sif_meta_data.columns:
            sif_meta_data["Histology"] = sif_meta_data[self.sif.histology].astype(str)

        self.rows = [c for c in self.rows if c in sif_meta_data.columns]
        self.rows_per_sample = [c for c in self.rows_per_sample if c in self.rows]

        if sif_meta_data.empty:
            return pd.DataFrame()

        per_patient_columns = [c for c in self.rows if c not in self.rows_per_sample]
        sif_meta_data = (
            sif_meta_data[[self.sif.sample, self.sif.patient] + list(self.rows)]
            .groupby(by=[self.by])
            .agg({k: list for k in self.rows_per_sample} | {k: lambda l: np.unique(l)[0] for k in per_patient_columns})
        )
        meta_data = sif_meta_data[self.rows] if self.rows is not None else pd.DataFrame()

        if meta_data.empty:
            return meta_data

        is_all_nan = meta_data.apply(
            lambda row: len([
                x
                for value in row
                for x in (value if isinstance(value, list) else [value])
                if not (isinstance(x, float) and np.isnan(x))
            ]),
            axis=0
        ) == 0
        is_all_nan |= meta_data.isna().all(axis=0)
        is_all_nan |= meta_data.isin(["unknown", "nan"]).all(axis=0)
        meta_data = meta_data[is_all_nan.loc[~is_all_nan].index]

        self.rows = [c for c in self.rows if c in meta_data.columns]
        self.rows_per_sample = [c for c in self.rows_per_sample if c in self.rows]

        return meta_data

    def reindex(self, index=None, columns=None):
        if index is not None:
            self.rows = index
            self.rows_per_sample = [r for r in index if r in self.rows_per_sample]
            self.df = self.df.reindex(columns=index)
        if columns is not None:
            self.sif.select({columns.name: columns}, inplace=True)
            self.df = self.df.reindex(index=columns)

    @property
    def legend_titles(self):
        return self.df.columns

    def get_tmb(self):
        if self.sif.tmb in self.sif.data.columns:
            tmb = self.sif.data.groupby(by=[self.by]).agg({self.sif.tmb: "mean"})
            if self.sif.tmb_error in self.sif.data.columns:
                tmb = tmb.join(self.sif.data.groupby(by=[self.by]).agg({self.sif.tmb_error: lambda e: np.sqrt(np.sum(e ** 2))}))
            if self.sif.n_vars in self.sif.data.columns:
                tmb = tmb.join(self.sif.data.groupby(by=[self.by]).agg({self.sif.n_vars: "sum"}))
            if self.sif.n_bases in self.sif.data.columns:
                tmb = tmb.join(self.sif.data.groupby(by=[self.by]).agg({self.sif.n_bases: "sum"}))
            return tmb
        return None
