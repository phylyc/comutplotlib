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
            sif_meta_data["Sample Type"] = sif_meta_data[self.sif.sample_type]
        if self.sif.material in sif_meta_data.columns:
            sif_meta_data["Material"] = sif_meta_data[self.sif.material]
        if "contamination" in sif_meta_data.columns:
            sif_meta_data["Contamination"] = sif_meta_data["contamination"]
        if "tumor_purity" in sif_meta_data.columns:
            sif_meta_data["Tumor Purity"] = sif_meta_data["tumor_purity"]
        if self.sif.platform_abv in sif_meta_data.columns:
            sif_meta_data["Platform"] = sif_meta_data[self.sif.platform_abv]
        if "Paired" in sif_meta_data.columns:
            sif_meta_data["has matched N"] = sif_meta_data["Paired"].replace(True, "yes").replace(False, "no")
        if self.sif.sex in sif_meta_data.columns:
            sif_meta_data["Sex"] = sif_meta_data[self.sif.sex]
        if self.sif.histology in sif_meta_data.columns:
            sif_meta_data["Histology"] = sif_meta_data[self.sif.histology]

        self.rows = [c for c in self.rows if c in sif_meta_data.columns]
        self.rows_per_sample = [c for c in self.rows_per_sample if c in self.rows]
        per_patient_columns = [c for c in self.rows if c not in self.rows_per_sample]
        sif_meta_data = (
            sif_meta_data[[self.sif.sample, self.sif.patient] + list(self.rows)]
            .groupby(by=[self.by])
            .agg({k: list for k in self.rows_per_sample} | {k: lambda l: np.unique(l)[0] for k in per_patient_columns})
        )
        meta_data = sif_meta_data[self.rows] if self.rows is not None else pd.DataFrame()

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
        meta_data = meta_data[is_all_nan.loc[~is_all_nan].index]

        return meta_data

    def reindex(self, index=None, columns=None):
        if columns is not None:
            self.sif.select({columns.name: columns}, inplace=True)
            self.df = self.df.reindex(index=columns)

    @property
    def legend_titles(self):
        return self.df.columns

    def get_tmb(self):
        if self.sif.tmb in self.sif.data.columns:
            tmb = self.sif.data.groupby(by=[self.by]).agg({self.sif.tmb: np.mean})
            if self.sif.tmb_error in self.sif.data.columns:
                tmb_error = self.sif.data.groupby(by=[self.by]).agg({self.sif.tmb_error: np.mean})
                tmb = tmb.join(tmb_error)
            return tmb
        return None
