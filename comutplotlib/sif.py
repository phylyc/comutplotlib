from copy import deepcopy
import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import Callable, Optional, Union

from comutplotlib.annotation_table import AnnotationTable
from comutplotlib.sample_annotation import SampleAnnotation
from comutplotlib.sample_classification import classify
from comutplotlib import pandas_util as pd_util


def join_sifs(sifs: list["SIF"], verbose=False):
    if len(sifs) == 0:
        return SIF()
    elif len(sifs) == 1:
        return sifs[0].copy()
    else:
        sif = sifs[0]
        for other in tqdm(sifs[1:], desc="Joining SIFs", total=len(sifs), initial=1, disable=not verbose):
            sif = sif.join(other=other)
        return sif


class SIF(SampleAnnotation, AnnotationTable):

    @classmethod
    def from_file(
        cls,
        path_to_file: str,
        selection: Union[Callable[..., bool], dict, None] = None,
        complement: bool = False,
        encoding: str = "utf8",
        **kwargs
    ):
        # TODO: refactor redundant double selection
        annot = super().from_file(
            path_to_file=path_to_file,
            selection=selection,
            complement=complement,
            encoding=encoding,
            **kwargs
        )
        return SIF(data=annot.data, selection=annot.selection, complement=complement)

    def __init__(self, data=None, selection=None, complement=False):
        super().__init__(data=data)
        if data is None:
            self.data = pd.DataFrame(data=None, columns=self.default_columns)
        self.enforce_dtype()
        self.select(selection=selection, complement=complement, inplace=True)

    def enforce_dtype(self):
        for key, dtype in self.column_dtype.items():
            if key in self.data.columns:
                with pd_util.PandasChainedAssignmentWarnHandler():
                    self.data[key] = self.data[key].astype(dtype)

    def join(self, other: "SIF"):
        joined = pd.concat([self.data, other.data], ignore_index=True)
        return SIF(data=joined)

    def copy(self) -> "SIF":
        sif = SIF()
        for name, value in self.__dict__.items():
            sif.__setattr__(name, deepcopy(value))
        return sif

    def assign_column(self, name: str, value, inplace: bool = False):
        if inplace:
            self.update_inplace(data=self.data.assign(**{name: value}))
        else:
            _self = self.copy()
            _self.update_inplace(data=_self.data.assign(**{name: value}))
            return _self

    def pool_annotations(
        self, pool_as: dict, regex: bool = False, inplace: bool = False
    ):
        if pool_as is None:
            return None if inplace else self
        if inplace:
            for pooled_value, to_replace in pool_as.items():
                self.data.replace(
                    to_replace=to_replace,
                    value=pooled_value,
                    regex=regex,
                    inplace=True,
                )
        else:
            data = self.data.copy()
            for pooled_value, to_replace in pool_as.items():
                data = data.replace(
                    to_replace=to_replace,
                    value=pooled_value,
                    regex=regex,
                )
            pooled_sif = self.copy()
            pooled_sif.data = data
            return pooled_sif

    def select(
        self,
        selection: Union[Callable[..., bool], dict, None],
        complement: bool = False,
        inplace: bool = False,
    ) -> Optional["SIF"]:
        if inplace:
            super().select(selection=selection, complement=complement, inplace=inplace)
        else:
            if selection is None:
                return self
            else:
                return SIF(data=self.data, selection=selection, complement=complement)

    def get_entry(
        self,
        sample: str | None = None,
        platform: str | None = None,
        data_type: str | None = None,
        histology: str | None = None,
        sample_type: str | None = None,
    ):
        mask = True
        if sample is not None:
            mask &= self.data[self.sample] == sample
        if data_type is not None:
            mask &= self.data[self.data_type] == data_type
        if platform is not None and self.platform_abv in self.data.columns:
            mask &= self.data[self.platform_abv] == platform
        if histology is not None and self.histology in self.data.columns:
            mask &= self.data[self.histology] == histology
        if sample_type is not None:
            mask &= self.data[self.sample_type] == sample_type
        return self.data.loc[mask]

    @property
    def empty(self) -> bool:
        return self.data.empty

    @property
    def num_patients(self) -> int:
        return self.patients.shape[0]

    @property
    def num_samples(self) -> int:
        return self.samples.shape[0]

    @property
    def patients(self) -> np.ndarray:
        patients = self.data.get(self.patient)
        if patients is not None:
            return patients.unique()
        else:
            return np.array([])

    @property
    def samples(self) -> np.ndarray:
        samples = self.data.get(self.sample)
        if samples is not None:
            return samples.unique()
        else:
            return np.array([])

    def get_entries(
        self,
        sample: str | None = None,
        patient: str | None = None,
        platform: str | None = None,
        data_type: str | None = None,
        histology: str | None = None,
        sample_type: str | None = None,
    ) -> pd.DataFrame:
        mask = True
        if sample is not None:
            mask &= self.data[self.sample] == sample
        if patient is not None:
            mask &= self.data[self.patient] == patient
        if data_type is not None:
            mask &= self.data[self.data_type] == data_type
        if platform is not None and self.platform_abv in self.data.columns:
            mask &= self.data[self.platform_abv] == platform
        if histology is not None and self.histology in self.data.columns:
            mask &= self.data[self.histology] == histology
        if sample_type is not None:
            mask &= self.data[self.sample_type] == sample_type
        return self.data.loc[mask]

    def get_patient(self, sample: str):
        patient_series = self.get_entries(sample=sample)[self.patient]
        if not patient_series.empty:
            return patient_series.to_numpy()[0]
        else:
            return None

    def get_samples(self, patient: str):
        return self.get_entries(patient=patient)[self.sample].unique()

    def get_matched_normal_sample(self, sample: str):
        patient = self.get_patient(sample=sample)
        if patient is not None:
            normal_samples = self.get_entries(patient=patient, sample_type="N")[
                self.sample
            ].to_numpy()
            if len(normal_samples) and sample not in normal_samples:
                return normal_samples[0]
        return None

    def add_annotations(self, force: bool = False, inplace: bool = False):
        sif = self.copy()
        sif.data = sif.data.loc[~sif.data[self.sample].isna()]
        cols_to_drop = sif.data.apply(
            lambda col: (
                len(col.unique()) == 1
                and (
                    np.isnan(col.unique()[0])
                    if isinstance(col.unique()[0], float)
                    else False
                )
            )
        )
        sif.data = sif.data.drop(cols_to_drop.loc[cols_to_drop].index, axis=1)

        # use sample ID for patients without clinical ID
        na_patients = sif.data[self.patient].isna()
        sif.data.loc[na_patients, self.patient] = sif.data.loc[na_patients, self.sample]

        # remove leading and trailing whitespace
        for column_name, column in sif.data.items():
            if column.dtype == object:
                sif.data[column_name] = column.map(
                    lambda v: v.strip() if isinstance(v, str) else v
                )

        if self.cancer_type in sif.data.columns:
            sif.data.loc[sif.data[self.cancer_type].isna(), self.cancer_type] = "NA"

        if self.histotype in sif.data.columns:
            sif.data.loc[sif.data[self.histotype].isna(), self.histotype] = "NA"

        if self.cancer_type in sif.data.columns:
            if self.histology in sif.data.columns and not force:
                pass
            else:
                sif.data[self.histology] = sif.data.apply(
                    lambda s: classify("histology", {
                        "cancer_type": s.get(self.cancer_type, ""),
                        "histotype": s.get(self.histotype, ""),
                    }),
                    axis=1,
                )

        if self.center in sif.data.columns:
            sif.data.loc[sif.data[self.center].isna(), self.center] = "NA"

        if self.platform in sif.data.columns:
            sif.data.loc[sif.data[self.platform].isna(), self.platform] = "NA"

        if self.platform in sif.data.columns or self.center in sif.data.columns:
            if self.platform_abv in sif.data.columns and not force:
                pass
            else:
                sif.data[self.platform_abv] = sif.data.apply(
                    lambda s: classify("platform", {
                        "platform": s.get(self.platform, ""),
                        "center": s.get(self.center, ""),
                    }),
                    axis=1,
                )

        if self.data_type in sif.data.columns:
            sif.data.loc[sif.data[self.data_type].isna(), self.data_type] = "NA"

        if self.sample_type_long in sif.data.columns:
            if self.sample_type in sif.data.columns and not force:
                pass
            else:
                sif.data[self.sample_type] = sif.data.apply(
                    lambda s: classify("sample_type", {
                        "sample_type": s.get(self.sample_type_long, ""),
                        "sample_description": s.get(self.sample_description, ""),
                    }),
                    axis=1,
                )

        if self.material in sif.data.columns:
            sif.data[self.material] = sif.data.apply(
                lambda s: classify("sample_material", {
                    "sample_material": s.get(self.material, ""),
                }),
                axis=1,
            )

        if self.sex in sif.data.columns:
            sif.data[SIF.sex_genotype] = (
                sif.data[self.sex]
                .str.replace("FEMALE", "XX")
                .str.replace("Female", "XX")
                .str.replace("female", "XX")
                .str.replace("MALE", "XY")
                .str.replace("Male", "XY")
                .str.replace("male", "XY")
                .fillna("unknown")
            )

        if self.sample_type in sif.data.columns:
            sif.data[SIF.is_paired] = sif.data[SIF.sample].apply(
                lambda s: sif.get_matched_normal_sample(sample=s) is not None
            )

        if self.sample_type in sif.data.columns:
            sif.data[SIF.has_metastasis] = sif.data[SIF.sample].apply(
                lambda s: sif.get_entries(patient=sif.get_patient(sample=s))[SIF.sample_type].isin(["BM", "EM"]).any()
            )

        if inplace:
            self.data = sif.data
        else:
            return sif
