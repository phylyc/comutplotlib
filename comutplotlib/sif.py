from copy import deepcopy
import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import Callable, Optional, Union

from comutplotlib.annotation_table import AnnotationTable
from comutplotlib.sample_annotation import SampleAnnotation


def join_sifs(sifs: list["SIF"]):
    if len(sifs) == 0:
        return SIF()
    elif len(sifs) == 1:
        return sifs[0].copy()
    else:
        sif = sifs[0]
        for other in tqdm(sifs[1:], desc="Joining SIFs", total=len(sifs), initial=1):
            sif = sif.join(other=other)
        return sif


class SIF(SampleAnnotation, AnnotationTable):

    @classmethod
    def from_file(
        cls,
        path_to_file: str,
        selection: Union[Callable[..., bool], dict] = None,
        complement: bool = False,
        encoding: str = "utf8",
    ):
        # TODO: refactor redundant double selection
        annot = super().from_file(
            path_to_file=path_to_file,
            selection=selection,
            complement=complement,
            encoding=encoding,
        )
        return SIF(data=annot.data, selection=annot.selection, complement=complement)

    def __init__(self, data=None, selection=None, complement=False):
        super().__init__(data=data)
        if data is None:
            self.data = pd.DataFrame(data=None, columns=self.default_columns)
        self.select(selection=selection, complement=complement, inplace=True)

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
        selection: Union[Callable[..., bool], dict],
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
        sample: str = None,
        platform: str = None,
        data_type: str = None,
        histology: str = None,
        sample_type: str = None,
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

    def add_annotations(self, inplace: bool = False):
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
            if isinstance(sif.data[column_name].dtype, str):
                sif.data[column_name] = column.str.strip()

        def get_histology(cancer_type, histotype):
            if "Breast" in cancer_type:
                return "BRCA"

            elif "Lung adenocarcinoma" in histotype or "Lung Adenocarcinoma" in histotype:
                return "LUAD"
            # elif "Adenosquamous lung carcinoma" in histotype:
            #     return "ALCA"
            # elif "Pleomorphic lung carcinoma" in histotype:
            #     return "PLCA"
            # elif "Squamous cell lung carcinoma" in histotype:
            #     return "SCLC"
            # elif "Non-small cell lung carcinoma" in histotype:
            #     return "NSCLC"
            # elif "Large cell lung carcinoma" in histotype:
            #     return "LCLC"
            elif "Lung cancer" in cancer_type:
                return "LUCA"

            elif "Melanoma" in cancer_type:
                return "MEL"

            elif "Renal cell carcinoma" in cancer_type or "Renal Clear Cell Carcinoma" in cancer_type:
                return "RCC"

            elif "Glioma" in cancer_type:
                return "GLIO"
            elif "Ependymoma" in cancer_type:
                return "GLIO"

            elif "Central Neurocytoma" in cancer_type:
                return "CNC"

            elif "Meningioma" in cancer_type:
                return "MEN"

            elif "Colon cancer" in cancer_type:
                return "GI"
            elif "Rectal cancer" in cancer_type:
                return "GI"
            elif "Colorectal cancer" in cancer_type:
                return "GI"
            elif "Duodenal cancer" in cancer_type:
                return "GI"
            elif "Esophageal cancer" in cancer_type:
                return "GI"
            elif "Pancreatic cancer" in cancer_type:
                return "GI"
            elif "GI cancer" in cancer_type:
                return "GI"
            elif "Gastric cancer" in cancer_type:
                return "GI"

            elif "Head & neck cancer" in cancer_type:
                return "HNC"
            elif "Oral cancer" in cancer_type:
                return "HNC"  # "ORCA"

            elif "Laryngeal cancer" in cancer_type:
                return "LACA"

            elif "Leukemia" in cancer_type:
                return "LEUK"

            elif "PCNSL" in histotype:  # Lymphoma
                return "PCNSL"
            # elif "Lymphoma" in cancer_type:
            #     return "LYM"

            elif "Esthesioneuroblastoma" in cancer_type:
                return "ESNB"

            elif "Sarcoma" in histotype or "Sarcoma" in cancer_type:
                return "SARC"

            elif "Pituitary tumor" in cancer_type:
                return "PITT"

            elif "Thyroid cancer" in cancer_type:
                return "THCA"

            elif "Cholangiocarcinoma (Bile duct cancer)" in cancer_type:
                return "CHOL"

            elif "Endocrine" in cancer_type:
                return "ECRI"

            elif "Urothelial cancer" in cancer_type:
                return "BLCA"

            elif "Prostate cancer" in cancer_type:
                return "PRAD"

            elif "Testicular cancer" in cancer_type:
                return "TECA"

            elif "Ovarian cancer" in cancer_type:
                return "GynOnc"  # "OVCA"
            elif "Ovarian  cancer" in cancer_type:
                return "GynOnc"  # "OVCA"
            elif "Serous ovarian carcinoma" in histotype:
                return "GynOnc"  # "OVCA"
            elif "Endometrial cancer" in cancer_type:
                return "GynOnc"  # "UCEC"
            elif "Uterine" in cancer_type:
                return "GynOnc"  # "UCEC"
            elif "Vaginal cancer" in cancer_type:
                return "GynOnc"  # "VACA"

            elif "Healthy donor" in cancer_type:
                return "hd"

            elif cancer_type == "N/A" and histotype == "N/A":
                return "unknown"
            elif cancer_type == "N/A":
                return histotype
            else:
                return cancer_type

        def get_platform(platform, center):
            if "Custom V2 Exome Bait, 48 RXN X 16 tubes" in platform:
                return "Agilent_Broad"
            elif "NimbleGen hg18 Exome v2" in platform:
                return "NimbleGen_hg18"
            elif "NimbleGen SeqCap EZ Human Exome Library v2.0" in platform:
                return "NimbleGen_v2"
            elif "NimbleGen SeqCap EZ Human Exome Library v3.0" in platform:
                return "NimbleGen_v3"
            elif "NimbleGen SeqCap EZ HGSC VCRome v2.1" in platform:
                return "NimbleGen_VCRome"
            elif "NimbleGen SeqCap EZ HGSC VCRome v2.1-PKv1" in platform:
                return "NimbleGen_VCRome-PKv1"
            elif "SureSelect Human All Exon 38 Mb v2" in platform:
                return "Agilent_SureSelect_38Mb"
            elif "SureSelect Human All Exon 50Mb Kit" in platform:
                return "Agilent_SureSelect_50Mb"
            elif platform == "Agilent":
                return "CCGD" if center == "CCGD" else "Agilent"
            elif platform == "ICE":
                return "ICE"
            elif platform == "TWIST":
                return "TWIST"
            elif "RNA" in platform:
                return "RNAseq"
            elif "Transcriptome" in platform:
                return "Transcr"
            elif "ULP" in platform:
                return "ULP"
            elif "HISEQ" in platform:
                return "HiSeq"
            elif "MiSeq" in platform:
                return "MiSeq"
            elif "Nova" in platform:
                return "NovaSeq"
            elif platform == "NA" and center == "NA":
                return "unknown"
            elif platform == "NA":
                return center
            else:
                return platform.replace(" ", "-")

        def get_sample_type(sample_type, sample_description):
            if isinstance(sample_type, float) or sample_type in ["Unknown", "unknown"]:
                return "unknown"
            elif sample_type == "Normal":
                return "N"
            elif "Primary" in sample_type:
                return "P"
            elif "Brain metastasis" in sample_type or "Brain Metastasis" in sample_type:
                return "BM"
            elif (
                sample_type == "Extracranial metastasis" or sample_type == "Metastatic"
            ):
                return "EM"
            elif sample_type == "Tumor":
                if isinstance(sample_description, str) and (
                    "cfDNA" in sample_description or sample_description == "cfDNa"
                ):
                    return "cfDNA"
                else:
                    return "T"
            else:
                return sample_type

        def get_sample_material(sample_material):
            if (
                isinstance(sample_material, float)
                or sample_material == "Unknown"
                or sample_material == "unknown"
            ):
                return "unknown"
            elif "FFPE" in sample_material:
                return "FFPE"
            elif "FF" in sample_material:
                return "FF"
            elif "Blood" in sample_material:
                return "Blood"
            elif "Buffycoat" in sample_material:
                return "Buffycoat"
            elif "Plasma" in sample_material:
                return "Plasma"
            elif "CSF" in sample_material:
                return "CSF"
            elif "Urine" in sample_material:
                return "Urine"
            else:
                return sample_material

        if self.cancer_type in sif.data.columns:
            sif.data.loc[sif.data[self.cancer_type].isna(), self.cancer_type] = "NA"

        if self.histotype in sif.data.columns:
            sif.data.loc[sif.data[self.histotype].isna(), self.histotype] = "NA"

        if self.histology not in sif.data.columns:
            sif.data[self.histology] = sif.data.apply(
                lambda s: get_histology(s.get(self.cancer_type, ""), s.get(self.histotype, "")),
                axis=1,
            )

        if self.center in sif.data.columns:
            sif.data.loc[sif.data[self.center].isna(), self.center] = "NA"

        if self.platform in sif.data.columns:
            sif.data.loc[sif.data[self.platform].isna(), self.platform] = "NA"

        if self.platform_abv not in sif.data.columns and (self.platform in sif.data.columns or self.center in sif.data.columns):
            sif.data[self.platform_abv] = sif.data.apply(
                lambda s: get_platform(s.get(self.platform, ""), s.get(self.center, "")),
                axis=1,
            )

        if self.data_type in sif.data.columns:
            sif.data.loc[sif.data[self.data_type].isna(), self.data_type] = "NA"

        if self.sample_type not in sif.data.columns:
            sif.data[self.sample_type] = sif.data.apply(
                lambda s: get_sample_type(
                    s.get(self.sample_type_long, ""), s.get(self.sample_description, "")
                ),
                axis=1,
            )

        sif.data[self.material] = sif.data.apply(
            lambda s: get_sample_material(s.get(self.material, "")), axis=1
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

        if inplace:
            self.data = sif.data
        else:
            return sif
