from copy import deepcopy
from functools import reduce
import numpy as np
import os
import pandas as pd
from tqdm import tqdm
from typing import Any, Callable, Dict, Optional, Union, Literal
import warnings

from src.annotation_table import AnnotationTable
from src.maf_encoding import MAFEncoding
from src.mutation_annotation import MutationAnnotation
import src.pandas_util as pd_util


def join_mafs(mafs: list["MAF"]):
    if len(mafs) == 0:
        return MAF()
    elif len(mafs) == 1:
        return mafs[0].copy()
    else:
        maf = mafs[0]
        for other in tqdm(mafs[1:], desc="Joining MAFs", total=len(mafs), initial=1):
            maf = maf.join(other=other)
        return maf


class MAF(MutationAnnotation, AnnotationTable):
    @classmethod
    def from_file(
        cls,
        path_to_file: str,
        selection: Union[Callable[..., bool], Dict[str, Any]] = None,
        complement: bool = False,
        usecols: list = None,
        source: str = MAFEncoding.funcotator,
        ignore_column_requirements: bool = False,
        warn_if_empty: bool = False,
    ):
        """
        :param path_to_file:
        :param selection:
        :param complement:
        :param usecols: list of MAF default column names
        :param source: {"Oncotator", "Funcotator", "Tumorportal"};
            brings data frame into default format defined in MutationAnnotation.
        :param ignore_column_requirements: bool
        """
        if not os.path.exists(path_to_file):
            raise FileNotFoundError(
                f"No MAF instance created. No such file or directory: '{path_to_file}'"
            )

        encoding = MAFEncoding(from_source=source)
        data = cls.__load(path_to_file=path_to_file, encoding=encoding, usecols=usecols)
        maf = MAF(
            data=data,
            selection=selection,
            complement=complement,
            file=path_to_file,
            name=os.path.splitext(os.path.basename(path_to_file))[0],
            encoding=encoding,
            ignore_column_requirements=ignore_column_requirements,
        )
        if warn_if_empty and maf.empty:
            warnings.warn(f"MAF loaded from {path_to_file} is empty!")
        return maf

    @classmethod
    def __load(
        cls, path_to_file: str, encoding: MAFEncoding, usecols: list = None
    ) -> pd.DataFrame:
        """Load file and return a pandas.DataFrame.

        Change source encoding to default encoding. Check requirements on columns.

        :param path_to_file:
        :param encoding:
        :param usecols:
        :return: loaded pandas.DataFrame
        """
        dtype = {
            encoding.column_names_from_default.get(column, column): column_dtype
            for column, column_dtype in cls.column_dtype.items()
        }
        data = pd.read_csv(
            filepath_or_buffer=path_to_file,
            sep="\t",
            engine="c",
            comment="#",
            usecols=usecols,
            dtype=dtype,
            low_memory=False,
        )
        return data

    def __init__(
        self,
        data: pd.DataFrame = None,
        selection: Union[Callable[..., bool], dict] = None,
        complement: bool = False,
        file: str = None,
        name: str = None,
        encoding: MAFEncoding = None,
        ignore_column_requirements: bool = False,
    ) -> None:
        """
        :param data:
        :param selection:
        :param complement:
        :param file: file path
        :param name: name of file or MAF instance
        :param encoding:
        :param ignore_column_requirements:
        """
        super().__init__(data=data)
        if data is None:
            self.data = pd.DataFrame(data=None, columns=self.default_columns)
        if encoding is not None:
            encoding.to_default(data=self.data)
        self.file = file
        self.name = name
        self.selection = selection
        self.select(selection=selection, complement=complement, inplace=True)
        self.ignore_column_requirements = ignore_column_requirements
        if not self.ignore_column_requirements:
            self.check_columns()

    def check_columns(self) -> None:
        absent_columns = set(self.default_columns) - set(self.data.columns)
        if absent_columns:
            raise KeyError(
                f"Columns expected but not found: The file '{self.file}' "
                f"is not complete. Please make sure that the file has all the "
                f"required columns.. The columns '{absent_columns}' are missing."
            )

    def add_required_columns(
        self, verbose: bool = True, inplace: bool = False
    ) -> Optional["MAF"]:
        _self = self.copy()
        if _self.sample not in _self.data.columns:
            _self.add_patient_as_sample()  # required column.
        if _self.end_pos not in _self.data.columns:
            _self.add_end_position()  # not a required column.
        if inplace:
            self.update_inplace(_self.data)
        else:
            return _self

    def assign_column(self, name: str, value, inplace: bool = False):
        if inplace:
            self.update_inplace(data=self.data.assign(**{name: value}))
        else:
            _self = self.copy()
            _self.update_inplace(data=_self.data.assign(**{name: value}))
            return _self

    def add_patient_as_sample(self) -> None:
        self.assign_column(self.sample, self.data.get(self.patient), inplace=True)

    def add_end_position(self) -> None:
        end_pos = self.data[self.start_pos] + self.data[self.ref_allele].apply(len) - 1
        self.assign_column(self.end_pos, end_pos, inplace=True)

    @property
    def empty(self) -> bool:
        return self.num_variants == 0

    @property
    def num_variants(self) -> int:
        return self.data.shape[0]

    @property
    def patients(self) -> np.ndarray:
        patients = self.data.get(self.patient)
        if patients is not None:
            return patients.unique()
        else:
            return np.array([])

    @property
    def num_patients(self) -> int:
        return self.patients.shape[0]

    @property
    def samples(self) -> np.ndarray:
        samples = self.data.get(self.sample)
        if samples is not None:
            return samples.unique()
        else:
            return np.array([])

    @property
    def num_samples(self) -> int:
        return self.samples.shape[0]

    @property
    def effects(self) -> np.ndarray:
        return self.data.get(self.effect, self.default_missing_column).unique()

    @property
    def num_effects(self) -> int:
        return self.effects.shape[0]

    def copy(self) -> "MAF":
        maf = MAF()
        for name, value in self.__dict__.items():
            maf.__setattr__(name, deepcopy(value))
        return maf

    def join(self, other: "MAF"):
        if other is None:
            return self
        else:
            merged_data = pd.concat(
                [self.data, other.data], ignore_index=True
            ).drop_duplicates(ignore_index=True)
            maf = MAF(
                data=merged_data,
                name=".".join([maf.name for maf in [self, other] if maf.name]),
                ignore_column_requirements=(
                    self.ignore_column_requirements and other.ignore_column_requirements
                ),
            )
            maf.selection = pd_util.merge_selections(self.selection, other.selection)
            return maf

    ####################################################################
    # VARIANT FILTERS
    ####################################################################

    def drop_patient_duplicates(
        self,
        subset=None,
        keep: Literal["first", "last", False] = "first",
        inplace: bool = False,
    ):
        default_columns = [
            self.patient,
            self.gene_id,
            self.gene_name,
            self.chromosome,
            self.start_pos,
            self.strand,
            self.ref_allele,
            self.alt_allele,
        ]
        _subset = [label for label in self.data.columns if label in default_columns]
        subset = _subset if subset is None else subset
        if inplace:
            return self.data.drop_duplicates(subset=subset, keep=keep, inplace=inplace)
        else:
            new_data = self.data.drop_duplicates(subset=subset, keep=keep)
            _maf = self.copy()
            _maf.data = new_data
            return _maf

    ####################################################################
    # POOL & SELECT & STRATIFY
    ####################################################################

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
            pooled_maf = self.copy()
            pooled_maf.data = data
            return pooled_maf

    def select(
        self,
        selection: Union[Callable[..., bool], dict],
        complement: bool = False,
        inplace: bool = False,
    ) -> Optional["MAF"]:
        """
        :param selection: (dict) keys are column names in self,
            values are (a) a value in that column, (b) a list of values
            in that column, or (c) a callable object that returns a boolean.
            The list of values can also contain boolean functions.
        :param complement: (bool) if true, select for complement of selection
        :param inplace: (bool) if true, update data to selected data in place
        :return: (VirtualMAF) with selection of data
        """
        if inplace:
            super().select(selection=selection, complement=complement, inplace=inplace)
        else:
            if selection is None:
                return self
            else:
                return MAF(
                    data=self.data,
                    selection=selection,
                    complement=complement,
                    file=self.file,
                    name=self.name,
                    ignore_column_requirements=self.ignore_column_requirements,
                )

    def stratify_mutation_counts(
        self,
        by: str | list[str],
        average: list[str] = (),
        expand_index: bool = True,
        counts_column_name: str | None = "count",
        batch_by_patient: bool = False,
    ) -> pd.DataFrame | pd.Series:
        """
        :param by: key or list of keys in self.data.columns
        :param average: (optional) list of keys in self.data.columns
            which are averaged over each stratum
        :param expand_index: if True, complete multi index to Cartesian product
        :param counts_column_name: column name of mutation counts. If None,
            returns a pandas.Series
        :param batch_by_patient: batch stratification over patients to save
            memory resources.
        :returns:
        """
        dtype = pd.SparseDtype(dtype="int64", fill_value=0)
        if (
            isinstance(by, list)
            and len(by) > 2
            and self.patient in by
            and batch_by_patient
        ):
            counts = []
            _by = [stratum for stratum in by if stratum != self.patient]
            index = None
            for patient in tqdm(self.patients, desc="Patients"):
                _counts = self.select(
                    selection={self.patient: patient}
                ).stratify_mutation_counts(
                    by=_by,
                    average=average,
                    expand_index=expand_index,
                    counts_column_name=counts_column_name,
                )
                _counts = pd.concat({patient: _counts}, names=[self.patient], axis=1)
                index = (
                    index.union(_counts.index, sort=False)
                    if index is not None
                    else _counts.index
                )
                counts.append(_counts)
            index = index.sort_values()
            counts_reindexed = [
                c.reindex(index=index).fillna(0).astype(dtype) for c in counts
            ]  # todo: this takes too much time for a lot of patients!
            stratification = reduce(lambda a, b: a.join(b), counts_reindexed)
            return stratification.stack(level=self.patient)

        if counts_column_name is None:
            counts_column_name = "None"

        stratification = self.data.groupby(by=by, sort=False, dropna=False)
        available_count_columns = [
            col for col in self.data.columns if col not in by and col not in average
        ]
        counts_column = available_count_columns[0]
        agg_dict = dict(
            **{counts_column_name: pd.NamedAgg(column=counts_column, aggfunc="size")},
            **{key: pd.NamedAgg(column=key, aggfunc="mean") for key in average},
        )
        mutation_counts = stratification.agg(**agg_dict)

        if mutation_counts.empty:
            if counts_column_name == "None":
                mutation_counts = mutation_counts[counts_column_name]
                mutation_counts.name = None
            return mutation_counts

        mutation_counts = mutation_counts.astype({counts_column_name: dtype})

        if expand_index and isinstance(mutation_counts.index, pd.MultiIndex):
            # expand multi index to full Cartesian product
            index = pd.MultiIndex.from_product(mutation_counts.index.levels)
            mutation_counts = mutation_counts.reindex(index=index)
            # replace NAN with 0 for the counts:
            pd_util.sparse_setitem(
                df=mutation_counts,
                index=mutation_counts[counts_column_name].isna(),
                columns=[counts_column_name],
                val=0
            )
            # cast back to lower precision int dtype
            mutation_counts = mutation_counts.astype({counts_column_name: dtype})
            if average:
                # if there are columns that the aggregation averaged over,
                # then broadcast missing values across each gene:
                mutation_counts = (
                    mutation_counts.unstack(level=-1)  # map gene stratum to columns
                    .fillna(method="ffill", axis=0)
                    .fillna(method="bfill", axis=0)
                    .stack(level=-1)  # map gene stratum back to index
                )

        if counts_column_name == "None" and not average:
            mutation_counts = mutation_counts[counts_column_name]
            mutation_counts.name = None

        return mutation_counts
