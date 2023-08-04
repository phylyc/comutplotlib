import numpy as np
import os
import pandas as pd
from typing import Callable, Optional, Union

import comutplotlib.pandas_util as pd_util


class AnnotationTable(object):
    @classmethod
    def from_file(
        cls,
        path_to_file: str,
        selection: Union[Callable[..., bool], dict] = None,
        complement: bool = False,
        encoding: str = None,
    ):
        if not os.path.exists(path_to_file):
            raise FileNotFoundError(
                f"No AnnotationTable instance created. "
                f"No such file or directory: '{path_to_file}'"
            )
        extension = os.path.splitext(path_to_file)[1]
        if extension == ".xls" or extension == ".xlsx":
            data = pd.read_excel(io=path_to_file)
        else:
            if encoding is not None:
                with open(path_to_file, encoding=encoding) as f:
                    data = pd.read_csv(
                        filepath_or_buffer=f, sep="\t", engine="c", comment="#"
                    )
            else:
                try:
                    with open(path_to_file, encoding="utf8") as f:
                        data = pd.read_csv(
                            filepath_or_buffer=f, sep="\t", engine="c", comment="#"
                        )
                except UnicodeError:
                    with open(path_to_file, encoding="latin1") as f:
                        data = pd.read_csv(
                            filepath_or_buffer=f, sep="\t", engine="c", comment="#"
                        )
        # Drop columns with no information:
        cols_to_drop = data.apply(
            lambda col: (
                len(col.unique()) == 1
                and (
                    np.isnan(col.unique()[0])
                    if isinstance(col.unique()[0], float)
                    else False
                )
            )
        )
        data.drop(cols_to_drop.loc[cols_to_drop].index, axis=1, inplace=True)
        return AnnotationTable(data=data, selection=selection, complement=complement)

    def __init__(
        self,
        data: pd.DataFrame,
        *args,
        selection: Union[Callable[..., bool], dict] = None,
        complement: bool = False,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.data = data
        self.selection = selection
        self.select(selection=selection, complement=complement, inplace=True)

    def constant_column(self, value):
        return pd.Series(data=value, index=self.data.index)

    @property
    def default_missing_column(self):
        return self.constant_column(value=pd.NA)

    def update_inplace(self, data: pd.DataFrame) -> None:
        self.data = data

    def select(
        self,
        selection: Union[Callable[..., bool], dict],
        complement: bool = False,
        inplace: bool = False,
    ) -> Optional["AnnotationTable"]:
        """
        :param selection: (dict) keys are column names in self,
            values are (a) a value in that column, (b) a list of values
            in that column, or (c) a callable object that returns a boolean.
            The list of values can also contain boolean functions.
        :param complement: (bool) if true, select for complement of selection
        :param inplace: (bool) if true, update data to selected data in place
        :return: (AnnotationTable) with selection of data
        """
        if inplace:
            selected_data = pd_util.select(
                data_frame=self.data, selection=selection, complement=complement
            )
            self.update_inplace(data=selected_data)
        else:
            if selection is None:
                return self
            else:
                return AnnotationTable(
                    data=self.data, selection=selection, complement=complement
                )

    def to_csv(self, path_to_file: str, **kwargs) -> None:
        # todo: add encodings
        self.data.to_csv(
            path_or_buf=path_to_file,
            header=True,
            index=False,
            sep="\t",
            na_rep="nan",
            mode="w+",
            **kwargs,
        )
