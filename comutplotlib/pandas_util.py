from functools import reduce
import numpy as np
import pandas as pd
from typing import Any, Callable, Dict, Iterable, List, Tuple, Union
import warnings


class Frame(object):
    """"""

    def __init__(self, *args: pd.Index):
        self.axes = list(args)

    @property
    def names(self):
        return [axis.name for axis in self.axes]


class PandasChainedAssignmentWarnHandler:
    """Context manager to temporarily set pandas chained assignment warning. Usage:

    with PandasChainedAssignmentWarnHandler():
         code

    with PandasChainedAssignmentWarnHandler('raise'):
         run my code and figure out which line causes the error!

    adapted from https://gist.github.com/notbanker/2be3ed34539c86e22ffdd88fd95ad8bc
    """

    def __init__(self, chained=None):
        acceptable = [None, "warn", "raise"]
        assert chained in acceptable, "chained must be in " + str(acceptable)
        self.swcw = chained

    def __enter__(self):
        self.saved_swcw = pd.options.mode.chained_assignment
        pd.options.mode.chained_assignment = self.swcw
        return self

    def __exit__(self, *args):
        pd.options.mode.chained_assignment = self.saved_swcw


def as_sparse_dtype(
    df: pd.DataFrame, dtype: str = None, precision: int = None, fill_value: Any = None
):
    """
    :param df:
    :param dtype: The dtype of the underlying array storing the non-fill
        value values. {int, float}
    :param precision: If None is give, the precision is determined by
        twice the maximum count.
    :param fill_value: The scalar value not stored in the SparseArray.
        By default, this depends on dtype.
            dtype       na_value
            float       np.nan
            int         0
            bool        False
            datetime64  pd.NaT
            timedelta64 pd.NaT
        The default value may be overridden by specifying a fill_value.
    """
    if dtype is None:
        if isinstance(df.dtypes, pd.Series):
            dtypes = df.dtypes.unique()
            if len(dtypes) > 1:
                return df
            else:
                _dtype = dtypes[0]
        else:
            _dtype = df.dtypes
    else:
        _dtype = dtype

    if "int" in _dtype.name:
        __dtype = "int"
    elif "float" in _dtype.name:
        __dtype = "float"
    else:
        raise TypeError

    if precision is None:
        if __dtype == "int":
            maximum = 2 * max(np.max(df.to_numpy()), 1)
            _precision = 2 ** int(
                np.ceil(np.log2(np.log2(maximum)))
            )  # get necessary precision
            _precision = max(8, min(_precision, 64))  # use nearest best available precision
        else:
            _precision = ""
    else:
        _precision = precision

    dtype = f"{__dtype}{_precision}"
    return df.astype(pd.SparseDtype(dtype=dtype, fill_value=fill_value))


def sparse_setitem(df, index, columns, val):
    """ https://stackoverflow.com/questions/49032856/assign-values-to-sparsearray-in-pandas
    Insert data in a DataFrame with SparseDtype format

    Only applicable for pandas version > 0.25

    Args
    ----
    df : DataFrame with series formatted with pd.SparseDtype
    index: str, or list, or slice object
        Same as one would use as first argument of .loc[]
    columns: str, list, or slice
        Same one would normally use as second argument of .loc[]
    val: insert values

    Returns
    -------
    df: DataFrame
        Modified DataFrame

    """

    # Save the original sparse format for reuse later
    spdtypes = df.dtypes[columns]
    # Convert concerned Series to dense format
    df[columns] = df[columns].sparse.to_dense()
    # Do a normal insertion with .loc[]
    df.loc[index, columns] = val
    # Back to the original sparse format
    df[columns] = df[columns].astype(spdtypes)
    return df


def stack(
    df: pd.DataFrame, level: Any = -1, dropna: Any = True, chunksize: int = 1000000
):
    if hasattr(df, "sparse") and df.columns.nlevels > 1:
        # Only dense dataframes are properly stacked. To preserve memory bounds,
        # stack dataframe in dense chunks.
        stacked = []
        batch_stratum = df.index
        num_batches = len(batch_stratum) // chunksize + 1
        if num_batches == 1:
            return stack(df.sparse.to_dense(), level=level, dropna=dropna)
        for batch in np.array_split(batch_stratum, num_batches):
            chunk = (
                subframe(df, batch=batch, batch_level=batch_stratum.names)
                if num_batches > 1
                else df
            )
            stacked.append(
                as_sparse_dtype(
                    df=chunk.sparse.to_dense().stack(level=level, dropna=dropna),
                    precision=64,
                )
            )
        return pd.concat(stacked)
    else:
        return df.stack(level=level, dropna=dropna)


def unstack(
    df: pd.DataFrame, level: Any = -1, fill_value: Any = None, chunksize: int = 1000000
):
    if hasattr(df, "sparse") and df.index.nlevels > 1:
        # Only dense dataframes are properly unstacked. To preserve memory bounds,
        # unstack dataframe in dense chunks.
        unstacked = []
        batch_stratum = df.index.droplevel(level=level).unique()
        num_batches = len(batch_stratum) // chunksize + 1
        for batch in np.array_split(batch_stratum, num_batches):
            chunk = (
                subframe(df, batch=batch, batch_level=batch_stratum.names)
                if num_batches > 1
                else df
            )
            unstacked.append(
                as_sparse_dtype(
                    df=chunk.sparse.to_dense().unstack(
                        level=level, fill_value=fill_value
                    ),
                    precision=64,
                )
            )
        return pd.concat(unstacked)
    else:
        return df.unstack(level=level, fill_value=fill_value)


def subframe(
    df: Union[pd.DataFrame, pd.Series],
    batch: Union[pd.Index, List[str]],
    batch_level: Union[str, List[str]],
    apply_along_batch_level=None,
) -> Union[pd.DataFrame, pd.Series]:
    if not isinstance(batch, pd.Index):
        batch = pd.Index(batch)

    # cast Index (of tuples) to 2darray
    try:
        _batch = pd.MultiIndex.from_tuples(batch).to_frame().to_numpy()
    except TypeError or ValueError:
        _batch = batch.to_frame().to_numpy()

    def get_batch_idx(names):
        batch_slice = []
        i = 0
        for level, name in enumerate(names):
            if name == batch_level or name in batch_level:
                batch_slice.append(np.unique(_batch[:, i]))
                i += 1
            else:
                batch_slice.append(slice(None))

        return tuple(batch_slice) if len(names) > 1 else batch

    def batch_level_is_in_levels(levels):
        if isinstance(batch_level, list):
            return all([level in levels for level in batch_level])
        else:
            return batch_level in levels

    batch_over_index = batch_level_is_in_levels(df.index.names)
    batch_over_columns = isinstance(df, pd.DataFrame) and batch_level_is_in_levels(
        df.columns.names
    )

    # todo: make selection faster
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=FutureWarning)
        # FutureWarning: The behavior of indexing on a MultiIndex with a nested
        # sequence of labels is deprecated and will change in a future version.
        # `series.loc[label, sequence]` will raise if any members of 'sequence'
        # or not present in the index's second level. To retain the old behavior,
        # use `series.index.isin(sequence, level=1)`
        if batch_over_index:
            batch_idx = get_batch_idx(df.index.names)
            _subframe = (
                df.loc[batch_idx, :]
                if isinstance(df, pd.DataFrame)
                else df.loc[batch_idx]
            )
            _axis = 0
        elif batch_over_columns:
            batch_idx = get_batch_idx(df.columns.names)
            _subframe = df.loc[:, batch_idx]
            _axis = 1
        else:
            _subframe = df
            _axis = -1

    return (
        _subframe.apply(apply_along_batch_level, axis=_axis)
        if apply_along_batch_level is not None and _axis in [0, 1]
        else _subframe
    )


def pool(df: pd.DataFrame, pool_as: Dict[str, Dict[str, List[str]]]) -> pd.DataFrame:
    if pool_as is None or len(pool_as) == 0:
        return df
    pooled_df = pd.DataFrame.from_dict(
        {
            key: pd.concat(
                [
                    subframe(
                        df=df,
                        batch_level=batch_level,
                        batch=batch,
                        apply_along_batch_level=sum,
                    )
                    # Assumption: There is only one pair in this dict.
                    for batch_level, batch in agg.items()
                ]
            )
            for key, agg in pool_as.items()
        }
    )
    if (
        df.index.names == pooled_df.columns.names
        or df.columns.names == pooled_df.index.names
    ):
        pooled_df = pooled_df.T
    pooled_df.index.names = df.index.names
    pooled_df.columns.names = df.columns.names
    return pooled_df


def to_instance(
    ndarray: np.ndarray, instance: Union[pd.DataFrame, pd.Series]
) -> Union[pd.Series, pd.DataFrame]:
    if np.prod(ndarray.shape) != np.prod(instance.shape):
        raise TypeError(f"got incompatible shapes: {ndarray.shape}, {instance.shape}.")

    if isinstance(instance, pd.Series):
        return pd.Series(ndarray.flatten(), index=instance.index)
    elif isinstance(instance, pd.DataFrame):
        return pd.DataFrame(
            ndarray.reshape((instance.shape[0], -1)),
            index=instance.index,
            columns=instance.columns,
        )
    else:
        raise TypeError(
            f"Instance of pandas.DataFrame or pandas.Series expected, "
            f"but {type(instance)} given."
        )


def to_ndarray(data_frame: Union[pd.DataFrame, pd.Series]) -> np.ndarray:
    if isinstance(data_frame.index, pd.MultiIndex):
        index_dimensions = [len(level) for level in data_frame.index.levels]
    else:
        index_dimensions = [data_frame.shape[0]]

    if isinstance(data_frame, pd.DataFrame):
        if isinstance(data_frame.columns, pd.MultiIndex):
            column_dimensions = [len(level) for level in data_frame.columns.levels]
        else:
            column_dimensions = [data_frame.shape[1]]
    else:
        column_dimensions = []

    return data_frame.to_numpy().reshape((*index_dimensions, *column_dimensions))


def select(
    data_frame: pd.DataFrame,
    selection: Union[
        Callable[[pd.Series], bool],
        Dict[
            str,
            Union[
                int,
                float,
                str,
                Callable[..., bool],
                Iterable[Union[int, float, str, Callable[..., bool]]],
            ],
        ],
        None,
    ],
    complement: bool = False,
) -> pd.DataFrame:
    """
    :param data_frame: (pandas.DataFrame)
    :param selection:
        (A) callable object on a row of the data_frame that returns a boolean
        (B) (dict) with keys and values specified in get_mask()
    :param complement: (bool) if true, select for complement of selection
    :return: (pandas.DataFrame) selection of dataframe
    """
    if selection is None:
        return data_frame

    if not len(data_frame.index):  # data_frame is empty
        return data_frame

    # (-)
    if selection:
        # (A)
        if callable(selection):
            mask = data_frame.apply(selection, axis=1)
        # (B)
        else:
            mask = pd.Series(data=True, index=data_frame.index)
            for _key, _value in selection.items():
                mask &= get_mask(data_frame, _key, _value)
    else:
        mask = pd.Series(data=False, index=data_frame.index)

    if complement:
        mask = ~mask
    return data_frame.loc[mask]


def get_mask(
    data_frame: pd.DataFrame,
    key: Union[str, Tuple[str]],
    value: Union[
        int,
        float,
        str,
        Callable[..., bool],
        Iterable[Union[int, float, str, Callable[..., bool]]],
    ],
) -> pd.Series:
    """
    :param data_frame: pandas DataFrame
    :param key: column name or tuple of column names in data_frame
    :param value: one of
        (a) a callable object that returns a boolean, or
        (b) a value in that column that is selected for, or
        (c) a list of (a) and (b).
    :return: (pandas.Series)
    """

    def _apply(data, columns, func):
        if isinstance(columns, list):  # apply acts on DataFrame
            return data[columns].apply(func, axis=1)
        else:  # apply acts on Series
            return data[columns].apply(func)

    if isinstance(key, tuple):
        key = list(key)
    # (a)
    if callable(value):
        mask = _apply(data=data_frame, columns=key, func=value)
    # (c)
    elif isinstance(value, (list, set, tuple, np.ndarray, pd.Index)):
        if len(value):
            mask = pd.Series(data=True, index=data_frame.index)

            values = [v for v in value if not callable(v)]
            if values:
                mask &= data_frame[key].isin(values)

            callables = [c for c in value if callable(c)]
            for c in callables:
                mask &= _apply(data=data_frame, columns=key, func=c)
        else:
            mask = pd.Series(data=False, index=data_frame.index)
    # (b)
    else:
        mask = data_frame[key] == value

    return mask


def merge_selections(*selections: Union[dict, None]) -> dict:
    return reduce(__merge_selections, selections) if selections else None


def __merge_selections(
    selection_1: Union[dict, None], selection_2: Union[dict, None]
) -> Dict[Any, list]:
    """
    :param selection_1: dict with lists as values, or None
    :param selection_2: dict with lists as values, or None
    :return: dict
    """
    if selection_1 is None:
        return selection_2
    if selection_2 is None:
        return selection_1

    def set_of_values(_values):
        if isinstance(_values, str):
            return {_values}
        try:
            return set(_values)
        except TypeError:
            return {_values}

    merged_selection = {}
    keys = set(selection_1.keys()) | set(selection_2.keys())
    for key in keys:
        values = set()
        values |= set_of_values(selection_1.get(key, values))
        values |= set_of_values(selection_2.get(key, values))
        merged_selection[key] = list(values)
    return merged_selection


def position_outside_intervals(
    positions: pd.Series,
    interval_starts: np.array,
    interval_ends: np.array,
    include=False,
) -> pd.Series:
    # More intuitive solution, but complexity scales worse than below:
    # (is 5x slower for ~5MB, 50x slower for ~2.5GB data)
    #
    # interval_index = pd.IntervalIndex.from_arrays(
    #     left=interval_starts, right=interval_ends, closed="both",
    # )
    # return positions.transform(
    #     lambda start_pos: not interval_index.contains(start_pos).any(axis=None)
    # )

    # This solution is of O(n) complexity:
    # https://stackoverflow.com/questions/38201057/efficiently-check-if-value-is-present-in-any-of-given-ranges
    #
    starts = pd.DataFrame(data={"start": 1}, index=interval_starts)
    ends = pd.DataFrame(data={"end": -1}, index=interval_ends)
    transitions = pd.merge(
        starts, ends, how="outer", left_index=True, right_index=True
    ).fillna(0)
    # This dataframe stores per position the type of transitions.
    # Now, we need to know at each position if we are in an interval or not.
    # This can be done by counting the opening & closing parenthesis.
    transitions["transition"] = (
        transitions.pop("end") + transitions.pop("start")
    ).cumsum()
    # If the value is strictly greater than 0, the position (index) in a
    # right-open interval. This handles overlapping intervals.

    positions_df = pd.DataFrame(index=positions.drop_duplicates())
    position_in_list = (
        pd.merge(
            transitions,
            positions_df,
            how="outer",
            left_index=True,
            right_index=True,
        )
        .fillna(method="ffill")
        .fillna(0)  # replace NANs for positions that are before first interval
        .loc[positions_df.index]
        .get("transition")
        .astype(bool)
    )
    positions_in_interval = position_in_list.loc[position_in_list].index.to_list()
    selected_positions = (
        positions_in_interval
        if include
        else positions_df.drop(positions_in_interval).index
    )
    return positions.apply(lambda pos: pos in selected_positions)
