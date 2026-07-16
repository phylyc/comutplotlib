from itertools import pairwise
import numpy as np
import pandas as pd
import scipy.stats as st


def decompose_rectangle_into_polygons(num, theta=0, pad=0):
    x_coords = [0, 1]
    y_coords = np.linspace(start=0, stop=1, num=num + 1)
    coords = []
    for pair in pairwise(y_coords):
        polygon_coords = []
        pair = (pair[0] + pad / 2, pair[1] - pad / 2)
        for y in pair:
            for x in x_coords:
                polygon_coords.append([x, y])
        corner_order = [0, 1, 3, 2]
        coords.append([polygon_coords[i] for i in corner_order])
    return np.array(coords)


def fold_change(n, N, m, M, alpha_ci=0.2, base=2, method="fisher"):
    # table = [[n, N - n],
    #          [m, M - m]]
    # if method == "fisher":
    #     _, p = st.fisher_exact(table, alternative="two-sided")
    #     return p
    # elif method == "barnard":
    #     p = st.barnard_exact(table, alternative="two-sided").pvalue
    # elif method == "boschloo":
    #     p = st.boschloo_exact(table, alternative="two-sided").pvalue
    # else:
    #     raise ValueError("method must be fisher/barnard/boschloo")

    # Jeffrey's prior
    p_n = (n + 0.5) / (N + 1)
    p_m = (m + 0.5) / (M + 1)
    mean = np.log(p_n / p_m) / np.log(base)
    var = ((1 / p_n - 1) / N + (1 / p_m - 1) / M + 1e-12) / np.log(base) ** 2

    scale = np.sqrt(var)
    lo_raw = st.norm.ppf(alpha_ci / 2, loc=mean, scale=scale)
    hi_raw = st.norm.ppf(1 - alpha_ci / 2, loc=mean, scale=scale)

    if isinstance(mean, pd.Series):
        lo = pd.Series(lo_raw, index=mean.index)
        hi = pd.Series(hi_raw, index=mean.index)
    elif isinstance(mean, pd.DataFrame):
        lo = pd.DataFrame(lo_raw, index=mean.index, columns=mean.columns)
        hi = pd.DataFrame(hi_raw, index=mean.index, columns=mean.columns)
    else:  # scalar or numpy array
        lo, hi = lo_raw, hi_raw

    is_present = (n > 0) & (m > 0) & pd.notna(n) & pd.notna(m)

    def _mask_absent(value):
        # set entries where the mutation is absent in either cohort to NaN
        if isinstance(value, (pd.Series, pd.DataFrame)):
            return value.mask(~is_present)
        return np.where(is_present, value, np.nan)

    return {
        "mean": _mask_absent(mean),
        "lo": _mask_absent(lo),
        "hi": _mask_absent(hi),
        "fdr": alpha_ci,
        "base": base,
        # "p": p
    }

