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

    if isinstance(mean, pd.Series):
        lo = pd.Series(st.norm.ppf(alpha_ci / 2, loc=mean, scale=np.sqrt(var)), index=mean.index)
        hi = pd.Series(st.norm.ppf(1 - alpha_ci / 2, loc=mean, scale=np.sqrt(var)), index=mean.index)
    elif isinstance(mean, pd.DataFrame):
        lo = pd.DataFrame(st.norm.ppf(alpha_ci / 2, loc=mean, scale=np.sqrt(var)), index=mean.index, columns=mean.columns)
        hi = pd.DataFrame(st.norm.ppf(1 - alpha_ci / 2, loc=mean, scale=np.sqrt(var)), index=mean.index, columns=mean.columns)

    is_present = (n > 0) & (m > 0) & ~n.isna() & ~m.isna()

    return {
        "mean": mean.mask(~is_present),
        "lo": lo.mask(~is_present),
        "hi": hi.mask(~is_present),
        "fdr": alpha_ci,
        "base": base,
        # "p": p
    }
