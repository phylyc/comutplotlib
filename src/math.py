from itertools import pairwise
import numpy as np


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
