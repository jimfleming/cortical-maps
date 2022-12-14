import numpy as np


def gaussian_fspecial(shape=(3, 3), sigma=0.5):
    m, n = [(ss - 1.0) / 2.0 for ss in shape]
    y, x = np.ogrid[-m : m + 1, -n : n + 1]
    h = np.exp(-(x * x + y * y) / (2.0 * sigma * sigma))
    h[h < np.finfo(h.dtype).eps * h.max()] = 0
    hsum = h.sum()
    if hsum != 0:
        h /= hsum
    return h
