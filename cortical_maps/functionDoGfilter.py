import numpy as np
from scipy.stats import multivariate_normal


def functionDoGfilter(c1, c2, s1, s2, CSRatio):
    # using multivariate 2D mask inhibitory and excitatory effect
    # mask of different angles affecting same input

    # input parameter
    #   Mask1
    center_width1 = c1
    center_width2 = c2
    #   Mask2
    surround_width1 = s1
    surround_width2 = s2

    orientation = 0  # orientation in degree, the masks parameters must be asymmetric

    mu = np.array([0, 0])  # center of mask
    rx1 = 20  # 15
    rx2 = 20  # 15
    x1 = np.arange(-rx1, rx1 + 1)
    x2 = np.arange(-rx2, rx2 + 1)
    X1, X2 = np.meshgrid(x1, -x2)

    variance1 = np.array(
        [
            [center_width1, 0],
            [0, center_width2],
        ]
    )
    theta1 = orientation * (np.pi / 180)
    rotation_mat1 = np.array(
        [
            [-np.cos(theta1), np.sin(theta1)],
            [np.sin(theta1), np.cos(theta1)],
        ]
    )
    covariance1 = rotation_mat1 @ variance1 @ rotation_mat1.T

    F1 = multivariate_normal.pdf(
        np.stack([X1.flat, X2.flat], axis=-1), mean=mu, cov=covariance1
    )
    F1 = np.reshape(F1, (len(x2), len(x1)))

    variance2 = np.array(
        [
            [surround_width1, 0],
            [0, surround_width2],
        ]
    )
    theta2 = orientation * (np.pi / 180)
    rotation_mat2 = np.array(
        [
            [-np.cos(theta2), np.sin(theta2)],
            [np.sin(theta2), np.cos(theta2)],
        ]
    )
    covariance2 = rotation_mat2 @ variance2 @ rotation_mat2.T

    F2 = multivariate_normal.pdf(
        np.stack([X1.flat, X2.flat], axis=-1), mean=mu, cov=covariance2
    )
    F2 = np.reshape(F2, (len(x2), len(x1)))
    mask1 = CSRatio * F1 - F2
    return mask1
