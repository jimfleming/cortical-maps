import numpy as np


def make_colormap(CrtxGridXX):
    nRGB1 = CrtxGridXX.shape[0]
    NColor = nRGB1
    RetinotopyRFspace_plot = np.zeros((nRGB1, nRGB1))
    RetinotopyRFspace_plot[:] = (np.arange(nRGB1 * nRGB1) + 1).reshape((nRGB1, nRGB1))

    RedGrid, _ = np.meshgrid(np.linspace(0, 255, NColor), np.linspace(0, 0, NColor))
    _, GreenGrid = np.meshgrid(
        np.linspace(230, 230, NColor), np.linspace(230, 51, NColor)
    )
    BlueGrid, _ = np.meshgrid(np.linspace(0, 255, NColor), np.linspace(0, 0, NColor))

    cmap = np.dstack((RedGrid, GreenGrid, BlueGrid)) / 255
    return cmap, RetinotopyRFspace_plot
