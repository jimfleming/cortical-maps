import numpy as np
import scipy.ndimage
import skimage.transform
from cortical_maps.fspecial import gaussian_fspecial


def smooth_map_app(appdata):
    # Self: Blurs
    ODCrtxPlt = appdata.ODCrtxPlt
    ONOFFCrtxPlt = appdata.ONOFFCrtxPlt
    n_interp = appdata.n_interpol

    # OD/ONOFF smoothing

    kern_std = 2  # The std of kernel for smoothing the OD map

    kern = gaussian_fspecial((kern_std * 6, kern_std * 6), kern_std)

    if appdata.aff_sampling_density < 1:
        kern = gaussian_fspecial((1, 1), 1)

    ODCrtxPlt_smooth = scipy.ndimage.convolve(ODCrtxPlt, kern)
    ODCrtxPlt_smooth = (ODCrtxPlt_smooth - np.amin(ODCrtxPlt_smooth)) / (
        np.amax(ODCrtxPlt_smooth) - np.amin(ODCrtxPlt_smooth)
    )
    ODCrtxPlt_interpolated = skimage.transform.rescale(
        ODCrtxPlt_smooth, n_interp, order=3
    )

    ONOFF_smoothed = scipy.ndimage.convolve(ONOFFCrtxPlt, kern)
    ONOFF_smoothed = (ONOFF_smoothed - np.amin(ONOFF_smoothed)) / (
        np.amax(ONOFF_smoothed) - np.amin(ONOFF_smoothed)
    )
    ONOFF_interpolated = skimage.transform.rescale(ONOFF_smoothed, n_interp, order=3)

    appdata.ONOFF_smoothed = ONOFF_smoothed
    appdata.ODCrtxPlt_smooth = ODCrtxPlt_smooth
    appdata.ONOFF_interpolated = ONOFF_interpolated
    appdata.ODCrtxPlt_interpolated = ODCrtxPlt_interpolated
