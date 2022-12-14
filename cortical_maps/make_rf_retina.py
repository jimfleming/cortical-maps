import numpy as np
import scipy.ndimage
from tqdm import trange
from cortical_maps.fspecial import gaussian_fspecial


def make_rf_retina(
    TotalRetinaCells, RFONCenter, RFOFFCenter, IndRetinaRfspace, onoff_rf_space
):
    # Make receptive field for each afferent coming to cortex based their
    # locationn and polarities in retina

    RetinaRF = []

    # 19 is the size of the filter it could be width of cortex matrix
    filt_rf_retina_on = gaussian_fspecial((60, 60), RFONCenter)
    filt_rf_retina_off = gaussian_fspecial((60, 60), RFOFFCenter)

    RFONfilt = filt_rf_retina_on
    RFOFFfilt = filt_rf_retina_off
    for jj in trange(TotalRetinaCells):
        tempGridFilter = IndRetinaRfspace == jj
        PolarityRFSpace = onoff_rf_space[tempGridFilter]
        import ipdb

        ipdb.set_trace()
        if PolarityRFSpace == 1:
            tmpRFON = scipy.ndimage.convolve(tempGridFilter, RFONfilt)
            tmpRFONOFF = tmpRFON / np.amax(np.abs(tmpRFON))
        elif PolarityRFSpace == -1:
            tmpRFOFF = scipy.ndimage.convolve(tempGridFilter, RFOFFfilt)
            tmpRFONOFF = (-1 * tmpRFOFF) / np.amax(np.abs(tmpRFOFF))
        RetinaRF.append(tmpRFONOFF)
    return RetinaRF
