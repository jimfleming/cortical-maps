import numpy as np


def modify_synaptic_weight(rfON, rfOFF, weight_rf, ONOFFCrtxPlt, norm=1):
    # Modifying synaptic weights
    max_resp_on = np.sum(rfON)
    max_resp_off = np.sum(np.abs(rfOFF))

    if max_resp_on != 0 and max_resp_off != 0:
        if ONOFFCrtxPlt == 1:
            LWON = 1
            if norm:
                LWOFF = weight_rf * (max_resp_on / max_resp_off)
            else:
                LWOFF = weight_rf
        elif ONOFFCrtxPlt == -1:
            if norm:
                LWON = weight_rf * (max_resp_off / max_resp_on)
            else:
                LWON = weight_rf
            LWOFF = 1
    else:
        LWON = 1
        LWOFF = 1
        zero_polarity = 1  # to measure the number of cortical locations that non-dominant polarity is zero

    return LWON, LWOFF, zero_polarity
