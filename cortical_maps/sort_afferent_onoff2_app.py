import numpy as np
import scipy.ndimage
from dataclasses import dataclass
from tqdm.notebook import trange


@dataclass
class MapStepHolder:
    OD_holder: np.ndarray
    ONOFF_holder: np.ndarray
    retinotopy_holder: np.ndarray
    retinotopy_plot_holder: np.ndarray


def sort_afferent_onoff2_app(appdata):
    filt_sort = appdata.SortFiltONOFF
    onoff_input = appdata.ONOFF_ODSort
    od_input = appdata.OD_ODSort
    retinotopy_input = appdata.RetODsorted
    retinotopyPlot = appdata.RetinotopyODSortedPlot
    N_repeat = appdata.NSortONOFF
    rng_trial = appdata.rng_trial

    # The same functionality as fun_aff_sort but the ON OFF can only move
    # within the same OD band
    rng = np.random.default_rng(seed=rng_trial)

    onoff_sorted = onoff_input.copy()
    od_sorted = od_input.copy()
    ret_sorted = retinotopy_input.copy()

    OD_holder = []
    ONOFF_holder = []
    retinotopy_plot_holder = []
    retinotopy_holder = []
    for n in trange(N_repeat):
        #   Selection of points should be randomized
        Affr, Affc = np.nonzero(onoff_sorted != 0)
        ind_rand2 = rng.permutation(len(Affr))
        for q in range(len(Affr)):
            i1 = Affr[ind_rand2[q]]
            j1 = Affc[ind_rand2[q]]

            temp = np.zeros_like(onoff_sorted)
            temp[i1, j1] = 1
            temp_mask1 = scipy.ndimage.convolve(temp, filt_sort, mode="constant")
            # Convolution Center Contra
            Conv_CC = (onoff_sorted > 0) * temp_mask1 * (temp_mask1 > 0)
            Sum_CC = np.abs(Conv_CC).sum()

            # Convolution Center Ipsi
            Conv_CI = (onoff_sorted < 0) * temp_mask1 * (temp_mask1 > 0)
            Sum_CI = np.abs(Conv_CI).sum()

            # Convolution Surround Contra
            Conv_SC = (onoff_sorted > 0) * temp_mask1 * (temp_mask1 < 0)
            Sum_SC = np.abs(Conv_SC).sum()

            # Convolution Surround Ipsi
            Conv_SI = (onoff_sorted < 0) * temp_mask1 * (temp_mask1 < 0)
            Sum_SI = np.abs(Conv_SI).sum()

            if onoff_sorted[i1, j1] == 1:
                Conv_Result = (Sum_CC - Sum_CI) - (Sum_SC - Sum_SI)
                mask_search = (
                    scipy.ndimage.convolve(temp, np.ones((3, 3)), mode="constant")
                    - temp
                )
                if od_sorted[i1, j1] == 1:
                    srch_r, srch_c = np.nonzero(
                        mask_search.astype(bool)
                        & (onoff_sorted != 1)
                        & (od_sorted == 1)
                    )
                elif od_sorted[i1, j1] == -1:
                    srch_r, srch_c = np.nonzero(
                        mask_search.astype(bool)
                        & (onoff_sorted != 1)
                        & (od_sorted == -1)
                    )
                Conv_temp = np.zeros(len(srch_r))
                Aff_to_be_replaced = np.zeros(len(srch_r))
                for k1 in range(len(srch_r)):
                    temp_result = onoff_sorted.copy()
                    temp_result[i1, j1] = onoff_sorted[srch_r[k1], srch_c[k1]]
                    temp_result[srch_r[k1], srch_c[k1]] = onoff_sorted[i1, j1]
                    temp_mask_Contra = np.roll(
                        temp_mask1, (srch_r[k1] - i1, srch_c[k1] - j1), axis=(0, 1)
                    )

                    # Convolution Center Contra
                    Conv_CC_temp = np.sum(
                        (temp_result > 0) * temp_mask_Contra * (temp_mask_Contra > 0)
                    )
                    # Convolution Center Ipsi
                    Conv_CI_temp = np.sum(
                        (temp_result < 0) * temp_mask_Contra * (temp_mask_Contra > 0)
                    )
                    # Convolution Surround Contra
                    Conv_SC_temp = np.sum(
                        (temp_result > 0) * temp_mask_Contra * (temp_mask_Contra < 0)
                    )
                    # Convolution Surround Ipsi
                    Conv_SI_temp = np.sum(
                        (temp_result < 0) * temp_mask_Contra * (temp_mask_Contra < 0)
                    )
                    Conv_temp[k1] = (Conv_CC_temp - Conv_CI_temp) - (
                        Conv_SC_temp - Conv_SI_temp
                    )
                    Aff_to_be_replaced[k1] = onoff_sorted[srch_r[k1], srch_c[k1]]
            elif onoff_sorted[i1, j1] == -1:
                Conv_Result = (Sum_CI - Sum_CC) - (Sum_SI - Sum_SC)
                mask_search = (
                    scipy.ndimage.convolve(temp, np.ones((3, 3)), mode="constant")
                    - temp
                )
                if od_sorted[i1, j1] == 1:
                    srch_r, srch_c = np.nonzero(
                        mask_search.astype(bool)
                        & (onoff_sorted != -1)
                        & (od_sorted == 1)
                    )
                elif od_sorted[i1, j1] == -1:
                    srch_r, srch_c = np.nonzero(
                        mask_search.astype(bool)
                        & (onoff_sorted != -1)
                        & (od_sorted == -1)
                    )
                Conv_temp = np.zeros(len(srch_r))
                Aff_to_be_replaced = np.zeros(len(srch_r))
                for k1 in range(len(srch_r)):
                    temp_result = onoff_sorted.copy()
                    temp_result[i1, j1] = onoff_sorted[srch_r[k1], srch_c[k1]]
                    temp_result[srch_r[k1], srch_c[k1]] = onoff_sorted[i1, j1]
                    temp_mask_Contra = np.roll(
                        temp_mask1, (srch_r[k1] - i1, srch_c[k1] - j1), axis=(0, 1)
                    )

                    # Convolution Center Contra
                    Conv_CC_temp = np.sum(
                        (temp_result > 0) * temp_mask_Contra * (temp_mask_Contra > 0)
                    )
                    # Convolution Center Ipsi
                    Conv_CI_temp = np.sum(
                        (temp_result < 0) * temp_mask_Contra * (temp_mask_Contra > 0)
                    )
                    # Convolution Surround Contra
                    Conv_SC_temp = np.sum(
                        (temp_result > 0) * temp_mask_Contra * (temp_mask_Contra < 0)
                    )
                    # Convolution Surround Ipsi
                    Conv_SI_temp = np.sum(
                        (temp_result < 0) * temp_mask_Contra * (temp_mask_Contra < 0)
                    )
                    Conv_temp[k1] = (Conv_CI_temp - Conv_CC_temp) - (
                        Conv_SI_temp - Conv_SC_temp
                    )
                    Aff_to_be_replaced[k1] = onoff_sorted[srch_r[k1], srch_c[k1]]

            if np.sum(Conv_temp > Conv_Result):
                (ind_replace_vec,) = np.nonzero(Conv_temp == Conv_temp.max())
                ind_replace = ind_replace_vec[0]
                onoff_sorted[srch_r[ind_replace], srch_c[ind_replace]] = onoff_sorted[
                    i1, j1
                ]
                onoff_sorted[i1, j1] = Aff_to_be_replaced[ind_replace]

                Crtx_plt_ONOFF_OD_temp = od_sorted[i1, j1]
                od_sorted[i1, j1] = od_sorted[srch_r[ind_replace], srch_c[ind_replace]]
                od_sorted[
                    srch_r[ind_replace], srch_c[ind_replace]
                ] = Crtx_plt_ONOFF_OD_temp

                ret_temp = ret_sorted[srch_r[ind_replace], srch_c[ind_replace]]
                ret_sorted[srch_r[ind_replace], srch_c[ind_replace]] = ret_sorted[
                    i1, j1
                ]
                ret_sorted[i1, j1] = ret_temp

                ret_temp2 = retinotopyPlot[srch_r[ind_replace], srch_c[ind_replace]]
                retinotopyPlot[
                    srch_r[ind_replace], srch_c[ind_replace]
                ] = retinotopyPlot[i1, j1]
                retinotopyPlot[i1, j1] = ret_temp2

        OD_holder.append(od_sorted)
        ONOFF_holder.append(onoff_sorted)
        retinotopy_holder.append(ret_sorted)
        retinotopy_plot_holder.append(retinotopyPlot)

    map_step_holder = MapStepHolder(
        OD_holder=OD_holder,
        ONOFF_holder=ONOFF_holder,
        retinotopy_holder=retinotopy_holder,
        retinotopy_plot_holder=retinotopy_plot_holder,
    )

    appdata.ONOFFCrtxPlt = onoff_sorted
    appdata.ODCrtxPlt = od_sorted
    appdata.RetONOFFsorted = ret_sorted
    appdata.RetinotopyONOFFSortedPlot = retinotopyPlot
    appdata.ONOFFmap_step_holder = map_step_holder
