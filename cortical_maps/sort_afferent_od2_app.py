import numpy as np
import scipy.signal
import scipy.ndimage
from dataclasses import dataclass
from tqdm.notebook import trange


@dataclass
class MapStepHolder:
    OD_holder: np.ndarray
    ONOFF_holder: np.ndarray
    retinotopy_holder: np.ndarray
    retinotopy_plot_holder: np.ndarray


def sort_afferent_od2_app(appdata):
    filt_sort = appdata.SortFiltOD
    Crtx_plt_OD = appdata.CrtxOD3m
    Crtx_plt_ONOFF = appdata.CrtxONOFF3m
    ret_initial = appdata.Retinotopy3mIndex
    retinotopyPlot = appdata.Retinotopy3mIndPlot
    N_repeat = appdata.NSortOD
    rng_trial = appdata.rng_trial

    rng = np.random.default_rng(seed=rng_trial)

    result_OD_sorted = Crtx_plt_OD.copy()
    Crtx_plt_ONOFF_OD_sorted = Crtx_plt_ONOFF.copy()
    ret_sorted = ret_initial.copy()

    OD_holder = []
    ONOFF_holder = []
    retinotopy_holder = []
    retinotopy_plot_holder = []
    for n in trange(N_repeat):
        # Selection of points should be randomized
        Affr, Affc = np.where(result_OD_sorted != 0)
        ind_rand2 = rng.permutation(len(Affr))

        # loop through each coordinate
        for q in range(len(Affr)):
            i1 = Affr[ind_rand2[q]]
            j1 = Affc[ind_rand2[q]]

            # use filt_sort to build a search area around the coordinate
            temp = np.zeros_like(result_OD_sorted)
            temp[i1, j1] = 1
            temp_mask1 = scipy.signal.convolve2d(temp, filt_sort, mode="same")

            # build mask of valid points
            # Convolution Center/Surround Contra/Ipsi
            Conv_CC = (result_OD_sorted > 0) * temp_mask1 * (temp_mask1 > 0)
            Conv_CI = (result_OD_sorted < 0) * temp_mask1 * (temp_mask1 > 0)
            Conv_SC = (result_OD_sorted > 0) * temp_mask1 * (temp_mask1 < 0)
            Conv_SI = (result_OD_sorted < 0) * temp_mask1 * (temp_mask1 < 0)

            # sum valid points
            Sum_CC = np.sum(np.abs(Conv_CC))
            Sum_CI = np.sum(np.abs(Conv_CI))
            Sum_SC = np.sum(np.abs(Conv_SC))
            Sum_SI = np.sum(np.abs(Conv_SI))

            # contra
            if result_OD_sorted[i1, j1] == 1:
                # score sums for contra
                Conv_Result = (Sum_CC - Sum_CI) - (Sum_SC - Sum_SI)

                # neighbor search area
                mask_search = (
                    scipy.signal.convolve2d(temp, np.ones((3, 3)), mode="same") - temp
                )
                # candidate neighbor coordinates
                srch_r, srch_c = np.where(
                    mask_search.astype(bool) & (result_OD_sorted != 1)
                )

                Conv_temp = np.zeros(len(srch_r))
                Aff_to_be_replaced = np.zeros(len(srch_r))

                # for each candidate
                for k1 in range(len(srch_r)):
                    # swap current coordinate with candidate
                    temp_result = result_OD_sorted.copy()
                    temp_result[i1, j1] = result_OD_sorted[srch_r[k1], srch_c[k1]]
                    temp_result[srch_r[k1], srch_c[k1]] = result_OD_sorted[i1, j1]

                    # roll to new position
                    temp_mask_Contra = np.roll(
                        temp_mask1, (srch_r[k1] - i1, srch_c[k1] - j1), axis=(0, 1)
                    )

                    # calculate result again
                    Conv_CC_temp = np.sum(
                        (temp_result > 0) * temp_mask_Contra * (temp_mask_Contra > 0)
                    )
                    Conv_CI_temp = np.sum(
                        (temp_result < 0) * temp_mask_Contra * (temp_mask_Contra > 0)
                    )
                    Conv_SC_temp = np.sum(
                        (temp_result > 0) * temp_mask_Contra * (temp_mask_Contra < 0)
                    )
                    Conv_SI_temp = np.sum(
                        (temp_result < 0) * temp_mask_Contra * (temp_mask_Contra < 0)
                    )
                    Conv_temp[k1] = (Conv_CC_temp - Conv_CI_temp) - (
                        Conv_SC_temp - Conv_SI_temp
                    )
                    Aff_to_be_replaced[k1] = result_OD_sorted[srch_r[k1], srch_c[k1]]

            # ipsi
            elif result_OD_sorted[i1, j1] == -1:
                Conv_Result = (Sum_CI - Sum_CC) - (Sum_SI - Sum_SC)

                mask_search = (
                    scipy.signal.convolve2d(temp, np.ones((3, 3)), mode="same") - temp
                )
                srch_r, srch_c = np.where(
                    mask_search.astype(bool) & (result_OD_sorted != -1)
                )

                Conv_temp = np.zeros(len(srch_r))
                Aff_to_be_replaced = np.zeros(len(srch_r))

                for k1 in range(len(srch_r)):
                    temp_result = result_OD_sorted.copy()
                    temp_result[i1, j1] = result_OD_sorted[srch_r[k1], srch_c[k1]]
                    temp_result[srch_r[k1], srch_c[k1]] = result_OD_sorted[i1, j1]

                    temp_mask_Ipsi = np.roll(
                        temp_mask1, (srch_r[k1] - i1, srch_c[k1] - j1), axis=(0, 1)
                    )

                    Conv_CC_temp = np.sum(
                        (temp_result > 0) * temp_mask_Ipsi * (temp_mask_Ipsi > 0)
                    )
                    Conv_CI_temp = np.sum(
                        (temp_result < 0) * temp_mask_Ipsi * (temp_mask_Ipsi > 0)
                    )
                    Conv_SC_temp = np.sum(
                        (temp_result > 0) * temp_mask_Ipsi * (temp_mask_Ipsi < 0)
                    )
                    Conv_SI_temp = np.sum(
                        (temp_result < 0) * temp_mask_Ipsi * (temp_mask_Ipsi < 0)
                    )
                    Conv_temp[k1] = (Conv_CI_temp - Conv_CC_temp) - (
                        Conv_SI_temp - Conv_SC_temp
                    )
                    Aff_to_be_replaced[k1] = result_OD_sorted[srch_r[k1], srch_c[k1]]

            # check if any candidate move is better than the current coordinate
            if np.sum(Conv_temp > Conv_Result) > 0:
                # use best candidate
                ind_replace = np.argmax(Conv_temp)

                # swap
                result_OD_sorted[
                    srch_r[ind_replace], srch_c[ind_replace]
                ] = result_OD_sorted[i1, j1]
                result_OD_sorted[i1, j1] = Aff_to_be_replaced[ind_replace]

                # swap
                Crtx_plt_ONOFF_OD_temp = Crtx_plt_ONOFF_OD_sorted[i1, j1]
                Crtx_plt_ONOFF_OD_sorted[i1, j1] = Crtx_plt_ONOFF_OD_sorted[
                    srch_r[ind_replace], srch_c[ind_replace]
                ]
                Crtx_plt_ONOFF_OD_sorted[
                    srch_r[ind_replace], srch_c[ind_replace]
                ] = Crtx_plt_ONOFF_OD_temp

                # swap
                ret_temp = ret_sorted[srch_r[ind_replace], srch_c[ind_replace]]
                ret_sorted[srch_r[ind_replace], srch_c[ind_replace]] = ret_sorted[
                    i1, j1
                ]
                ret_sorted[i1, j1] = ret_temp

                # swap
                ret_temp2 = retinotopyPlot[srch_r[ind_replace], srch_c[ind_replace]]
                retinotopyPlot[
                    srch_r[ind_replace], srch_c[ind_replace]
                ] = retinotopyPlot[i1, j1]
                retinotopyPlot[i1, j1] = ret_temp2

        OD_holder.append(result_OD_sorted.copy())
        ONOFF_holder.append(Crtx_plt_ONOFF_OD_sorted.copy())
        retinotopy_holder.append(ret_sorted.copy())
        retinotopy_plot_holder.append(retinotopyPlot.copy())

    map_step_holder = MapStepHolder(
        OD_holder=OD_holder,
        ONOFF_holder=ONOFF_holder,
        retinotopy_holder=retinotopy_holder,
        retinotopy_plot_holder=retinotopy_plot_holder,
    )

    appdata.OD_ODSort = result_OD_sorted
    appdata.ONOFF_ODSort = Crtx_plt_ONOFF_OD_sorted
    appdata.RetODsorted = ret_sorted
    appdata.RetinotopyODSortedPlot = retinotopyPlot
    appdata.map_step_holder = map_step_holder
