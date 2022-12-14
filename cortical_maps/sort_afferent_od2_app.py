import numpy as np
import scipy.ndimage
from dataclasses import dataclass
from tqdm import trange


@dataclass
class MapStepHolder:
    OD_holder: np.ndarray
    ONOFF_holder: np.ndarray
    retinotopy_holder: np.ndarray


def sort_afferent_od2_app(appdata):
    filt_sort = appdata.SortFiltOD
    Crtx_plt_OD = appdata.CrtxOD3m
    Crtx_plt_ONOFF = appdata.CrtxONOFF3m
    ret_initial = appdata.Retinotopy3mIndex
    retinotopyPlot = appdata.Retinotopy3mIndPlot
    N_repeat = appdata.NSortOD
    rng_trial = appdata.rng_trial

    rng = np.random.default_rng(seed=rng_trial)

    result_OD_sorted = Crtx_plt_OD
    Crtx_plt_ONOFF_OD_sorted = Crtx_plt_ONOFF
    ret_sorted = ret_initial

    OD_holder = []
    ONOFF_holder = []
    retinotopy_holder = []
    retinotopy_plot_holder = []
    for n in trange(N_repeat):
        #   Selection of points should be randomized
        Affr, Affc = np.nonzero(result_OD_sorted != 0)
        ind_rand2 = rng.permutation(len(Affr))
        for q in range(len(Affr)):
            ind_replace = []
            i1 = Affr[ind_rand2[q]]
            j1 = Affc[ind_rand2[q]]
            temp = np.zeros_like(result_OD_sorted)
            temp[i1, j1] = 1
            temp_mask1 = scipy.ndimage.convolve(temp, filt_sort, mode="constant")

            # Convolution Center Contra
            Conv_CC = (result_OD_sorted > 0) * temp_mask1 * (temp_mask1 > 0)
            Sum_CC = np.sum(np.abs(Conv_CC))

            # Convolution Center Ipsi
            Conv_CI = (result_OD_sorted < 0) * temp_mask1 * (temp_mask1 > 0)
            Sum_CI = np.sum(np.abs(Conv_CI))

            # Convolution Surround Contra
            Conv_SC = (result_OD_sorted > 0) * temp_mask1 * (temp_mask1 < 0)
            Sum_SC = np.sum(np.abs(Conv_SC))

            # Convolution Surround Ipsi
            Conv_SI = (result_OD_sorted < 0) * temp_mask1 * (temp_mask1 < 0)
            Sum_SI = np.sum(np.abs(Conv_SI))

            if result_OD_sorted[i1, j1] == 1:
                Conv_Result = (Sum_CC - Sum_CI) - (Sum_SC - Sum_SI)
                mask_search = (
                    scipy.ndimage.convolve(temp, np.ones((3, 3)), mode="constant")
                    - temp
                )
                srch_r, srch_c = np.nonzero(
                    mask_search.astype(bool) & (result_OD_sorted != 1)
                )

                Conv_temp = np.zeros(len(srch_r))
                Aff_to_be_replaced = np.zeros(len(srch_r))

                for k1 in range(len(srch_r)):
                    temp_result = result_OD_sorted.copy()
                    temp_result[i1, j1] = result_OD_sorted[srch_r[k1], srch_c[k1]]
                    temp_result[srch_r[k1], srch_c[k1]] = 1
                    temp_mask_Contra = np.roll(
                        temp_mask1, ([srch_r[k1] - i1, srch_c[k1] - j1])
                    )

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
            elif result_OD_sorted[i1, j1] == -1:
                Conv_Result = (Sum_CI - Sum_CC) - (Sum_SI - Sum_SC)

                mask_search = (
                    scipy.ndimage.convolve(temp, np.ones((3, 3)), mode="constant")
                    - temp
                )
                srch_r, srch_c = np.nonzero(
                    mask_search.astype(bool) & (result_OD_sorted != -1)
                )

                Conv_temp = np.zeros(len(srch_r))
                Aff_to_be_replaced = np.zeros(len(srch_r))

                for k1 in range(len(srch_r)):
                    temp_result = result_OD_sorted.copy()
                    temp_result[i1, j1] = result_OD_sorted[srch_r[k1], srch_c[k1]]
                    temp_result[srch_r[k1], srch_c[k1]] = 1
                    temp_mask_Contra = np.roll(
                        temp_mask1, ([srch_r[k1] - i1, srch_c[k1] - j1])
                    )

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
                    Conv_temp[k1] = (Conv_CI_temp - Conv_CC_temp) - (
                        Conv_SI_temp - Conv_SC_temp
                    )

                    Aff_to_be_replaced[k1] = result_OD_sorted[srch_r[k1], srch_c[k1]]

            if np.sum(Conv_temp > Conv_Result):
                (ind_replace_vec,) = np.nonzero(Conv_temp == Conv_temp.max())
                ind_replace = ind_replace_vec[0]
                result_OD_sorted[
                    srch_r[ind_replace], srch_c[ind_replace]
                ] = result_OD_sorted[i1, j1]
                result_OD_sorted[i1, j1] = Aff_to_be_replaced[ind_replace]

                Crtx_plt_ONOFF_OD_temp = Crtx_plt_ONOFF_OD_sorted[i1, j1]
                Crtx_plt_ONOFF_OD_sorted[i1, j1] = Crtx_plt_ONOFF_OD_sorted[
                    srch_r[ind_replace], srch_c[ind_replace]
                ]
                Crtx_plt_ONOFF_OD_sorted[
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

        OD_holder.append(result_OD_sorted)
        ONOFF_holder.append(Crtx_plt_ONOFF_OD_sorted)
        retinotopy_holder.append(ret_sorted)
        retinotopy_plot_holder.append(retinotopyPlot)

    map_step_holder = MapStepHolder(
        OD_holder=OD_holder,
        ONOFF_holder=ONOFF_holder,
        retinotopy_holder=retinotopy_plot_holder,
    )

    appdata.OD_ODSort = result_OD_sorted
    appdata.ONOFF_ODSort = Crtx_plt_ONOFF_OD_sorted
    appdata.RetODsorted = ret_sorted
    appdata.RetinotopyODSortedPlot = retinotopyPlot
    appdata.map_step_holder = map_step_holder
