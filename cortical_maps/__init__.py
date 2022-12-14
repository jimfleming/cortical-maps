import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
import skimage.morphology

from cortical_maps.make_colormap import make_colormap
from cortical_maps.make_retina_grid_rf_space_app import make_retina_grid_rf_space_app
from cortical_maps.make_cortex_app import make_cortex_app
from cortical_maps.functionDoGfilter import functionDoGfilter
from cortical_maps.sort_afferent_od2_app import sort_afferent_od2_app
from cortical_maps.sort_afferent_onoff2_app import sort_afferent_onoff2_app
from cortical_maps.smooth_map_app import smooth_map_app
from cortical_maps.make_rf_retina import make_rf_retina
from cortical_maps.generate_rf_cortex_primord_app_sn import (
    generate_rf_cortex_primord_app_sn,
)
from cortical_maps.label_ononff_od_island_gui import label_ononff_od_island_gui


def update_thalamic_aff_UI(appdata):
    fig, (ax, ax_density) = plt.subplots(1, 2, figsize=(6, 3))

    # update thalamic afferent sampling density
    CrtxONOFF3m = appdata.CrtxONOFF3m
    Retinotopy3mIndex = appdata.Retinotopy3mIndex
    rf_space_x = appdata.rf_space_x

    sum_aff = np.zeros_like(rf_space_x)
    disk_element = skimage.morphology.disk(round(appdata.RFONCenter))

    w_cortex = CrtxONOFF3m.shape[0]
    hh = round(w_cortex / 2)

    cortical_radius = 16
    cortical_radius_half = cortical_radius // 2
    cortical_area = cortical_radius**2

    row_aff_all = np.zeros(cortical_area)
    col_aff_all = np.zeros(cortical_area)
    cc = 0

    # sample 8 cortical pixels along each axis (64 afferents in total)
    for ii in range(hh - cortical_radius_half, hh + cortical_radius_half):
        for jj in range(hh - cortical_radius_half, hh + cortical_radius_half):
            indRf = Retinotopy3mIndex[ii, jj]
            row_aff, col_aff = np.nonzero(indRf == appdata.IndRetinaRfspace)

            if CrtxONOFF3m[ii, jj] > 0:
                circle = matplotlib.patches.Circle(
                    (col_aff[0], row_aff[0]),
                    appdata.RFONCenter,
                    color="r",
                    linewidth=0.5,
                    fill=False,
                )
                ax.add_patch(circle)
            elif CrtxONOFF3m[ii, jj] < 0:
                circle = matplotlib.patches.Circle(
                    (col_aff[0], row_aff[0]),
                    appdata.RFOFFCenter,
                    color="b",
                    linewidth=0.5,
                    fill=False,
                )
                ax.add_patch(circle)

            row_aff_all[cc] = row_aff
            col_aff_all[cc] = col_aff
            cc += 1

            # to calculate the number of overlapping receptive fields in visual space
            temp_sum_aff = np.zeros_like(rf_space_x)
            temp_sum_aff[row_aff, col_aff] = 1
            sum_aff += skimage.morphology.dilation(temp_sum_aff, disk_element)

    ax.set_xlim(col_aff_all.min() - 10, col_aff_all.max() + 10)
    ax.set_ylim(row_aff_all.min() - 10, row_aff_all.max() + 10)
    ax.set_aspect("equal")
    ax.axis("off")

    max_number_overlap_rf = sum_aff.max()
    ax.set_title(f"Afferent Sampling Density")

    cax_density = ax_density.imshow(sum_aff, cmap="turbo", interpolation="bicubic")
    ax_density.axis("off")
    fig.colorbar(cax_density, ax=ax_density)
    ax_density.set_title(f"{max_number_overlap_rf} RFs Per Point")

    plt.tight_layout()


def show_maps_segregation(appdata):
    ONOFF_smoothed = appdata.ONOFF_smoothed
    ODCrtxPlt_smooth = appdata.ODCrtxPlt_smooth
    map_retinotopy = appdata.RetinotopyONOFFSortedPlot

    _, (axes_od, axes_od_onoff, axes_retinotopy) = plt.subplots(1, 3, figsize=(9, 3))

    axes_od.imshow(ODCrtxPlt_smooth > 0.5)
    axes_od.contour(
        ODCrtxPlt_smooth, appdata.od_contour_levels, colors="k", linewidths=1
    )
    axes_od.axis("off")
    axes_od.set_title("Cortex OD")

    z3 = (((ODCrtxPlt_smooth > 0.5)).astype(float) + 1) * 1
    z4 = ((ONOFF_smoothed < 0.5).astype(float) + 0) * 2
    axes_od_onoff.imshow(z3 + z4)
    axes_od_onoff.contour(
        ODCrtxPlt_smooth, appdata.od_contour_levels, colors="k", linewidths=1
    )
    axes_od_onoff.axis("off")
    axes_od_onoff.set_title("Cortex ON/OFF")

    axes_retinotopy.imshow(map_retinotopy)
    axes_retinotopy.axis("off")
    axes_retinotopy.set_title("Cortex Retinotopy")

    plt.tight_layout()


def update_asf_figure(ax, appdata):
    CenterRadiusX = appdata.CenterRadiusX
    CenterRadiusY = appdata.CenterRadiusY
    SurroundRadiusX = appdata.SurroundRadiusX
    SurroundRadiusY = appdata.SurroundRadiusY

    CrtxLength = appdata.CrtxLength
    hh = round(CrtxLength / 2)  # half cortical length

    center = matplotlib.patches.Ellipse(
        (hh, hh),
        CenterRadiusX,
        CenterRadiusY,
        linewidth=2,
        fill=False,
    )
    surround = matplotlib.patches.Ellipse(
        (hh, hh),
        SurroundRadiusX,
        SurroundRadiusY,
        linewidth=2,
        fill=False,
    )

    ax.add_patch(center)
    ax.add_patch(surround)

    ax.set_title("Afferent Sorting Filter")
    ax.set_aspect("equal")
    ax.set_xlim(0, CrtxLength)
    ax.set_ylim(0, CrtxLength)
    ax.invert_yaxis()


def update_primord_orimap_UI(appdata):
    primord_orimap_contra = appdata.data_primord_contra.OriPreferred
    primord_orimap_ipsi = appdata.data_primord_ipsi.OriPreferred

    _, (axes1, axes2) = plt.subplots(1, 2, figsize=(6, 3))
    axes1.imshow(primord_orimap_contra)
    axes1.set_xlabel("contra")
    axes1.title.set_color(appdata.contra_font_color)

    axes2.imshow(primord_orimap_ipsi)
    axes2.set_xlabel("ipsi")
    axes2.title.set_color(appdata.ipsi_font_color)
    plt.tight_layout()
