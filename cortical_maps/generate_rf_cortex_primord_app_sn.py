# from cortical_maps.generate_rf_cortex_primordial_initial_gui import generate_rf_cortex_primordial_initial_gui
# from cortical_maps.smooth_ori_euler import smooth_ori_euler


def generate_rf_cortex_primord_app_sn(appdata):
    RetinaRF = appdata.RetinaRF
    r_covered_aff = appdata.r_covered_aff
    acf_parameters_coverage = appdata.acf_parameters_coverage
    acf_parameters_smooth = appdata.acf_parameters_smooth
    rowRange = appdata.rowRange
    colRange = appdata.colRange
    ODCrtxPlt = appdata.ODCrtxPlt
    ONOFFCrtxPlt = appdata.ONOFFCrtxPlt
    RetONOFFsorted = appdata.RetONOFFsorted
    rAffSpread = appdata.rAffSpread
    IndRetinaRfspace = appdata.IndRetinaRfspace
    weight_map_island = appdata.weight_map_island
    ODCrtxPlt_smooth = appdata.ODCrtxPlt_smooth
    ONOFF_smoothed = appdata.ONOFF_smoothed
    pix2deg = appdata.pix2deg
    sf_lSamp = appdata.sf_lSamp
    sigma_LHI_2d = appdata.sigma_LHI_2d

    ONOFFODLabelSorted = appdata.ONOFFODLabelSorted
    NumPinwheel = appdata.NumPinwheel
    CenterPinwheelSorted = appdata.CenterPinwheelSorted
    n_ori_smooth = appdata.n_ori_smooth

    # using onoff_smooth instead of onoff may result in some off aff in rf_on and vice versa, because their polarity will change by smoothing at some location but the selected RGC does not change
    ONOFF_smoothed_binary = (ONOFF_smoothed > 0.5) * 1 + (ONOFF_smoothed <= 0.5) * -1

    #  make primordial cortical receptive fields

    rf_primord_contra, data_primord_contra = generate_rf_cortex_primordial_initial_gui(
        "contra",
        RetinaRF,
        r_covered_aff,
        rowRange,
        colRange,
        ODCrtxPlt,
        ONOFF_smoothed_binary,
        RetONOFFsorted,
        rAffSpread,
        IndRetinaRfspace,
        weight_map_island,
        ODCrtxPlt_smooth,
        pix2deg,
        sf_lSamp,
        n_ori_smooth,
        sigma_LHI_2d,
    )
    rf_primord_ipsi, data_primord_ipsi = generate_rf_cortex_primordial_initial_gui(
        "ipsi",
        RetinaRF,
        r_covered_aff,
        rowRange,
        colRange,
        ODCrtxPlt,
        ONOFF_smoothed_binary,
        RetONOFFsorted,
        rAffSpread,
        IndRetinaRfspace,
        weight_map_island,
        ODCrtxPlt_smooth,
        pix2deg,
        sf_lSamp,
        n_ori_smooth,
        sigma_LHI_2d,
    )

    # to make the reference orientation map, the orientation value at the dominant eye is considered
    primord_orimap_contra = data_primord_contra.OriPreferred
    primord_orimap_contra_smooth = data_primord_contra.ori_map_smooth
    primord_orimap_ipsi = data_primord_ipsi.OriPreferred
    primord_orimap_ipsi_smooth = data_primord_ipsi.ori_map_smooth

    primord_orimap_dominant_eye = primord_orimap_contra * (
        ODCrtxPlt == 1
    ) + primord_orimap_ipsi * (ODCrtxPlt == -1)
    primord_orimap_smooth, _ = smooth_ori_euler(
        primord_orimap_dominant_eye, [], n_ori_smooth, 0, 0
    )
    show_orimap_contra_ipsi(
        primord_orimap_contra_smooth,
        primord_orimap_ipsi_smooth,
        primord_orimap_smooth,
        ODCrtxPlt_smooth,
    )

    # Modify coverage of ori map for each island
    ori_result_coverage = modify_ori_coverage2(
        primord_orimap_smooth,
        ONOFFODLabelSorted,
        NumPinwheel,
        CenterPinwheelSorted,
        acf_parameters_coverage,
    )  # maximum coverage within an island
    ori_map_coverage_smooth1, _ = smooth_ori_euler(
        ori_result_coverage, [], n_ori_smooth, 0, 0
    )
    reference_ori_map, _, _, _ = orientation_map_development2(
        ori_map_coverage_smooth1, acf_parameters_smooth, 10, 0
    )  # swindale

    appdata.rf_primord_contra = rf_primord_contra
    appdata.data_primord_contra = data_primord_contra
    appdata.rf_primord_ipsi = rf_primord_ipsi
    appdata.data_primord_ipsi = data_primord_ipsi
    appdata.reference_ori_map = reference_ori_map
    appdata.ori_map_coverage_smooth1 = ori_map_coverage_smooth1

    appdata.primord_orimap_smooth = primord_orimap_smooth
