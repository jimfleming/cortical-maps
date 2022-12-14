import numpy as np
import skimage.measure
import skimage.morphology


def label_ononff_od_island_gui(appdata):
    ODCrtxPlt = appdata.ODCrtxPlt
    ONOFFCrtxPlt = appdata.ONOFFCrtxPlt

    # Labeling each ONOFF OD feature separately
    # in pixel (islands with area smaller than this number are removed)
    FeatureRemoveCriteria = 8
    disk_openning_size = 0

    ONOFFODLabel = label_onoff_unsorted(ODCrtxPlt, ONOFFCrtxPlt, disk_openning_size)
    FeatureCenter, FeatureArea, DistToCenter, CenterPWImage = calculate_center_island(
        ONOFFODLabel, FeatureRemoveCriteria
    )

    # Resort the pinwheels index based on the distance to the center of image
    # in pixel, removing small islands
    FeatureIncluded = FeatureArea > FeatureRemoveCriteria
    NumPinwheel = np.sum(FeatureIncluded)
    indPinwheelSorted = np.argsort(DistToCenter)
    ONOFFODLabelSorted = np.zeros_like(ONOFFODLabel)

    if appdata.aff_sampling_density >= 1:
        CenterPinwheelSorted = []
        cc = 0
        for bb in range(len(indPinwheelSorted)):  #  NumPinwheel
            iPw = indPinwheelSorted[bb]
            if FeatureArea(iPw) > FeatureRemoveCriteria:
                ONOFFODLabelSorted = ONOFFODLabelSorted + (ONOFFODLabel == iPw) * cc
                CenterPinwheelSorted.append(FeatureCenter[iPw])
                cc = cc + 1
    else:
        CenterPinwheelSorted = []
        CenterPWImage = []
    return ONOFFODLabelSorted, CenterPinwheelSorted, NumPinwheel, CenterPWImage


def label_onoff_unsorted(ODCrtxPlt, ONOFFCrtxPlt, disk_opening_size):
    SEOpenning = skimage.morphology.disk(disk_opening_size)
    BWNumConnection = 4

    ONOFFODTemp = ONOFFCrtxPlt * 2 + ODCrtxPlt
    IpsiOFFOpened = skimage.morphology.binary_opening(ONOFFODTemp == -3, SEOpenning)
    IpsiOFFLabeled, IpsiOFFNum = skimage.measure.label(
        IpsiOFFOpened, BWNumConnection, return_num=True
    )
    ContraOFFOpened = skimage.morphology.binary_opening(ONOFFODTemp == -1, SEOpenning)
    ContraOFFLabeled, ContraOFFNum = skimage.measure.label(
        ContraOFFOpened, BWNumConnection, return_num=True
    )
    ContraONOpened = skimage.morphology.binary_opening(ONOFFODTemp == 3, SEOpenning)
    ContraONLabeled, ContraONNum = skimage.measure.label(
        ContraONOpened, BWNumConnection, return_num=True
    )
    IpsiONOpened = skimage.morphology.binary_opening(ONOFFODTemp == 1, SEOpenning)
    IpsiONLabeled, _ = skimage.measure.label(
        IpsiONOpened, BWNumConnection, return_num=True
    )

    # Labelling again giving each ON OFF Ipsi Contra a number
    ONOFFODLabel1 = IpsiOFFLabeled
    ONOFFODLabel2 = ContraOFFLabeled + (ContraOFFLabeled != 0) * (IpsiOFFNum)
    ONOFFODLabel3 = ContraONLabeled + (ContraONLabeled != 0) * (
        IpsiOFFNum + ContraOFFNum
    )
    ONOFFODLabel4 = IpsiONLabeled + (IpsiONLabeled != 0) * (
        IpsiOFFNum + ContraOFFNum + ContraONNum
    )
    ONOFFODLabel = ONOFFODLabel1 + ONOFFODLabel2 + ONOFFODLabel3 + ONOFFODLabel4
    return ONOFFODLabel


def calculate_center_island(ONOFFODLabel, RemoveCriteria):
    # Finding the Center of each ONOFF OD feature
    # Calculating distance from center of Islands (for Spatial Frequency measurement)

    ONOFFODNum = np.amax(ONOFFODLabel)
    FeatureCenter = []
    FeatureArea = np.zeros(ONOFFODNum)
    DistToCenter = np.zeros(ONOFFODNum)
    RowHalfSize = ONOFFODLabel.shape[0] / 2
    ColHalfSize = ONOFFODLabel.shape[1] / 2
    CenterPWImage = np.zeros_like(ONOFFODLabel)
    for ii in range(ONOFFODNum):
        Feature = (ONOFFODLabel == (ii + 1)).astype(np.uint8)
        Center = skimage.measure.regionprops(Feature)[0].centroid

        # The result of regionprops could be out of the region so the closest point to the central line (skeletion) is found
        RowCenterMidLine, ColCenterMidLine = Center
        FeatureSkel = skimage.morphology.skeletonize(Feature)
        r1, c1 = np.nonzero(FeatureSkel)
        dist = (c1 - ColCenterMidLine) ** 2 + (r1 - RowCenterMidLine) ** 2
        aa = np.argmin(dist)
        FeatureCenter.append([r1[aa], c1[aa]])
        FeatureArea[ii] = np.sum(Feature)
        DistToCenter[ii] = np.sqrt(
            (c1[aa] - ColHalfSize) ** 2 + (r1[aa] - RowHalfSize) ** 2
        )

        if FeatureArea[ii] > RemoveCriteria:
            CenterPWImage[r1[aa], c1[aa]] = 1

    return FeatureCenter, FeatureArea, DistToCenter, CenterPWImage
