import numpy as np


def make_retina_grid_rf_space_app(appdata):
    rf_space_x = appdata.rf_space_x
    rf_space_y = appdata.rf_space_y

    TotalRetinaCells = appdata.TotalRetinaCells
    RetinaAllX = appdata.RetinaAllX
    RetinaAllY = appdata.RetinaAllY
    RetinaAllOD = appdata.RetinaAllOD
    RetinaAllONOFF = appdata.RetinaAllONOFF
    rng_trial = appdata.rng_trial

    # Putting the cells in Matrix  without changing the location (cortical plate)
    # plotting the RGCs on Retinotopic map
    rng = np.random.default_rng(seed=rng_trial)

    TakenSpots = np.zeros_like(rf_space_x)
    od_rf_space = np.zeros_like(rf_space_x)
    onoff_rf_space = np.zeros_like(rf_space_x)
    IndRetinaRfspace = np.full_like(rf_space_x, -1)

    indSelectRand = rng.permutation(TotalRetinaCells)
    for ii in range(TotalRetinaCells):
        indTmp = indSelectRand[ii]
        distLocRetCrtx = (rf_space_x.flat - RetinaAllX[indTmp]) ** 2 + (
            rf_space_y.flat - RetinaAllY[indTmp]
        ) ** 2
        indSelectLoc = np.argsort(distLocRetCrtx)
        indSelectLoc = indSelectLoc[TakenSpots.flat[indSelectLoc] != 1]
        od_rf_space.flat[indSelectLoc[0]] = RetinaAllOD[indTmp]
        onoff_rf_space.flat[indSelectLoc[0]] = RetinaAllONOFF[indTmp]
        IndRetinaRfspace.flat[indSelectLoc[0]] = indTmp
        TakenSpots.flat[indSelectLoc[0]] = 1

    appdata.od_rf_space = od_rf_space
    appdata.onoff_rf_space = onoff_rf_space
    appdata.IndRetinaRfspace = IndRetinaRfspace
