import numpy as np
from scipy.stats import multivariate_normal
from cortical_maps.modify_synaptic_weight import modify_synaptic_weight
from cortical_maps.smooth_ori_euler import smooth_ori_euler
from cortical_maps.imdialte_grid import imdialte_grid

def generate_rf_cortex_primordial_initial_gui(pref_eye,RetinaRF,r_covered_aff,rowRange,colRange,ODCrtxPlt,ONOFFCrtxPlt,
    RetONOFFsorted,sdspread,i_aff_rf_space,synaptic_weight_factor,ODCrtxPlt_smooth,pix2deg,sf_lSamp,n_ori_smooth,sigma_LHI_2d,debug,show_fig):
    # The spread Gaussian function is applied based on the distance in Cortrex not RF space
    # Measure receptive field in cortex after afferent spread in cortical plate
    # Cortical RFs : Gaussian weighted sum of afferents
    #
    # the RF is measured both for one eye (ipsi or contra)
    #
    #
    # r_covered_aff   # Average radius of area covered by an afferent in visual space (diameter ~= 1 deg)
    #
    # Sohrab Note       :     the old name was "make_rf_cortex_same_eye2"

    # Initialization (pre assigning the matrices)
    allCXrfnorm{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = []
    allCXrfON{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = []
    allCXrfOFF{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = []

    rRFcenter = np.zeros_like(ONOFFCrtxPlt)
    cRFcenter = np.zeros_like(ONOFFCrtxPlt)

    OriPreferred = np.zeros_like(ONOFFCrtxPlt)
    dist_center_on_off = np.zeros_like(ONOFFCrtxPlt)
    SFPref = np.zeros_like(ONOFFCrtxPlt)

    max_on_response = np.zeros_like(ONOFFCrtxPlt)
    max_off_response = np.zeros_like(ONOFFCrtxPlt)
    max_on_response_norm = np.zeros_like(ONOFFCrtxPlt)
    max_off_response_norm = np.zeros_like(ONOFFCrtxPlt)

    CxOnOffBal = np.zeros_like(ONOFFCrtxPlt)
    SF50_pref = np.zeros_like(ONOFFCrtxPlt)
    LPI = np.zeros_like(ONOFFCrtxPlt)
    CV = np.zeros_like(ONOFFCrtxPlt)
    ODI = np.zeros_like(ONOFFCrtxPlt)

    OriHist{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = []
    sf_tuning_all{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = []
    sf_tuning_dog_all{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = []

    val_rf_dom_at_r_coverage = np.zeros_like(ONOFFCrtxPlt)
    rf_dom_pow2 = np.zeros_like(ONOFFCrtxPlt)

    TAA = initialize_TAA(ONOFFCrtxPlt, rowRange, colRange)

    # Main
    mask_covered_aff  = strel('disk', r_covered_aff, 0)
    mask_covered_aff2 = strel('disk', r_covered_aff-1 , 0)

    weight_arbor_spread = make_cortical_spread_dist(sdspread);
    weight_coef_RF = synaptic_weight_factor / np.amax(synaptic_weight_factor) # % 0 to 1 Cycle per degree

    if pref_eye =='contra':
        ODMap = ODCrtxPlt == 1
    elif pref_eye =='ipsi':
        ODMap = ODCrtxPlt == -1;

    col_grid, row_grid = np.meshgrid(1:ODMap.shape[0])
    width_rf_space = RetinaRF{1}.shape[0]
    col_grid_rf, row_grid_rf = np.meshgrid(1:width_rf_space)

    for ii in  rowRange:
        for jj in colRange:
            WTemp = np.zeros_like(ODCrtxPlt)
            WTemp[ii,jj] = 1;
            tmpweights = filter2(weight_arbor_spread,WTemp)

            mask_spread = imdialte_grid(row_grid, col_grid, ii, jj, sdspread)
            SelectionCondition = ODMap == 1 & mask_spread == 1
            inda, indb = np.where( SelectionCondition )

            ODI[ii,jj] = np.sum(SelectionCondition*tmpweights)

            # dominant RF
            aff_center_weight_plot = np.full((len(inda),4), np.nan)
            dominant_pol = ONOFFCrtxPlt[ii, jj]
            non_dominant_pol = -dominant_pol
            allCXrfTempDominant, aff_center_weight_plot, TAA = gen_rf(RetinaRF, RetONOFFsorted, ONOFFCrtxPlt, i_aff_rf_space, inda, indb, tmpweights, dominant_pol, 1, aff_center_weight_plot, TAA, ii, jj, debug)

            # this method worked better to find the cneter of RFs this center in mature_rf code is used to rotate this RF and find new weight for afferents
            rf_norm = allCXrfTempDominant / np.amax(np.abs(allCXrfTempDominant))

            #  center of rfs = max resp of dominant polarity ( reason : rotation of non_dom rf around the max response  )

            # finding the center of rfs by averaging all the pixels
            center_rf_row, center_rf_col = find_center_rf(rf_norm)
            if np.isnan(center_rf_row):
                raise

            rf_norm2, val_rf_dom_at_r_coverage[ii,jj], rf_dom_pow2[ii,jj] = add_power_rf_dominant2(rf_norm, row_grid_rf, col_grid_rf, r_covered_aff, center_rf_row, center_rf_col, .1)

            # non_dominant RF
            temp_mask1 = imdialte_grid(row_grid_rf, col_grid_rf, center_rf_row, center_rf_col, r_covered_aff)

            weight_other_pol = (1 - np.abs(rf_norm2)) * temp_mask1
            allCXrfTempNonDominant, aff_center_weight_plot, TAA = gen_rf(RetinaRF, RetONOFFsorted, ONOFFCrtxPlt, i_aff_rf_space, inda, indb, weight_other_pol, non_dominant_pol, 0, aff_center_weight_plot, TAA, ii, jj, debug)

            if dominant_pol == 1:
                allCXrfTempON = allCXrfTempDominant
                allCXrfTempOFF = allCXrfTempNonDominant
            elif dominant_pol == -1:
                allCXrfTempOFF = allCXrfTempDominant
                allCXrfTempON = allCXrfTempNonDominant

            row_center_on,col_center_on =  np.where(np.amax(allCXrfTempON.flat) == allCXrfTempON)
            row_center_off,col_center_off =  np.where(np.amax(np.abs(allCXrfTempOFF.flat)) == np.abs(allCXrfTempOFF))

            # in case there are two separate subregions
            try
                dist_pol = norm(row_center_on-row_center_off, col_center_on-col_center_off)
                dist_center_on_off[ii,jj] = dist_pol
            except Exception:
                dist_center_on_off[ii,jj] = np.nan

            # Modifying the weight and tuning curve measurement
            LWON, LWOFF = modify_synaptic_weight(allCXrfTempON, allCXrfTempOFF, weight_coef_RF[ii,jj], ONOFFCrtxPlt(ii,jj))

            CXrfON = allCXrfTempON * LWON
            CXrfOFF = allCXrfTempOFF * LWOFF
            CxOnOffBal[ii,jj] = ONOFF_balance(CXrfON,CXrfOFF)

            RFSpaceSim = CXrfON + CXrfOFF ;
            SingleRFNorm = RFSpaceSim / np.amax(np.abs(RFSpaceSim))

            max_on_response[ii, jj] = max(allCXrfTempON)
            max_off_response[ii, jj] = max(np.abs(allCXrfTempOFF))

            [Max_on, Max_off] = find_max_resp_norm(SingleRFNorm)
            max_on_response_norm[ii, jj] = Max_on
            max_off_response_norm[ii, jj] = Max_off

            ori_tuning_resp, sf_tuning_resp, sf_tuning_resp_dog, cv, ori_pref, sf_pref_deg, sf50_deg, lpi, angles_bin, sf_bin_deg = \
               fft_ori_tuning_uf(SingleRFNorm,sf_lSamp)

            # using the center of dominant polarity (in making mature RFs, it might not work because after rotation, there might not be enough aff available for the other polarity)
            rfcenONOFFrow = center_rf_row(1)
            rfcenONOFFcol = center_rf_col(1)

            OriPreferred[ii, jj] = ori_pref
            SFPref[ii, jj] = sf_pref_deg
            SF50_pref[ii, jj] = sf50_deg
            LPI[ii, jj] = lpi
            CV[ii, jj] = cv

            allCXrfnorm[ii, jj] = RFSpaceSim
            allCXrfON[ii, jj] = CXrfON
            allCXrfOFF[ii, jj] = CXrfOFF

            rRFcenter[ii, jj] = rfcenONOFFrow
            cRFcenter[ii, jj] = rfcenONOFFcol

            OriHist[ii, jj] = ori_tuning_resp
            sf_tuning_all[ii, jj] = sf_tuning_resp
            sf_tuning_dog_all[ii, jj] = sf_tuning_resp_dog

    # ori smooth
    ori_map_smooth, ori_map_interpolated = smooth_ori_euler(OriPreferred,ODCrtxPlt_smooth,n_ori_smooth,0,0)

    # LHI measurement
    LHI_map = measure_LHI(OriPreferred,sigma_LHI_2d)
    LHI_map_smooth = measure_LHI(ori_map_smooth,sigma_LHI_2d)

    output_rf.allCXrfON = allCXrfON
    output_rf.allCXrfOFF = allCXrfOFF

    output_data.rRFcenter = rRFcenter
    output_data.cRFcenter = cRFcenter
    output_data.CxOnOffBal = CxOnOffBal

    output_data.max_on_response = max_on_response
    output_data.max_off_response = max_off_response
    output_data.max_on_response_norm = max_on_response_norm
    output_data.max_off_response_norm = max_off_response_norm

    output_data.OriPreferred = OriPreferred
    output_data.ori_map_smooth = ori_map_smooth
    output_data.ori_map_interpolated = ori_map_interpolated
    output_data.SFPref = SFPref 
    output_data.SF50Pref = SF50_pref
    output_data.LPI = LPI 
    output_data.CV = CV
    output_data.LHI_map = LHI_map
    output_data.LHI_map_smooth = LHI_map_smooth

    ODI = ODI / np.amax(ODI)

    output_data.ODI = ODI
    output_data.OriHist = OriHist
    output_data.OriBin = angles_bin
    output_data.sf_bin_deg = sf_bin_deg
    output_data.SFHist = sf_tuning_all
    output_data.SFHist_dog = sf_tuning_dog_all
    output_data.dist_center_on_off = dist_center_on_off
    output_data.val_rf_dom_at_r_coverage = val_rf_dom_at_r_coverage
    output_data.rf_dom_pow2 = rf_dom_pow2
    output_data.TAA = TAA

    return output_rf, output_data


def make_cortical_spread_dist(sdspread):
    r2 = int(math.ceil(sdspread*3))
    mu = [0, 0] # Center of mask
    rx1 = r2-1
    x1 = np.arange(-rx1, rx1+1)
    x2 = np.arange(-rx1, rx1+1)
    X1, X2 = np.meshgrid(x1, -x2)
    covariance1 = [[sdspread, 0], [0, sdspread]]
    w = multivariate_normal.pdf(np.hstack((X1.flat[:, None], X2.flat[:, None])), mu, covariance1)
    w = w.reshape(len(x2), len(x1))
    w = w / np.amax(w)
    return w

def gen_rf(RetinaRF, RetONOFFsorted, ONOFFCrtxPlt, i_aff_rf_space, row_ind, col_ind, weights, pol, is_dominant, aff_center_weight_plot, TAA, ii, jj, debug):
    allCXrf = np.zeros_like(i_aff_rf_space)
    for pp in range(len(row_ind)):
        indRf = RetONOFFsorted[row_ind[pp], col_ind[pp]]
        addon = RetinaRF[indRf]

        # using onoff_smooth instead of onoff may result in some off aff in rf_on and vice versa so just to avoid this condition the following line is added
        sgn_rf_retina = np.sign(np.sum(addon))

        if ONOFFCrtxPlt[row_ind[pp], col_ind[pp]] == pol and sgn_rf_retina == pol:
            if is_dominant:
                tmp_w = weights[row_ind[pp], col_ind[pp]]
                allCXrf += tmp_w * addon
            else:
                rCirc, cCirc = np.where(indRf == i_aff_rf_space)
                tmp_w = weights[rCirc, cCirc]
                allCXrf += tmp_w * addon   # check!

            if debug == 1:
                rCirc, cCirc = np.where(indRf == i_aff_rf_space)
                aff_center_weight_plot[pp, 0] = cCirc
                aff_center_weight_plot[pp, 1] = rCirc
                aff_center_weight_plot[pp, 2] = tmp_w
                aff_center_weight_plot[pp, 3] = ONOFFCrtxPlt[row_ind[pp], col_ind[pp]]

        if ONOFFCrtxPlt[row_ind[pp], col_ind[pp]] == pol:
        # using onoff_smooth instead of onoff may result in some off aff in rf_on and vice versa so just to avoid this condition the following line is added
        # there are conditions that ONOFFCrtxPlt[row_ind[pp], col_ind[pp]
        # == pol is met but sgn_rf_retina == pol is not true, because
        # smoothing makes some cortical locations to switch the polarity
            if is_dominant:
                tmp_w = weights[row_ind[pp], col_ind[pp]]
            else:
                rCirc, cCirc = np.where(indRf == i_aff_rf_space)
                tmp_w = weights[rCirc, cCirc]
            TAA[row_ind[pp]][col_ind[pp]][ii, jj] = tmp_w
    return allCXrf, aff_center_weight_plot, TAA


def ONOFF_balance(CXRFON,CXRFOFF):
    tmpCXRFON_norm = np.amax(np.abs(CXRFON ))
    tmpCXRFOFF_norm = np.amax(np.abs(CXRFOFF))
    CxOnOffBal =  (tmpCXRFON_norm - tmpCXRFOFF_norm) /(tmpCXRFON_norm + tmpCXRFOFF_norm)
    return CxOnOffBal

def find_max_resp_norm(SingleRFNorm):
    # Flatten the array to a 1D array
    SingleRFNorm_1d = SingleRFNorm.flatten()

    # Find the maximum response in the on RF region
    ind_on_rf = SingleRFNorm_1d > 0
    max_on_response_norm = np.amax(SingleRFNorm_1d[ind_on_rf])
    if max_on_response_norm == 0:
        max_on_response_norm = 0

    # Find the maximum response in the off RF region
    ind_off_rf = SingleRFNorm_1d < 0
    max_off_response_norm = np.amax(np.abs(SingleRFNorm_1d[ind_off_rf]))
    if max_off_response_norm == 0:
        max_off_response_norm = 0

    return max_on_response_norm, max_off_response_norm

def add_power_rf_dominant2(rf_norm, row_grid_rf, col_grid_rf, r_covered_aff, center_rf_row, center_rf_col, thresh_rf):
    # The function calculates the average value of rf_norm at the r_covered_aff
    # it adds power function to the rf_norm if it is above the thresh_rf
    # to bring the average
    #
    # Inputs
    # thresh_rf = .1;   the desired value of the sum of RF at r afferent coeverage value below thresh_rf
    # mask_covered_aff  = strel('disk',r_covered_aff,0);
    # mask_covered_aff2 = strel('disk',r_covered_aff-1 ,0);
    #
    # Outputs
    # rf_norm : The new rf_nom, if the value is under the threshold, this does not change
    # val_rf_dom_at_r_coverage : the average value of rf_norm at the r_covered_aff
    # rf_dom_pow2              : the power used to calculate the new rf_norm

    mask_spread1 = imdialte_grid(row_grid_rf, col_grid_rf, center_rf_row, center_rf_col, r_covered_aff-1)
    mask_spread2 = imdialte_grid(row_grid_rf, col_grid_rf, center_rf_row, center_rf_col, r_covered_aff)
    mask_coverage_circle = mask_spread2 - mask_spread1
    rf_dom_at_r_coverage = (mask_coverage_circle * rf_norm)
    val_rf_dom_at_r_coverage = np.abs(np.sum(rf_dom_at_r_coverage) / np.sum(mask_coverage_circle))

    rf_dom_pow2 = 0
    if val_rf_dom_at_r_coverage > thresh_rf:
        rf_dom_pow2 = np.log(thresh_rf)/np.log(val_rf_dom_at_r_coverage)
        rf_norm2 = abs(rf_norm ** rf_dom_pow2)
        rf_dom_at_r_coverage = (mask_coverage_circle * rf_norm2)
        val_rf_dom_at_r_coverage = np.sum(rf_dom_at_r_coverage) / np.sum(mask_coverage_circle)
    else:
        rf_norm2 = rf_norm
    return rf_norm2,val_rf_dom_at_r_coverage,rf_dom_pow2

def find_center_rf(rf_norm):
    # Finding the center of RFs by averaging the rf value for all the pixels
    rf_abs = np.abs(rf_norm)
    rf_size = rf_abs.shape[0]
    XX, YY = np.meshgrid(np.arange(1, rf_size+1), np.arange(1, rf_size+1))

    # Measure rf center by averaging the rf_abs
    sum_rf = np.sum(rf_abs)
    rf_xx = np.sum(rf_abs * XX) / sum_rf
    rf_yy = np.sum(rf_abs * YY) / sum_rf

    rf_center_row = int(np.round(rf_yy))
    rf_center_col = int(np.round(rf_xx))

    return rf_center_row, rf_center_col

def initialize_TAA(ONOFFCrtxPlt, row_range, col_range):
    num_rows = ONOFFCrtxPlt.shape[0]
    num_cols = ONOFFCrtxPlt.shape[1]
    TAA = [None] * num_rows
    for i in row_range:
        TAA[i] = [None] * num_cols
        for j in col_range:
            TAA[i][j] = np.zeros((num_rows, num_cols))
    return TAA
