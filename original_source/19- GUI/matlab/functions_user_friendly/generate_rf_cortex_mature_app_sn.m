function [appdata] = generate_rf_cortex_mature_app_sn(appdata)
%  Genrating mature cortical RFs 
%
%  RF_data1,OriSFdata1  :  they were saved to debug 
%  Note : the previous name of this function was generate_rf_eye1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eyepref = appdata.eyepref;
rf_contra_primord = appdata.rf_primord_contra;
data_contra_primord = appdata.data_primord_contra;
rf_ipsi_primord = appdata.rf_primord_ipsi;
data_ipsi_primord = appdata.data_primord_ipsi;
reference_ori_map = appdata.reference_ori_map;
RetinaRF = appdata.RetinaRF;
rowRange = appdata.rowRange;
colRange = appdata.colRange;
IndRetinaRfspace = appdata.IndRetinaRfspace;
rAffSpread = appdata.rAffSpread;
RetONOFFsorted = appdata.RetONOFFsorted;
ONOFFCrtxPlt = appdata.ONOFFCrtxPlt;
ODCrtxPlt = appdata.ODCrtxPlt;
ODCrtxPlt_smooth = appdata.ODCrtxPlt_smooth;
electrode_position = appdata.electrode_position;
sigma_LHI = appdata.sigma_LHI_2d;
n_ori_smooth = appdata.n_ori_smooth;
sf_lSamp = appdata.sf_lSamp;
pix2deg = appdata.pix2deg;
weight_coef_RF = appdata.weight_map_rf_reference_orimap;
alpha = appdata.alpha;
ONOFFODLabelSorted = appdata.ONOFFODLabelSorted;
ONOFF_smoothed = appdata.ONOFF_smoothed;
r_covered_aff = appdata.r_covered_aff;
debug = 0;
%%
    % initial rf mature generation 
    [RF_data1,OriSFdata1] = generate_rf_cortex_mature_initial2_uf(eyepref,rf_contra_primord,data_contra_primord,rf_ipsi_primord,data_ipsi_primord,reference_ori_map,RetinaRF,rowRange,colRange,...
        IndRetinaRfspace,rAffSpread,RetONOFFsorted,ONOFFCrtxPlt,ODCrtxPlt,ODCrtxPlt_smooth,electrode_position,sigma_LHI,n_ori_smooth,sf_lSamp,pix2deg,weight_coef_RF,alpha,debug);

    % mature rf generation (adjusting rf weights with pw locations)
    ori_map_eye = OriSFdata1.ori_map_smooth; 
    [weight_map2_pw] = synaptic_weight_island2(ori_map_eye,ONOFFODLabelSorted,ODCrtxPlt_smooth,ONOFF_smoothed,0);
    weight_map2_pw = (weight_map2_pw * 1.5);  % 1.5 resuls in less lpi = 1 
%     [RF_data2,OriSFdata2] = generate_rf_cortex_mature_final_weight(eyepref,RF_data1,OriSFdata1,reference_ori_map,RetinaRF,rowRange,colRange,...
%         IndRetinaRfspace,rAffSpread,RetONOFFsorted,ONOFFCrtxPlt,ODCrtxPlt,ODCrtxPlt_smooth,electrode_position,sigma_LHI,n_ori_smooth,sf_lSamp,pix2deg,weight_map2_pw,alpha,debug);
    [RF_data2,OriSFdata2] = generate_rf_cortex_mature_final_weight2(eyepref,RF_data1,OriSFdata1,reference_ori_map,RetinaRF,rowRange,colRange,...
        IndRetinaRfspace,rAffSpread,RetONOFFsorted,ONOFFCrtxPlt,ODCrtxPlt,ODCrtxPlt_smooth,electrode_position,sigma_LHI,n_ori_smooth,sf_lSamp,pix2deg,weight_map2_pw,alpha,debug);
    
    
    OriSFdata2.weight_map2_pw = weight_map2_pw;
    OriSFdata2.ODI = data_contra_primord.ODI; 
    OriSFdata2.num_projected_aff = OriSFdata1.num_projected_aff; 
    OriSFdata2.R_iso_retinotopy = OriSFdata1.R_iso_retinotopy; 
    OriSFdata2.rRFcenter = OriSFdata1.rRFcenter;
    OriSFdata2.cRFcenter = OriSFdata1.cRFcenter;
    
    if strcmp(eyepref,'contra')
        appdata.rf_contra_mature = RF_data2;
        appdata.data_contra_mature = OriSFdata2;
    elseif strcmp(eyepref,'ipsi') 
        appdata.rf_ipsi_mature = RF_data2;
        appdata.data_ipsi_mature = OriSFdata2;
    end
end 



