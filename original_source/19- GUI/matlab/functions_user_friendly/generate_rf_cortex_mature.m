function [RF_data2,OriSFdata2,RF_data1,OriSFdata1] = generate_rf_cortex_mature(eyepref,rf_contra_primord,data_contra_primord,rf_ipsi_primord,data_ipsi_primord,reference_ori_map,RetinaRF,rowRange,colRange,...
    IndRetinaRfspace,rAffSpread,RetONOFFsorted,ONOFFCrtxPlt,ODCrtxPlt,ODCrtxPlt_smooth,electrode_position,sigma_LHI,n_ori_smooth,sf_lSamp,pix2deg,weight_coef_RF,alpha,ONOFFODLabelSorted,ONOFF_smoothed,r_covered_aff,debug)
%  Genrating mature cortical RFs 
%
%  RF_data1,OriSFdata1  :  they were saved to debug 
%  Note : the previous name of this function was generate_rf_eye1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
end 



