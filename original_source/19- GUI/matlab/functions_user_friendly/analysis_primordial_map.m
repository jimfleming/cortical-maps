function [output_primord1,reference_ori_map] = analysis_primordial_map(RetinaRF,rowRange,colRange,ODCrtxPlt,ONOFFCrtxPlt,distAPON, ... 
        RetONOFFsorted,rAffSpread,IndRetinaRfspace,weight_map_island,ODCrtxPlt_smooth,sigma_LHI,pix2deg,sf_lSamp,debug,show_fig)

    output_primord1 = make_rf_cortex(RetinaRF,rowRange,colRange,ODCrtxPlt,ONOFFCrtxPlt,distAPON, ... 
        RetONOFFsorted,rAffSpread,IndRetinaRfspace,weight_map_island,ODCrtxPlt_smooth,pix2deg,sf_lSamp,debug,show_fig);

    % Development of Primordial Orientation Map (Kai Lun)
    primordial_ori_map = output_primord1.OriPreferred; 
    [reference_ori_map,~,~,~] = orientation_map_development2(primordial_ori_map, 10, debug);

    % Modify RF reference synaptic weight based on LHI 
%     LHIdevelopedmap = measure_LHI(reference_ori_map,sigma_LHI); 
%     output_primord2 = make_rf_cortex(RetinaRF,rowRange,colRange,ODCrtxPlt,ONOFFCrtxPlt,distAPON, ... 
%         RetONOFFsorted,rAffSpread,IndRetinaRfspace,LHIdevelopedmap,ODCrtxPlt_smooth,pix2deg,sf_lSamp,debug,show_fig);

end 
