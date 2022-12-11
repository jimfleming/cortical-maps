function [rf_primord_contra,data_primord_contra,rf_primord_ipsi,data_primord_ipsi,reference_ori_map,ori_map_coverage_smooth1] = generate_rf_cortex_primord(RetinaRF,r_covered_aff,acf_parameters_coverage,acf_parameters_smooth,rowRange,colRange,ODCrtxPlt,ONOFFCrtxPlt, ...
    RetONOFFsorted,rAffSpread,IndRetinaRfspace,weight_map_island,ODCrtxPlt_smooth,ONOFF_smoothed,pix2deg,sf_lSamp,sigma_LHI_2d,debug,ONOFFODLabelSorted,NumPinwheel,CenterPinwheelSorted,n_ori_smooth,show_fig,save_file)

    % using onoff_smooth instead of onoff may result in some off aff in rf_on and vice versa, because their polarity will change by smoothing at some location but the selected RGC does not change 
    ONOFF_smoothed_binary = (ONOFF_smoothed > .5) * 1 + (ONOFF_smoothed <= .5) * -1; 

    %%  make primordial cortical receptive fields 
%     [rf_primord_contra,data_primord_contra] = generate_rf_cortex_primordial_initial3('contra',RetinaRF,r_covered_aff,rowRange,colRange,ODCrtxPlt,ONOFF_smoothed_binary, ... 
%     RetONOFFsorted,rAffSpread,IndRetinaRfspace,weight_map_island,ODCrtxPlt_smooth,pix2deg,sf_lSamp,n_ori_smooth,sigma_LHI_2d,debug,show_fig);

    [rf_primord_contra,data_primord_contra] = generate_rf_cortex_primordial_initial3_2_uf('contra',RetinaRF,r_covered_aff,rowRange,colRange,ODCrtxPlt,ONOFF_smoothed_binary, ... 
    RetONOFFsorted,rAffSpread,IndRetinaRfspace,weight_map_island,ODCrtxPlt_smooth,pix2deg,sf_lSamp,n_ori_smooth,sigma_LHI_2d,debug,show_fig);

%     [rf_primord_ipsi,data_primord_ipsi] = generate_rf_cortex_primordial_initial3('ipsi',RetinaRF,r_covered_aff,rowRange,colRange,ODCrtxPlt,ONOFF_smoothed_binary, ... 
%     RetONOFFsorted,rAffSpread,IndRetinaRfspace,weight_map_island,ODCrtxPlt_smooth,pix2deg,sf_lSamp,n_ori_smooth,sigma_LHI_2d,debug,show_fig);

    [rf_primord_ipsi,data_primord_ipsi] = generate_rf_cortex_primordial_initial3_2_uf('ipsi',RetinaRF,r_covered_aff,rowRange,colRange,ODCrtxPlt,ONOFF_smoothed_binary, ... 
    RetONOFFsorted,rAffSpread,IndRetinaRfspace,weight_map_island,ODCrtxPlt_smooth,pix2deg,sf_lSamp,n_ori_smooth,sigma_LHI_2d,debug,show_fig);

    if save_file == 1 
        save('rf_primord_contra.mat','rf_primord_contra')
        save('data_primord_contra.mat','data_primord_contra')
        save('rf_primord_ipsi.mat','rf_primord_ipsi')
        save('data_primord_ipsi.mat','data_primord_ipsi')
    end 
    %% to make the reference orientation map, the orientation value at the dominant eye is considered
    primord_orimap_contra = data_primord_contra.OriPreferred;
    primord_orimap_contra_smooth = data_primord_contra.ori_map_smooth;
    primord_orimap_ipsi = data_primord_ipsi.OriPreferred;
    primord_orimap_ipsi_smooth = data_primord_ipsi.ori_map_smooth; 

    primord_orimap_dominant_eye = primord_orimap_contra .* (ODCrtxPlt == 1) + primord_orimap_ipsi .* (ODCrtxPlt == -1);
    [primord_orimap_smooth,~] = smooth_ori_euler(primord_orimap_dominant_eye,[],n_ori_smooth,0,0);
    show_orimap_contra_ipsi(primord_orimap_contra_smooth,primord_orimap_ipsi_smooth,primord_orimap_smooth,ODCrtxPlt_smooth,show_fig)
 
    %% Modify coverage of ori map for each island 
    [ori_result_coverage] = modify_ori_coverage2(primord_orimap_smooth,ONOFFODLabelSorted,NumPinwheel,CenterPinwheelSorted,acf_parameters_coverage); % maximum coverage within an island
    [ori_map_coverage_smooth1,~] = smooth_ori_euler(ori_result_coverage,[],n_ori_smooth,0,0);
    [reference_ori_map,~,~,~] = orientation_map_development2(ori_map_coverage_smooth1,acf_parameters_smooth, 10, 0); % swindale 
    show_ori_comparison(primord_orimap_smooth,ori_map_coverage_smooth1,reference_ori_map,ODCrtxPlt_smooth,show_fig)
  
end 

%% functions plots 
function show_onoff_balance(output_primord1_contra,output_primord1_ipsi,ODCrtxPlt_smooth,show_fig)
    if show_fig == 1   
        figure
        set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
        subplot(121),imagesc(output_primord1_contra.CxOnOffBal),colormap(gca,'jet'),axis('square'),colorbar,title('ONOFF balance Contra','fontsize',20),axis off,caxis([-.3 .3])
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5); caxis(.5*[-1 1])
        %   xlim([rect_x1 rect_x2]), ylim([rect_y1 rect_y2])
        %   hold on, rectangle('Position',[rect_x1 rect_y1 length_crtx_present length_crtx_present],'LineWidth',3)
        subplot(122), imagesc(output_primord1_ipsi.CxOnOffBal), colormap(gca,'jet'), axis('square'), colorbar, title('ONOFF balance Ipsi','fontsize',20),axis off,caxis([-.3 .3])
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5); caxis(.5*[-1 1])
        %   xlim([rect_x1 rect_x2]), ylim([rect_y1 rect_y2])
        %   hold on, rectangle('Position',[rect_x1 rect_y1 length_crtx_present length_crtx_present],'LineWidth',3)
    end 
end 

function show_orimap_contra_ipsi(primordial_ori_contra_smooth,primordial_ori_ipsi_smooth,primordial_ori_map_smooth,ODCrtxPlt_smooth,show_fig)
    if show_fig == 1 
        figure
        subplot(131),imagesc(primordial_ori_contra_smooth),colormap('hsv'), freezeColors , axis('square'),colorbar, title('Contra','fontsize',20), axis off
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5);
        subplot(132),imagesc(primordial_ori_ipsi_smooth),colormap('hsv'),freezeColors , axis('square'),colorbar, title('Ipsi','fontsize',20), axis off
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5);
        subplot(133),imagesc(primordial_ori_map_smooth),colormap('hsv'),freezeColors , axis('square'),colorbar, title('Primordial map','fontsize',20), axis off
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5);
    end 
end 


function show_ori_comparison(primordial_ori_map_smooth,ori_map_coverage_smooth1,reference_ori_map,ODCrtxPlt_smooth,show_fig)
    if show_fig == 1 
        figure
        subplot(131),imagesc(primordial_ori_map_smooth),colormap('hsv'), freezeColors , axis('square'),colorbar, title('Primordial map','fontsize',20), axis off
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5);
        subplot(132),imagesc(ori_map_coverage_smooth1),colormap('hsv'),freezeColors , axis('square'),colorbar, title('Coverage Modified','fontsize',20), axis off
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5);
        subplot(133),imagesc(reference_ori_map),colormap('hsv'),freezeColors , axis('square'),colorbar, title('Swindale','fontsize',20), axis off
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5);
    end 
end 


