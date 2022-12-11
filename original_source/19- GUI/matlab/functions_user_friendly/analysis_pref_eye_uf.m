function analysis_pref_eye_uf(eyepref,OriSFdata,cortex_RF,ODCrtxPlt_smooth,ODCrtxPlt_interpolated,electrode_position,boundary,NumRandPoint,show_fig)

%   rf_analysis, save_address, save_data, save_RF
[rect1,rect2,size_rect1,size_rect2] = remove_boundary(OriSFdata,boundary); 
% selecting points randomly
rand_points =   randperm(size_rect1, NumRandPoint);
rand_points_interpol = randperm(size_rect2, NumRandPoint);

%% figures  
% annotation('textbox',[0.05 .8 .1 .1],'String',{eyepref},'FitBoxToText','on','FontSize',24);

if show_fig     
    %   show_ori_map(ori_map,ori_map_interpolated,ODCrtxPlt_smooth)
    %   show_retinotopy_map(OriSFdata,ODCrtxPlt_smooth)
    %   show_max_onoff_response(max_on_response,max_off_response,electrode_position)
    %   show_num_projected_afferents(OriSFdata)
    %   [distSameOD,distOPPOD] = RetinotopyOdSkeletDist(ODCrtxPlt,RetONOFFsorted,IndRetinaRfspace,CenterRadiusX)
    
    % show maps 
    show_raw_map(OriSFdata,ODCrtxPlt_smooth,eyepref)
    show_interpolated_map(OriSFdata,ODCrtxPlt_interpolated,eyepref)
    % show_onoff_balance(OriSFdata)
    
    % show correlations 
    % show_correlation_CV_LHI(OriSFdata,rect1,rect2,rand_points,rand_points_interpol,eyepref)
    % show_all_correlation(OriSFdata,rect1,rect2,rand_points,rand_points_interpol,eyepref,1);
    % old test_functionCorrCVSFLHILPIPlot(OriSFdata,rect1,rect2,rand_points,rand_points_interpol,eyepref)
    % all_correlation_with_repeat(OriSFdata,boundary,NumRandPoint,100,eyepref,1);
    
    % Correlation for linear penetrations 
    % show_all_correlation_linear_penetration(OriSFdata,electrode_position,electrode_position*4,eyepref,1)
    all_corr_linear_penetration_repeat(OriSFdata,eyepref,1);
    
    % plot receptive field
    plot_receptive_field(OriSFdata,cortex_RF,electrode_position,eyepref)
    % rf_plot_with_iso_retinotopic_patch(OriSFdata,cortex_RF,electrode_position)
    % rf_plot_with_iso_retinotopic_patch2(OriSFdata,cortex_RF,electrode_position,'hor')
%      plot_rf_eye_dominance(OriSFdata,cortex_RF,ODCrtxPlt_smooth,electrode_position,eyepref)
    
    %plot_rf_linear_lhi_test2(OriSFdata,cortex_RF,ODCrtxPlt_smooth,30,180,'contra')

end 


end 


function show_retinotopy_map(OriSFdata,ODCrtxPlt_smooth)
%   show retintopic grid based on the center of receptive fields 

    rRFcenter = OriSFdata.rRFcenter ;
    cRFcenter = OriSFdata.cRFcenter ;
    
    num_cols = size(ODCrtxPlt_smooth,2)+1; 
    num_rows = size(ODCrtxPlt_smooth,1)+1; 
    ref_cen_row = rRFcenter; 
    ref_cen_col = cRFcenter; 
    ref_cen_col = ref_cen_col - min(ref_cen_col(:));
    ref_cen_row = ref_cen_row - min(ref_cen_row(:));
    ref_cen_col = (ref_cen_col/max(ref_cen_col(:)))*num_cols;
    ref_cen_row = (ref_cen_row/max(ref_cen_row(:)))*num_rows;

    figure
    z = zeros(size(ODCrtxPlt_smooth));
    s_ON = surf(ref_cen_col,ref_cen_row,z);
    s_ON.FaceColor = 'w';
    s_ON.EdgeColor = 'k';
    xlim([1 num_cols])
    ylim([1 num_rows])
    axis square
    view(0, -90);

    od_contour_levels = 1; 
    hold on, contour(ODCrtxPlt_smooth, od_contour_levels, 'k', 'LineWidth', 5);
    title('Retinotopic map','fontsize',20)
end 

function show_max_onoff_response(max_on_res,max_off_res,electrode_position)
    %% Plot ON/OFF max response 
    
    figure, 
    set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
    plot(max_on_res(electrode_position,:),'r')
    hold on 
    plot(max_off_res(electrode_position,:),'b')

end 

function show_onoff_balance(OriSFdata)
    CxOnOffBal = OriSFdata.CxOnOffBal;
    figure, 
    imagesc(CxOnOffBal), colormap('jet'), axis square
    title('ONOFF balance','fontsize',20)
end 

function show_num_projected_afferents(OriSFdata)
    num_projected_aff = OriSFdata.num_projected_aff; 
   figure, 
   imagesc(num_projected_aff), colormap('jet')
   title('Number of projected afferents','fontsize',20)
end 


function show_ori_map(ori_map,ori_map_interpolated,ODCrtxPlt_smooth)

    n_interpolation = 4; 

    figure
    set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);

    subplot(2,2,1);imagesc(ori_map);colormap(gca, 'hsv'); colorbar;axis square;axis off;
    annotation('textbox',[0.03 0.87 0.08 0.08],'String','Raw','EdgeColor','none','fontsize',32)
    subplot(2,2,2);imagesc(ori_map_interpolated);colormap(gca, hsv); colorbar;axis square;axis off;
    annotation('textbox',[0.55 0.87 0.08 0.08],'String','Interpolated','EdgeColor','none','fontsize',32)

    % with OD
    subplot(2,2,3);
    imagesc(ori_map);
    colormap(gca, hsv); colorbar;axis square;axis off;
    hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5);

    subplot(2,2,4);
    imagesc(ori_map_interpolated)
    colormap(gca, 'hsv'); colorbar;axis square;axis off;
    hold on, contour(imresize(ODCrtxPlt_smooth,n_interpolation), 1, 'k', 'LineWidth', 5);
    %title('Ipsi','fontsize',20); 
end 



function show_raw_map(OriSFdata,ODCrtxPlt_smooth,eyepref)
    % Raw Data OriLHI(Local homegeniety index) CV(circular variance) SF(spatial frequency)
    % fucntionRAWLHIOriCVSFPlot(OriLHICVSFreference)
    %OriLHICVSFContra = OriLHICVSFreference; 
    OriMapPref = OriSFdata.ori_map_smooth;
    LHIRaw = OriSFdata.LHI_map;
    CVRaw = OriSFdata.CV_map;
    SFRaw = OriSFdata.SF_map;
    LPIRaw = OriSFdata.LPI_map;
    SF50Raw = OriSFdata.SF50_map;
    %ODCrtxPlt = OriSFdata.ODCrtxPlt;

    figure; clf;
    set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
    %subplot(2,3,1);imagesc(squeeze(OriMapPref(:,:)));title('Orientation');colormap(gca,'hsv'); colorbar;axis square;
    subplot(2,3,1);imagesc(OriMapPref);title('Orientation');colormap(gca,'hsv'); colorbar;axis square;
    hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5);
    subplot(2,3,2); imagesc(LHIRaw); title('LHI'); axis square; colormap(gca,'jet');colorbar
    hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5);
    subplot(2,3,3); imagesc(CVRaw); title('CV');axis square;colormap(gca,'jet');colorbar
    hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5);
    subplot(2,3,4);  imagesc(SFRaw);axis square; title('SFPref');colorbar
    hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5);
    subplot(2,3,5);  imagesc(LPIRaw);axis square; title('LPI');colorbar
    hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5);
    subplot(2,3,6);  imagesc(SF50Raw);axis square; title('SF50'); colorbar; colormap(gca,'jet')
    hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 5); %caxis([.4 1.9])
    annotation('textbox',[0.03 0.87 0.08 0.08],'String',['Raw Data ' eyepref],'EdgeColor','none','fontsize',16)

end 

function show_interpolated_map(OriSFdata,ODCrtxPlt_interpolated,eyepref)
    %% interpolated LHI Ori CV SF
    %   fucntionLHIOriCVSFPlot(OriLHICVSFreference)

    Ori2 = OriSFdata.ori_map_interpolated;
    LHI2 = OriSFdata.LHI_interpolated;
    CV2 = OriSFdata.cv_map_interpolated;
    SF2 = OriSFdata.sf_map_interpolated;
    LPI2 = OriSFdata.LPI_intepolated;
    SF50_2 = OriSFdata.sf50_map_interpolated;
    
    od_contour_levels2 = 1;
    od_contour = 1; 
    od_contour_w = 2; 
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(231),imagesc(Ori2); axis off; title('Orientation map'); axis square; colormap(gca, 'hsv');
    %hold on, contour(imresize(ODCrtxPlt,4), od_contour_levels2, 'k', 'LineWidth', 5);
    hold on, contour(ODCrtxPlt_interpolated, od_contour, 'k', 'LineWidth', od_contour_w);

    subplot(232),imagesc(LHI2);axis off; title('LHI');axis square;colormap(gca, 'jet');
    %hold on, contour(imresize(ODCrtxPlt,4), od_contour_levels2, 'k', 'LineWidth', 5);
    hold on, contour(ODCrtxPlt_interpolated, od_contour, 'k', 'LineWidth', od_contour_w);
    
    subplot(233),imagesc(CV2);axis square;axis off; title('CV');colormap(gca, 'jet');
    %hold on, contour(imresize(ODCrtxPlt,4), od_contour_levels2, 'k', 'LineWidth', 5);
    hold on, contour(ODCrtxPlt_interpolated, od_contour, 'k', 'LineWidth', od_contour_w);
    
    subplot(234),imagesc(SF2);axis off; title('SF');axis square;colormap(gca, 'jet');
    hold on, contour(ODCrtxPlt_interpolated, od_contour, 'k', 'LineWidth', od_contour_w);
    
    subplot(235),imagesc(LPI2);axis off; title('LPI');axis square;colormap(gca, 'jet');
    hold on, contour(ODCrtxPlt_interpolated, od_contour, 'k', 'LineWidth', od_contour_w);
    
    subplot(236),imagesc(SF50_2);axis off; title('SF50');axis square;colormap(gca, 'jet');
    hold on, contour(ODCrtxPlt_interpolated, od_contour, 'k', 'LineWidth', od_contour_w);
    
    
    annotation('textbox',[0.03 0.87 0.08 0.08],'String',eyepref,'EdgeColor','none','fontsize',32)

end 

function show_correlation_CV_LHI(OriSFdata,rect1,rect2,rand_points,rand_points_interpol,eyepref)
    %% CV Vs. LHI
    % fucntionCorrCVLHIPlot(OriLHICVSFreference)
    CV = OriSFdata.CV_map;
    LHI = OriSFdata.LHI_map;
    CV_interpolated = OriSFdata.cv_map_interpolated;
    LHI_interpolated = OriSFdata.LHI_interpolated;

    cvcrop = imcrop(CV,rect1);%figure,imagesc(cvcrop)
    %cvcrop = cvcrop(:);
    cv = cvcrop(rand_points);
    lhicrop = imcrop(LHI,rect1); 
    %lhicrop = lhicrop(:);
    lhi = lhicrop(rand_points);

    cvcrop2 = imcrop(CV_interpolated,rect2);%figure,imagesc(cvcrop)
    %cvcrop = cvcrop(:);
    cv2 = cvcrop2(rand_points_interpol);
    lhicrop2 = imcrop(LHI_interpolated,rect2); 
    %lhicrop = lhicrop(:);
    lhi2 = lhicrop2(rand_points_interpol);

    figure; clf;
    %set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
    subplot(2,2,1);
    %temp = CV; %imcrop(OriCV2,rect);
    plot(lhi, cv,'ko'); lsline;
    [r,p]=corr(lhi', cv');
    title(['Contra  R = ', num2str(r), ' p = ', num2str(p)]);
    xlabel('LHI'), ylabel('CV');   ylim([0 1]);   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')

    subplot(2,2,2);
    %temp = CV_interpolated; %imcrop(newCVmap,rect2);
    plot(lhi2, cv2,'ko'); lsline;
    [r,p]=corr(lhi2', cv2');
    title(['R = ', num2str(r), ' p = ', num2str(p)]);
    xlabel('LHI'), ylabel('CV'); ylim([0 1]);   xlim([0 1]);
    axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    annotation('textbox',[0.03 0.87 0.08 0.08],'String',eyepref,'EdgeColor','none','fontsize',32)
    annotation('textbox',[0.55 0.87 0.08 0.08],'String','Smoothed','EdgeColor','none','fontsize',14)

end 


function [rr_average,pp_average] = all_correlation_with_repeat(OriSFdata,boundary,NumRandPoint,num_repeat,eyepref,show_fig)
    %%  Correlation CV/SF LHI/SF
    
    [rect1,rect2,size_rect1,size_rect2] = remove_boundary(OriSFdata,boundary); 
    
    rr_all = zeros(12,num_repeat);
    pp_all = zeros(12,num_repeat);
    for ii = 1 : num_repeat
        rand_points =   randperm(size_rect1, NumRandPoint);
        rand_points_interpol = randperm(size_rect2, NumRandPoint);
        
        [rr,pp] = show_all_correlation(OriSFdata,rect1,rect2,rand_points,rand_points_interpol,eyepref,0);
    
        rr_all(:,ii) = rr; 
        pp_all(:,ii) = pp; 
    end 
    
    rr_average = mean(rr_all,2); 
    pp_average = mean(pp_all,2); 
    
    %%
    if show_fig == 1
        figure;
        set(gcf , 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
        
        n_bin_pval = 20 ; 
        subplot(2,6,1); histogram(pp_all(1,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(1), pp_average(1)));
        xlabel('CV-SF'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,2); histogram(pp_all(2,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(2), pp_average(2)));
        xlabel('CV-LPI'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,3); histogram(pp_all(3,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(3), pp_average(3)));
        xlabel('CV-SF50'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,4); histogram(pp_all(4,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(4), pp_average(4)));
        xlabel('LHI-SF'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,5); histogram(pp_all(5,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(5), pp_average(5)));
        xlabel('LHI-LPI'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,6); histogram(pp_all(6,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(6), pp_average(6)));
        xlabel('LHI-SF50'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        %   ****** Correlation for Interpolated version ****************
        subplot(2,6,7); histogram(pp_all(7,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(7), pp_average(7)));
        xlabel('CV-SF'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,8); histogram(pp_all(8,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(8), pp_average(8)));
        xlabel('CV-LPI'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,9); histogram(pp_all(9,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(9), pp_average(9)));
        xlabel('CV-SF50'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,10); histogram(pp_all(10,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(10), pp_average(10)));
        xlabel('LHI-SF'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,11); histogram(pp_all(11,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(11), pp_average(11)));
        xlabel('LHI-LPI'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,12); histogram(pp_all(12,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(12), pp_average(12)));
        xlabel('LHI-SF50'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        annotation('textbox',[0.03 0.87 0.5 0.08],'String',[eyepref ' ' num2str(num_repeat) ' random repetitions'],'EdgeColor','none','fontsize',20)
    end
end 

function [rr_average,pp_average] = all_corr_linear_penetration_repeat(OriSFdata,eyepref,show_fig)
    %%  Correlation CV/SF LHI/SF
    num_repeat = 40;
    rr_all = zeros(12,num_repeat);
    pp_all = zeros(12,num_repeat);
    sz_crtx = size(OriSFdata.CV_map,1); 
    ini_point = round(sz_crtx * .1);
    end_point = round(sz_crtx * .9);
    for ii = ini_point : end_point
        electrode_row = ii ;
        electrode_row_inperp = ii * 4; 
        [rr,pp] = show_all_correlation_linear_penetration(OriSFdata,electrode_row,electrode_row_inperp,eyepref,0);
    
        rr_all(:,ii) = rr; 
        pp_all(:,ii) = pp; 
    end 
    
    rr_average = mean(rr_all,2); 
    pp_average = mean(pp_all,2); 
    
    %%
    if show_fig == 1
        figure;
        set(gcf , 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
        
        n_bin_pval = 20 ; 
        subplot(2,6,1); histogram(pp_all(1,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(1), pp_average(1)));
        xlabel('CV-SF'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,2); histogram(pp_all(2,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(2), pp_average(2)));
        xlabel('CV-LPI'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,3); histogram(pp_all(3,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(3), pp_average(3)));
        xlabel('CV-SF50'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,4); histogram(pp_all(4,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(4), pp_average(4)));
        xlabel('LHI-SF'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,5); histogram(pp_all(5,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(5), pp_average(5)));
        xlabel('LHI-LPI'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,6); histogram(pp_all(6,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(6), pp_average(6)));
        xlabel('LHI-SF50'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        %   ****** Correlation for Interpolated version ****************
        subplot(2,6,7); histogram(pp_all(7,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(7), pp_average(7)));
        xlabel('CV-SF'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,8); histogram(pp_all(8,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(8), pp_average(8)));
        xlabel('CV-LPI'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,9); histogram(pp_all(9,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(9), pp_average(9)));
        xlabel('CV-SF50'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,10); histogram(pp_all(10,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(10), pp_average(10)));
        xlabel('LHI-SF'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,11); histogram(pp_all(11,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(11), pp_average(11)));
        xlabel('LHI-LPI'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,6,12); histogram(pp_all(12,:),n_bin_pval)
        title(sprintf(' R = %.4f \nP = %.4f', rr_average(12), pp_average(12)));
        xlabel('LHI-SF50'), axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        annotation('textbox',[0.03 0.87 0.48 0.08],'String',[eyepref ' ' num2str(num_repeat) ' linear repetitions'],'EdgeColor','none','fontsize',20)
    end
end 


function [rect1,rect2,size_rect1,size_rect2] = remove_boundary(OriSFdata,boundary)
    % remove boundary
    LHI = OriSFdata.LHI_map; 
    LHI_interpolated = OriSFdata.LHI_interpolated; 
    bound1 = round(boundary * size(LHI,1));
    bound2 = round(boundary * size(LHI_interpolated,1));
    rect1 = [bound1 bound1 size(LHI,1)-bound1 size(LHI,2)-bound1];
    rect2 = [bound2 bound2 size(LHI_interpolated,1)-bound2 size(LHI_interpolated,2)-bound2];
    
    size_rect1 = (size(LHI,1)- 2*bound1)^2 ;
    size_rect2 = (size(LHI_interpolated,1)- 2*bound2)^2 ;
end


%%         - Figure, Receptive Field Reference
function plot_receptive_field(OriSFdata,cortex_RF,electrodePos,eyepref)

% LHIRaw = OriSFdata.LHI_map;
% CVRaw = OriSFdata.CV_map;
    
allCXrfONReference = cortex_RF.allCXrfON ;
allCXrfOFFReference = cortex_RF.allCXrfOFF ;

OriHistRFReference = OriSFdata.OriTunPolarPlot;
OriBinRFReference = OriSFdata.angles;

SFHistRFReference = OriSFdata.allsfcurve_dog;
sfpref = OriSFdata.SF_map; 
sf_bin_deg = OriSFdata.sf_bin_deg; 
SF50Raw = OriSFdata.SF50_map;

figure;clf
set(gcf,'position',[10         300        1200         350])
%set(fg11,'position',[10         677        1673         300])
xs = linspace(0.02,0.99-(1/32),32);
width = 1/(32*1.3);
height = 1/(3*1.35);
caxisVal = 1;

for ii = 2 : 2 : size(sfpref, 1)
    RFrefON = allCXrfONReference{electrodePos,ii}; 
    RFrefOFF = allCXrfOFFReference{electrodePos,ii}; 
    RFreONOFF = RFrefON + RFrefOFF; 
    MaxValrf = max(max(abs(RFrefON(:))),max(abs(RFrefOFF(:)))); 
    
    axes('Position',[xs(ii/2+1),0.7,width,height])
    RFreONOFFnorm = RFreONOFF / MaxValrf; 
    imagesc(RFreONOFFnorm)
    %xlim([100 300]),ylim([100 300])
    title(num2str(ii))
    caxis([-caxisVal caxisVal])
    axis off
    axis square
    
    axes('Position',[xs(ii/2+1),0.55,width,height])
    RFrefONnorm = RFrefON / MaxValrf; 
    imagesc(RFrefONnorm)
    caxis([-caxisVal caxisVal])
    axis off
    axis square
    
    axes('Position',[xs(ii/2+1),0.4,width,height])
    RFrefOFFnorm = RFrefOFF / MaxValrf; 
    imagesc(RFrefOFFnorm)
        caxis([-caxisVal caxisVal])
    axis off
    axis square
    
    axes('Position',[xs(ii/2+1),0.2,width,height])
    ori_tuning = OriHistRFReference{electrodePos,ii};
    angles = OriBinRFReference{electrodePos,ii};%{electrodePos,ii};
    ph = polar1(angles,ori_tuning','k');
    %title(sprintf('cv:%.2f \n lhi:%.2f',CVRaw(electrodePos,ii),LHIRaw(electrodePos,ii)))

    axes('Position',[xs(ii/2+1),0.00,width,height])
    SFPlot = SFHistRFReference{electrodePos,ii};
    %SFPlot = (SFHistRFReference(electrodePos,ii,:));
    plot(SFPlot,'lineWidth',2);
    axis square
    set(gca, 'XScale', 'log');
    set(gca,'box','off','Tickdir','out','YTickLabel',[],'XTickLabel',[])
    title(sprintf('%.2f',SF50Raw(electrodePos,ii)))
    
    annotation('textbox',[0.01 0.88 0.1 0.1],'String',[eyepref ' ,row = ' num2str(electrodePos)],'EdgeColor','none','fontsize',8)
end

colormap('jet')
end 


function rf_plot_with_iso_retinotopic_patch(OutputRFCrtx,cortex_RF,electrodePos)

allCXrfONReference = cortex_RF.allCXrfON ;
allCXrfOFFReference = cortex_RF.allCXrfOFF ;

OriHistRFReference = OutputRFCrtx.OriTunPolarPlot;
OriBinRFReference = OutputRFCrtx.angles;

R_iso_retinotopy =OutputRFCrtx.R_iso_retinotopy;

rRFcenter = OutputRFCrtx.rRFcenter;
cRFcenter = OutputRFCrtx.cRFcenter;

%
figure;clf;
set(gcf,'position',[10         677        1673         300])
xs = linspace(0.02,0.99-(1/32),32);
width = 1/(32*1.3);
height = 1/(3*1.35);
caxisVal = 1;

width_rf = size(allCXrfONReference{1},1); 

for ii = 2:2:60
    RFrefON = allCXrfONReference{electrodePos,ii}; 
    RFrefOFF = allCXrfOFFReference{electrodePos,ii}; 
    RFreONOFF = RFrefON + RFrefOFF; 
    MaxValrf = max(max(abs(RFrefON(:))),max(abs(RFrefOFF(:)))); 
    
    axes('Position',[xs(ii/2+1),0.7,width,height])
    RFreONOFFnorm = RFreONOFF / MaxValrf; 
    imagesc(RFreONOFFnorm)
    title(num2str(ii))
    caxis([-caxisVal caxisVal])
    axis off
    axis square
    
    axes('Position',[xs(ii/2+1),0.55,width,height])
    RFrefONnorm = RFrefON / MaxValrf; 
    imagesc(RFrefONnorm)
    caxis([-caxisVal caxisVal])
    axis off
    axis square
    
    axes('Position',[xs(ii/2+1),0.4,width,height])
    RFrefOFFnorm = RFrefOFF / MaxValrf; 
    imagesc(RFrefOFFnorm)
        caxis([-caxisVal caxisVal])
    axis off
    axis square
    
    axes('Position',[xs(ii/2+1),0.2,width,height])
    tuning = OriHistRFReference{electrodePos,ii};
    angles = OriBinRFReference{electrodePos,ii};%{electrodePos,ii};
    ph = polar1(angles,tuning','k');

    axes('Position',[xs(ii/2+1),0.00,width,height])
    axis off, axis square
    radius = R_iso_retinotopy(electrodePos,ii);
    x_center = cRFcenter(electrodePos,ii);
    y_center = rRFcenter(electrodePos,ii);
    hold on, viscircles([x_center y_center],radius,'Color','black')

    hold on
    w=[1 1];
    x=[1 width_rf];
    y=[width_rf width_rf];
    z=[width_rf 1];
    Points = [w;x;y;z;w];   % in desired order
    plot( Points(:,1), Points(:,2), 'k');   % draw bounding box
end
colormap('jet')

end 



%   - Figure, Receptive Field Reference vertical/horizontal 
function rf_plot_with_iso_retinotopic_patch2(OutputRFCrtxReference,cortex_RF,electrodePos,dir)

allCXrfONReference = cortex_RF.allCXrfON ;
allCXrfOFFReference = cortex_RF.allCXrfOFF ;

OriHistRFReference = OutputRFCrtxReference.OriTunPolarPlot;
OriBinRFReference = OutputRFCrtxReference.angles;

R_iso_retinotopy = OutputRFCrtxReference.R_iso_retinotopy;
rRFcenter = OutputRFCrtxReference.rRFcenter;
cRFcenter = OutputRFCrtxReference.cRFcenter;

figure;clf
set(gcf,'position',[10         677        1673         300])
%set(fg11,'position',[10         677        1673         225])
%set(fg11,'position',[10         300        1673         450])
xs = linspace(0.02,0.99-(1/32),32);
width = 1/(32*1.3);
height = 1/(3*1.35);
caxisVal = 1;

width_rf = size(allCXrfONReference{1},1); 


for ii = 2:2:60
    
    if strcmp(dir,'hor')
        ind_row = electrodePos;
        ind_col = ii; 
    elseif strcmp(dir,'ver')
        ind_row = ii; 
        ind_col = electrodePos; 
    end 
    
    RFrefON = allCXrfONReference{ind_row,ind_col}; 
    RFrefOFF = allCXrfOFFReference{ind_row,ind_col}; 
    RFreONOFF = RFrefON + RFrefOFF; 
    MaxValrf = max(max(abs(RFrefON(:))),max(abs(RFrefOFF(:)))); 
    
    axes('Position',[xs(ii/2+1),0.7,width,height])
    RFreONOFFnorm = RFreONOFF / MaxValrf; 
    imagesc(RFreONOFFnorm)
    %xlim([100 300]),ylim([100 300])
    title(num2str(ii))
    caxis([-caxisVal caxisVal])
    axis off
    axis square
    
    axes('Position',[xs(ii/2+1),0.55,width,height])
    RFrefONnorm = RFrefON / MaxValrf; 
    imagesc(RFrefONnorm)
    caxis([-caxisVal caxisVal])
    axis off
    axis square
    
    axes('Position',[xs(ii/2+1),0.4,width,height])
    RFrefOFFnorm = RFrefOFF / MaxValrf; 
    imagesc(RFrefOFFnorm)
        caxis([-caxisVal caxisVal])
    axis off
    axis square
    
    axes('Position',[xs(ii/2+1),0.2,width,height])
    ori_tuning = OriHistRFReference{ind_row,ind_col};
    angles = OriBinRFReference{ind_row,ind_col};%{electrodePos,ii};
    ph = polar1(angles,ori_tuning','k');

%     r = sfpref(electrodePos,ii);
%     CyclePix = size(RFreONOFF,1) / r; % Cycle in pixel
%     CycleDeg = CyclePix / deg2pix;
%     CyclePerDeg = 1/CycleDeg;
    
    axes('Position',[xs(ii/2+1),0.00,width,height])
    axis off, axis square
    radius = R_iso_retinotopy(ind_row,ind_col);
    x_center = cRFcenter(ind_row,ind_col);
    y_center = rRFcenter(ind_row,ind_col);
    hold on, viscircles([x_center y_center],radius,'Color','black')

    hold on
    w=[1 1];
    x=[1 width_rf];
    y=[width_rf width_rf];
    z=[width_rf 1];
    Points = [w;x;y;z;w];   % in desired order
    plot( Points(:,1), Points(:,2), 'k');   % draw bounding box
    set(gca,'ydir','reverse')
end

colormap('jet')
end 


function plot_rf_eye_dominance(OriSFdata, cortex_RF, ODCrtxPlt_smooth, electrodePos, eyepref)

LHIRaw = OriSFdata.LHI_map;
CVRaw = OriSFdata.CV_map;
% ODI = OriSFdata.ODI; 
LPI_map = OriSFdata.LPI_map;

allCXrfONReference = cortex_RF.allCXrfON ;
allCXrfOFFReference = cortex_RF.allCXrfOFF ;

OriHistRFReference = OriSFdata.OriTunPolarPlot;
OriBinRFReference = OriSFdata.angles;

SFHistRFReference = OriSFdata.allsfcurve_dog;
sfpref = OriSFdata.SF_map; 
sf_bin_deg = OriSFdata.sf_bin_deg; 
SF50Raw = OriSFdata.SF50_map;

figure;clf
set(gcf,'position',[10         677        1673         300])
%set(fg11,'position',[10         677        1673         225])
%set(fg11,'position',[10         300        1673         450])
xs = linspace(0.02,0.99-(1/32),32);
width = 1/(32*1.3);
height = 1/(3*1.35);
caxisVal = 1;

for ii = 2:2:60
    RFrefON = allCXrfONReference{electrodePos,ii}; 
    RFrefOFF = allCXrfOFFReference{electrodePos,ii}; 
    RFreONOFF = RFrefON + RFrefOFF; 
    MaxValrf = max(max(abs(RFrefON(:))),max(abs(RFrefOFF(:)))); 
    
    axes('Position',[xs(ii/2+1),0.6,width,height])
    RFreONOFFnorm = RFreONOFF / MaxValrf; 
    imagesc(RFreONOFFnorm)
    %xlim([100 300]),ylim([100 300])
    
    caxis([-caxisVal caxisVal])
    axis off
    axis square
    if ODCrtxPlt_smooth(electrodePos,ii) < .5 % ipsi
        title(num2str(ii),'color',[0.9290 0.6940 0.1250]) % orange 
    else 
        title(num2str(ii),'color','k')
    end 

    axes('Position',[xs(ii/2+1),0.35,width,height])
    ori_tuning = OriHistRFReference{electrodePos,ii};
    angles = OriBinRFReference{electrodePos,ii};%{electrodePos,ii};
    ph = polar1(angles,ori_tuning','k');
    title(sprintf('cv:%.2f \n lhi:%.2f',CVRaw(electrodePos,ii),LHIRaw(electrodePos,ii)))
%     title(sprintf('cv:%.2f \n lhi:%.2f \n od:%.2f ',...
%         CVRaw(electrodePos,ii),LHIRaw(electrodePos,ii),ODI(electrodePos,ii)))
    
    axes('Position',[xs(ii/2+1),0.1,width,height])
    SFPlot = SFHistRFReference{electrodePos,ii};
    %SFPlot = (SFHistRFReference(electrodePos,ii,:));
    plot(SFPlot,'lineWidth',2);
    axis square
    set(gca, 'XScale', 'log');
    set(gca,'box','off','Tickdir','out','YTickLabel',[],'XTickLabel',[])
    %title(sprintf('sf:%.2f \n %.2f',sfpref(electrodePos,ii), SF50Raw(electrodePos,ii)))
    title(sprintf('sf:%.2f \n 50:%.2f\n lpi:%.2f',sfpref(electrodePos,ii), SF50Raw(electrodePos,ii),LPI_map(electrodePos,ii)))

    annotation('textbox',[0.01 0.88 0.1 0.1],'String',[eyepref ',row = ' num2str(electrodePos)],'EdgeColor','none','fontsize',14)
end

colormap('jet')
end 


