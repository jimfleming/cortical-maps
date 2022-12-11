function [RF_data,OriSFdata] = measure_rf_stats(eyepref,rf_input,data_input,reference_ori_map,RetinaRF,rowRange,colRange,...
    IndRetinaRfspace,rAffSpread,RetONOFFsorted,ONOFFCrtxPlt,ODCrtxPlt,ODCrtxPlt_smooth,electrode_position,sigma_LHI,n_ori_smooth,sf_lSamp,pix2deg,weight_coef_RF,alpha,debug)

%   measure the statistics of input RFs for both eyes 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
sdspread = rAffSpread;
%idxAffRfspace = IndRetinaRfspace;
SelectionCondition = 0; % should come from input data
radius_rf = 0; 
% data from both eyes is needed to modify the rfs based on the dominant eye in that cortical location

rf_on_all = rf_input.allCXrfON ;
rf_off_all = rf_input.allCXrfOFF ;
% row_rf_center_contra = data_input.rRFcenter;
% col_rf_center_contra = data_input.cRFcenter ;
% pref_orimap_contra = data_input.OriPreferred ;
% ODI_contra = data_contra_primord.ODI;


%%
%   clear rfAllPrefEye
od = zeros(size(RetONOFFsorted));
Ori_map = zeros(size(od));
CV_map = zeros(size(od));
SF_map = zeros(size(Ori_map));
SF50_map = zeros(size(Ori_map));
LPI_map = zeros(size(Ori_map));

allsfcurve = zeros([size(od) length(sf_lSamp)]);
alloricurve = zeros([size(od) length(sf_lSamp)]);
allangles = zeros([size(od) length(sf_lSamp)]);
OriTunPolarPlot{size(od,1),size(od,2)} = [];
sf_tuning_curve_all{size(od,1),size(od,2)} = [];
sf_tuning_dog_all{size(od,1),size(od,2)} = [];
ori_bin_plot_all{size(od,1),size(od,2)} = [];

rfAllPrefEye{size(od,1),size(od,2)} = []; 
rfAllPrefEyeON{size(od,1),size(od,2)} = []; 
rfAllPrefEyeOFF{size(od,1),size(od,2)} = []; 
OriPrefEye = zeros(size(od));
num_projected_aff = zeros(size(od));  
max_on_response = zeros(size(od)); 
max_off_response = zeros(size(od));
max_on_response_norm = zeros(size(od)); 
max_off_response_norm = zeros(size(od)); 
R_iso_retinotopy = zeros(size(od));
CxOnOffBal = zeros(size(od));
zero_polarity = zeros(size(od));

%% Main 
 
for ii = rowRange 
    for jj = colRange 
        %
        allCXrfTempON1  = rf_on_all{ii,jj};
        allCXrfTempOFF1 = rf_off_all{ii,jj};

        %
        RFSpaceSim =   allCXrfTempON1 + allCXrfTempOFF1 ;
        SingleRFNorm = RFSpaceSim / max(abs(RFSpaceSim(:)));
%         [ori_tuning_resp,sf_tuning_resp,sf_tuning_resp_dog,cv,ori_pref,sf_pref_deg,sf50_deg,...
%             LPI,angles_bin,num_cycles_bin_deg] = fft_ori_sf_tuning4(SingleRFNorm,sf_lSamp,pix2deg,0);
%           [ori_tuning_resp,sf_tuning_resp,sf_tuning_resp_dog,cv,ori_pref,sf_pref_deg,sf50_deg,LPI,angles_bin,num_cycles_bin_deg] = ... 
%                fft_ori_sf_tuning5(SingleRFNorm,sf_lSamp,pix2deg,0); 
        [ori_tuning_resp,sf_tuning_resp,sf_tuning_resp_dog,cv,ori_pref,sf_pref_deg,sf50_deg,LPI,angles_bin,num_cycles_bin_deg] = ... 
               fft_ori_sf_tuning6(SingleRFNorm,sf_lSamp, 0, pix2deg, 0); 
           
       if debug
%            debug_RF(SingleRFNorm,angles_bin,ori_tuning_resp,ori_pref,sf_tuning_resp_dog,sf_lSamp,sf_pref_deg,pix2deg)
           % debug_RF2(allCXrfTempON,allCXrfTempOFF,SingleRFNorm,aff_center_weight_plot,angles_bin,ori_tuning_resp,ori_pref,sf_bin_deg,sf_tuning_resp_dog,sf_pref_deg,sf50_deg)
       end
            
        max_on_response(ii,jj) = max(allCXrfTempON1(:)); 
        max_off_response(ii,jj) = max(abs(allCXrfTempOFF1(:))); 
        
        
        ind_on_rf = SingleRFNorm(:) > 0 ; 
%         if ~isempty(ind_on_rf)
        try
            max_on_response_norm(ii,jj) = max(SingleRFNorm(ind_on_rf)); 
        catch 
            max_on_response_norm(ii,jj) = 0; 
        end 
        
        ind_off_rf = SingleRFNorm(:) < 0 ;
%         if ~isempty(ind_off_rf)
        try
            max_off_response_norm(ii,jj) = max(abs(SingleRFNorm(ind_off_rf)));
        catch
            max_off_response_norm(ii,jj) = 0;
        end
        
        % for presentation purposes 
        CxOnOffBal(ii,jj) = ONOFF_balance(allCXrfTempON1,allCXrfTempOFF1); 
        
        % iso retinotopic patch 
        % alpha = .5; % 50% of the maximum value is set to find the 
        % radius_rf = measure_radius_rf(rowRFcenter,colRFcenter,allCXrfTempON1,allCXrfTempOFF1,alpha); 
        
        rfAllPrefEyeON{ii,jj} = allCXrfTempON1; 
        rfAllPrefEyeOFF{ii,jj} = allCXrfTempOFF1; 
        rfAllPrefEye{ii,jj} = SingleRFNorm; 
        
        OriPrefEye(ii,jj) = ori_pref; 
        CV_map(ii,jj) = cv;
        SF_map(ii,jj) = sf_pref_deg;
        SF50_map(ii,jj) = sf50_deg;
        LPI_map(ii,jj) = LPI; 
        
        OriTunPolarPlot{ii,jj} = ori_tuning_resp;
        %alloricurve(ii,jj,:)= ori_tuning_resp;
        sf_tuning_curve_all{ii,jj} = sf_tuning_resp; 
        %allsfcurve(ii,jj,:) = sf_tuning_resp;
        sf_tuning_dog_all{ii,jj} = sf_tuning_resp_dog;
        ori_bin_plot_all{ii,jj} = angles_bin; 
        %allangles(ii,jj,:) = angles_bin;
        num_projected_aff(ii,jj) = sum(SelectionCondition(:)); 
        R_iso_retinotopy(ii,jj) = radius_rf; 
        
    end
end

%% ori smooth  
% n_ori_smooth = 4;
[ori_map_smooth,ori_map_interpolated] = smooth_ori_euler(OriPrefEye,ODCrtxPlt_smooth,n_ori_smooth,0,0); 

% LHI measurement 
LHI_map = measure_LHI(ori_map_smooth,sigma_LHI); 
sigma_LHI_smooth = sigma_LHI * n_ori_smooth;
LHI_interpolated = measure_LHI(ori_map_interpolated,sigma_LHI_smooth); 

%% Map interpolation 
% n_interpolation = 2; 
% cv_map_interpolated = interp2(CV_map,n_interpolation, 'linear');
% sf_map_interpolated = interp2(SF_map,n_interpolation,'linear');
% LPI_intepolated = interp2(LPI_map,n_interpolation,'linear');
% LPI_intepolated(LPI_intepolated>1) = 1;
% sf50_map_interpolated = interp2(SF50_map,n_interpolation,'linear');
cv_map_interpolated = imresize(CV_map,n_ori_smooth);
sf_map_interpolated = imresize(SF_map,n_ori_smooth);
LPI_intepolated = imresize(LPI_map,n_ori_smooth);
%LPI_intepolated(LPI_intepolated>1) = 1;
sf50_map_interpolated = imresize(SF50_map,n_ori_smooth);
%% 
OriSFdata.ori_map_smooth = ori_map_smooth;
OriSFdata.ori_pref_rf = OriPrefEye;
OriSFdata.LHI_map = LHI_map; 
OriSFdata.LPI_map = LPI_map; 
OriSFdata.SF_map = SF_map;
OriSFdata.SF50_map = SF50_map;
OriSFdata.CV_map = CV_map;

OriSFdata.ori_map_interpolated = ori_map_interpolated;
OriSFdata.LHI_interpolated = LHI_interpolated; 
OriSFdata.cv_map_interpolated = cv_map_interpolated; 
OriSFdata.sf_map_interpolated = sf_map_interpolated; 
OriSFdata.LPI_intepolated = LPI_intepolated; 
OriSFdata.sf50_map_interpolated = sf50_map_interpolated; 

OriSFdata.allsfcurve = sf_tuning_curve_all;
OriSFdata.allsfcurve_dog = sf_tuning_dog_all;
%OriSFdata.alloricurve = alloricurve;
%OriSFdata.allangles = allangles;
OriSFdata.angles = ori_bin_plot_all;
OriSFdata.sf_bin_deg = num_cycles_bin_deg; 
OriSFdata.OriTunPolarPlot = OriTunPolarPlot;
OriSFdata.sdspread = sdspread; 
OriSFdata.sf_lSamp = sf_lSamp; 

OriSFdata.num_projected_aff = num_projected_aff; 
OriSFdata.max_on_response = max_on_response; 
OriSFdata.max_off_response = max_off_response;
OriSFdata.max_on_response_norm = max_on_response_norm; 
OriSFdata.max_off_response_norm = max_off_response_norm;
OriSFdata.CxOnOffBal = CxOnOffBal;
% OriSFdata.R_iso_retinotopy = R_iso_retinotopy; 
% OriSFdata.rRFcenter = row_rf_center_contra; % the center does not change comparing to the reference RF
% OriSFdata.cRFcenter = col_rf_center_contra;
% OriSFdata.zero_polarity = zero_polarity;

% OriSFdata.ODI = ODI_contra ; 
% OriSFdata.ODI_contra = ODI_contra ; 
% OriSFdata.ODI_ipsi = ODI_ipsi ; 


% RF_data.rfAllPrefEye = rfAllPrefEye; 
RF_data.allCXrfON = rfAllPrefEyeON; 
RF_data.allCXrfOFF = rfAllPrefEyeOFF;

    if debug 
%         plot_receptive_field(OriSFdata,RF_data,electrode_position,reference_ori_map)
        plot_rf_eye_dominance_mature_cortex(OriSFdata,RF_data,ODCrtxPlt_smooth,electrode_position,eyepref)
        
%         figure, histogram(SF50_map(30,:),10), xlabel('SF50 range','fontsize',20)
%         set(gca,'box','off','Tickdir','out','FontSize',20)
        
        cv = OriSFdata.CV_map(electrode_position,:); 
        onoff_bal = OriSFdata.CxOnOffBal(electrode_position,:); 
        figure, plot(onoff_bal, cv,'ko'); lsline;
        [r,p]=corr(onoff_bal', cv');
        title(['R = ', num2str(r), ' p = ', num2str(p)]);
        axis square, xlabel('ONOFF balance'), ylabel('CV(OS)'), 
        set(gca,'box' ,'off','TickDir','OUT')
    end 

end


%% debug
function debug_RF(rf_input,angles_bin,ori_tuning_resp,ori_pref,sf_tuning_resp,sf_lSamp,sf_pref_deg,pix2deg)
    %debug
    figure;
    set(gcf,'Units','Normalized','Outerposition',[0 0 .8 .8])
    subplot(131),imagesc(rf_input),colormap jet;colorbar,axis square
    
    subplot(132),polar1(angles_bin,ori_tuning_resp','k');hold on;axis square
    title(sprintf('Ori : %.2f ',ori_pref))
    
    sfallnorm = sf_tuning_resp - min(sf_tuning_resp(:)) ;
    sfallnorm = sfallnorm / max(sfallnorm(:));
    YCyclePix = size(rf_input,1) ./ (1:max(sf_lSamp)); % Cycle in pixel
    YCyclePerDeg = YCyclePix ./ pix2deg;
    YCyclePerDeg = 1 ./YCyclePerDeg;
    subplot(133),plot(YCyclePerDeg,sfallnorm,'lineWidth',2);axis square
    title(sprintf('SF : %.2f CPD',sf_pref_deg))
    xlabel('Cyc/Deg')
    set(gca,'box','off','Tickdir','out');
end
