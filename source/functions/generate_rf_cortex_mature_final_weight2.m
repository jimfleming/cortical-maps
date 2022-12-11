function [RF_data,OriSFdata] = generate_rf_cortex_mature_final_weight2(eyepref,rf_input,data_input,reference_ori_map,RetinaRF,rowRange,colRange,...
    idxAffRfspace,sdspread,RetONOFFsorted,ONOFFCrtxPlt,ODCrtxPlt,ODCrtxPlt_smooth,electrode_position,sigma_LHI,n_ori_smooth,sf_lSamp,pix2deg,weight_coef_RF,alpha,debug)
%   modifying receptive field weights  based on the orientaiton map 
%
%   eyepref = 'contra'/'ipsi'; % the eye for that new rf will be calculated
%   developed_reference_ori_map : the reference orientation map 
%   alpha : for iso retinotopic measurement : percentage of maximum reponse in RF space to look for from the center of RF
%
%   Sohrab note     :   The old name of this function was "rf_modify_eye_weight1" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

rf_on_all = rf_input.allCXrfON ;
rf_off_all = rf_input.allCXrfOFF ;

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
max_on_response = zeros(size(od)); 
max_off_response = zeros(size(od)); 
max_on_response_norm = zeros(size(od)); 
max_off_response_norm = zeros(size(od)); 
CxOnOffBal = zeros(size(od));
zero_polarity = zeros(size(od));

%% Main 
%%
for ii = rowRange 
    for jj = colRange 
        %
        rfON = rf_on_all{ii,jj};
        rfOFF = rf_off_all{ii,jj};
        [LWON, LWOFF ] = modify_synaptic_weight(rfON, rfOFF, weight_coef_RF(ii,jj), ONOFFCrtxPlt(ii,jj), 0); 
        
        %   tuning curve measurements 
        allCXrfTempON1 = rfON * LWON; 
        allCXrfTempOFF1 = rfOFF * LWOFF; 
        
        RFSpaceSim =   allCXrfTempON1 + allCXrfTempOFF1 ;
        SingleRFNorm = RFSpaceSim / max(abs(RFSpaceSim(:)));
        [ori_tuning_resp,sf_tuning_resp,sf_tuning_resp_dog,cv,ori_pref,sf_pref_deg,sf50_deg,LPI,angles_bin,num_cycles_bin_deg] = ...
            fft_ori_sf_tuning6(SingleRFNorm, sf_lSamp, 0, pix2deg, 0);
        
        max_on_response(ii,jj) = max(allCXrfTempON1(:)); 
        max_off_response(ii,jj) = max(abs(allCXrfTempOFF1(:))); 
        
        [Max_on, Max_off] = find_max_resp_norm(SingleRFNorm); 
        max_on_response_norm(ii,jj) = Max_on; 
        max_off_response_norm(ii,jj) = Max_off;  
        CxOnOffBal(ii,jj) = ONOFF_balance(allCXrfTempON1,allCXrfTempOFF1); 

        rfAllPrefEyeON{ii,jj} = allCXrfTempON1; 
        rfAllPrefEyeOFF{ii,jj} = allCXrfTempOFF1; 
        rfAllPrefEye{ii,jj} = SingleRFNorm; 
        
        OriPrefEye(ii,jj) = ori_pref; 
        CV_map(ii,jj) = cv;
        SF_map(ii,jj) = sf_pref_deg;
        SF50_map(ii,jj) = sf50_deg;
        LPI_map(ii,jj) = LPI; 
        
        OriTunPolarPlot{ii,jj} = ori_tuning_resp;
        % alloricurve(ii,jj,:)= ori_tuning_resp;
        sf_tuning_curve_all{ii,jj} = sf_tuning_resp; 
        % allsfcurve(ii,jj,:) = sf_tuning_resp;
        sf_tuning_dog_all{ii,jj} = sf_tuning_resp_dog;
        ori_bin_plot_all{ii,jj} = angles_bin; 
        % allangles(ii,jj,:) = angles_bin;
        
    end
end

%% ori smooth  
[ori_map_smooth,ori_map_interpolated] = smooth_ori_euler(OriPrefEye,ODCrtxPlt_smooth,n_ori_smooth,0,0); 
    % LHI measurement 
    LHI_map = measure_LHI(ori_map_smooth,sigma_LHI); 
    sigma_LHI_smooth = sigma_LHI * n_ori_smooth;
    LHI_interpolated = measure_LHI(ori_map_interpolated,sigma_LHI_smooth); 

    %% Map interpolation 
    cv_map_interpolated = imresize(CV_map,n_ori_smooth);
    sf_map_interpolated = imresize(SF_map,n_ori_smooth);
    LPI_intepolated = imresize(LPI_map,n_ori_smooth);
    sf50_map_interpolated = imresize(SF50_map,n_ori_smooth);
    %% outputs 
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
    OriSFdata.angles = ori_bin_plot_all;
    OriSFdata.sf_bin_deg = num_cycles_bin_deg; 
    OriSFdata.OriTunPolarPlot = OriTunPolarPlot;
    OriSFdata.sdspread = sdspread; 
    OriSFdata.sf_lSamp = sf_lSamp; 

    OriSFdata.max_on_response = max_on_response; 
    OriSFdata.max_off_response = max_off_response;
    OriSFdata.max_on_response_norm = max_on_response_norm; 
    OriSFdata.max_off_response_norm = max_off_response_norm;
    OriSFdata.CxOnOffBal = CxOnOffBal;
    OriSFdata.zero_polarity = zero_polarity;

    RF_data.allCXrfON = rfAllPrefEyeON; 
    RF_data.allCXrfOFF = rfAllPrefEyeOFF;
end


%% function
function [max_on_response_norm, max_off_response_norm] = find_max_resp_norm(SingleRFNorm)
    ind_on_rf = SingleRFNorm(:) > 0 ;
    max_on_response_norm = max(SingleRFNorm(ind_on_rf));
    if isempty(max_on_response_norm)
        max_on_response_norm = 0;
    end

    ind_off_rf = SingleRFNorm(:) < 0 ;
    max_off_response_norm = max(abs(SingleRFNorm(ind_off_rf)));
    if isempty(max_off_response_norm)
        max_off_response_norm = 0;
    end
end 

