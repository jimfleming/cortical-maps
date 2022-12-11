function [RF_data,OriSFdata] = generate_rf_cortex_mature_final_weight(eyepref,rf_input,data_input,reference_ori_map,RetinaRF,rowRange,colRange,...
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
% row_rf_center_contra = data_input.rRFcenter;
% col_rf_center_contra = data_input.cRFcenter ;

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
% num_projected_aff = zeros(size(od));  
max_on_response = zeros(size(od)); 
max_off_response = zeros(size(od)); 
max_on_response_norm = zeros(size(od)); 
max_off_response_norm = zeros(size(od)); 
% R_iso_retinotopy = zeros(size(od));
CxOnOffBal = zeros(size(od));
zero_polarity = zeros(size(od));

%% Main 
 
for ii = rowRange 
    for jj = colRange 
        %%
        rfON  = rf_on_all{ii,jj};
        rfOFF = rf_off_all{ii,jj};
        
        % Modifying synaptic weights based on pw locations  
        max_resp_on = sum(rfON(:));         % max(allCXrfTempON(:));
        max_resp_off = sum(abs(rfOFF(:)));  % max(abs(allCXrfTempOFF(:)));

        if max_resp_on ~= 0 && max_resp_off ~= 0
            weight_rf = weight_coef_RF(ii,jj); % .9 to make it lower than the dominant polarity
            if ONOFFCrtxPlt(ii,jj) == 1
                LWON =  1;
                LWOFF =  weight_rf ;%* (max_resp_on/max_resp_off);
            elseif ONOFFCrtxPlt(ii,jj) == -1
                LWON  =  weight_rf ;%* (max_resp_off/max_resp_on);
                LWOFF =  1;
            end
        else 
            LWON = 1; 
            LWOFF =  1;
            zero_polarity(ii,jj) = 1;  % to measure the number of cortical locations that non-dominant polarity is zero
        end 
        %%   tuning curve measurements 
        allCXrfTempON1 = rfON * LWON; 
        allCXrfTempOFF1 = rfOFF * LWOFF; 
        
        RFSpaceSim =   allCXrfTempON1 + allCXrfTempOFF1 ;
        SingleRFNorm = RFSpaceSim / max(abs(RFSpaceSim(:)));
        [ori_tuning_resp,sf_tuning_resp,sf_tuning_resp_dog,cv,ori_pref,sf_pref_deg,sf50_deg,LPI,angles_bin,num_cycles_bin_deg] = ...
            fft_ori_sf_tuning6(SingleRFNorm,sf_lSamp,1,pix2deg,0);
        
        if debug
            % debug_RF(SingleRFNorm,angles_bin,ori_tuning_resp,ori_pref,sf_tuning_resp_dog,sf_lSamp,sf_pref_deg,pix2deg)
            % debug_RF2(allCXrfTempON,allCXrfTempOFF,SingleRFNorm,aff_center_weight_plot,angles_bin,ori_tuning_resp,ori_pref,sf_bin_deg,sf_tuning_resp_dog,sf_pref_deg,sf50_deg)
        end
        %%
        max_on_response(ii,jj) = max(allCXrfTempON1(:)); 
        max_off_response(ii,jj) = max(abs(allCXrfTempOFF1(:))); 
        
        ind_on_rf = SingleRFNorm(:) > 0 ; 
        try
            max_on_response_norm(ii,jj) = max(SingleRFNorm(ind_on_rf)); 
        catch
            max_on_response_norm(ii,jj) = 0; 
        end 
        
        ind_off_rf = SingleRFNorm(:) < 0 ;
        try
            max_off_response_norm(ii,jj) = max(abs(SingleRFNorm(ind_off_rf))); 
        catch
            max_off_response_norm(ii,jj) = 0; 
        end 
        
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
%OriSFdata.alloricurve = alloricurve;
%OriSFdata.allangles = allangles;
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
% OriSFdata.rRFcenter = row_rf_center_contra; % the center does not change comparing to the reference RF
% OriSFdata.cRFcenter = col_rf_center_contra;
OriSFdata.zero_polarity = zero_polarity;

% RF_data.rfAllPrefEye = rfAllPrefEye; 
RF_data.allCXrfON = rfAllPrefEyeON; 
RF_data.allCXrfOFF = rfAllPrefEyeOFF;

    if debug 
        % plot_receptive_field(OriSFdata,RF_data,electrode_position,reference_ori_map)
        plot_rf_eye_dominance(OriSFdata,RF_data,ODCrtxPlt_smooth,electrode_position,eyepref)
        
        % debug ori selectivity vs onoff balance 
        cv = OriSFdata.CV_map(electrode_position,:); 
        onoff_bal = OriSFdata.CxOnOffBal(electrode_position,:); 
        figure, plot(onoff_bal, cv,'ko'); lsline;
        [r,p]=corr(onoff_bal', cv');
        title(['R = ', num2str(r), ' p = ', num2str(p)]);
        axis square, xlabel('ONOFF balance'), ylabel('CV(OS)'), 
        set(gca,'box' ,'off','TickDir','OUT')
    end 

end


%% function

function CxOnOffBal = ONOFF_balance(CXRFON,CXRFOFF)
    tmpCXRFON_norm = max(max(abs(CXRFON )));
    tmpCXRFOFF_norm = max(max(abs(CXRFOFF)));
    CxOnOffBal =  (tmpCXRFON_norm - tmpCXRFOFF_norm) /(tmpCXRFON_norm + tmpCXRFOFF_norm);
    CxOnOffBal = abs(CxOnOffBal); 
end
 
function plot_receptive_field(OriSFdata,cortex_RF,electrodePos,reference_ori_map)
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

ori_pref_rf = OriSFdata.ori_pref_rf;

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
    title(sprintf('%.2f \n %.2f',ori_pref_rf(electrodePos,ii),reference_ori_map(electrodePos,ii)))
    
    axes('Position',[xs(ii/2+1),0.00,width,height])
    SFPlot = SFHistRFReference{electrodePos,ii};
    %SFPlot = (SFHistRFReference(electrodePos,ii,:));
    plot(SFPlot,'lineWidth',2);
    axis square
    set(gca, 'XScale', 'log');
    set(gca,'box','off','Tickdir','out','YTickLabel',[],'XTickLabel',[])
    title(sprintf('%.2f',SF50Raw(electrodePos,ii)))
    
    annotation('textbox',[0.01 0.88 0.1 0.1],'String',['row = ' num2str(electrodePos)],'EdgeColor','none','fontsize',14)
end

colormap('jet')
end 

function plot_rf_eye_dominance(OriSFdata,cortex_RF,ODCrtxPlt_smooth,electrodePos,eyepref)

% ODI = OriSFdata.ODI; 
LHIRaw = OriSFdata.LHI_map;
CVRaw = OriSFdata.CV_map;
LPI_map = OriSFdata.LPI_map;
    
allCXrfONReference = cortex_RF.allCXrfON ;
allCXrfOFFReference = cortex_RF.allCXrfOFF ;

OriHistRFReference = OriSFdata.OriTunPolarPlot;
OriBinRFReference = OriSFdata.angles;

SFHistRFReference = OriSFdata.allsfcurve_dog;
sfpref = OriSFdata.SF_map; 
sf_bin_deg = OriSFdata.sf_bin_deg; 
SF50Raw = OriSFdata.SF50_map;

CxOnOffBal = OriSFdata.CxOnOffBal;

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
%     title(sprintf('cv:%.2f \n lhi:%.2f \n od:%.2f ',...
%         CVRaw(electrodePos,ii),LHIRaw(electrodePos,ii),ODI(electrodePos,ii)))
    title(sprintf('cv:%.2f \n lhi:%.2f \n bal:%.2f ',...
        CVRaw(electrodePos,ii),LHIRaw(electrodePos,ii),CxOnOffBal(electrodePos,ii)))
    
    axes('Position',[xs(ii/2+1),0.1,width,height])
    SFPlot = SFHistRFReference{electrodePos,ii};
    %SFPlot = (SFHistRFReference(electrodePos,ii,:));
    plot(SFPlot,'lineWidth',2);
    axis square
    set(gca, 'XScale', 'log');
    set(gca,'box','off','Tickdir','out','YTickLabel',[],'XTickLabel',[])
    title(sprintf('sf:%.2f \n 50:%.2f\n lpi:%.2f',sfpref(electrodePos,ii), SF50Raw(electrodePos,ii),LPI_map(electrodePos,ii)))
    
    annotation('textbox',[0.01 0.88 0.1 0.1],'String',[eyepref ',row = ' num2str(electrodePos)],'EdgeColor','none','fontsize',14)
end

colormap('jet')
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
