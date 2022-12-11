function [RF_data,OriSFdata] = generate_rf_cortex_mature_initial(eyepref,rf_contra_primord,data_contra_primord,rf_ipsi_primord,data_ipsi_primord,reference_ori_map,RetinaRF,rowRange,colRange,...
    idxAffRfspace,sdspread,RetONOFFsorted,ONOFFCrtxPlt,ODCrtxPlt,ODCrtxPlt_smooth,electrode_position,sigma_LHI,n_ori_smooth,sf_lSamp,pix2deg,weight_coef_RF,alpha,debug)
%   modifying receptive field for contra/ipsi eye based on the reference
%   orientaiton map 
%   this function rotates the rf to match the preferred ori based on the mature orientation map

%   eyepref           : 'contra'/'ipsi'; % the eye for that new rf will be calculated
%   reference_ori_map : the reference orientation map 
%   alpha             : for iso retinotopic measurement : percentage of maximum reponse in RF space to look for from the center of RF
%
%   Sohrab Note              : The old name of this function was rf_modify_eye5_test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
rf_on_contra = rf_contra_primord.allCXrfON ;
rf_off_contra = rf_contra_primord.allCXrfOFF ;
row_rf_center_contra = data_contra_primord.rRFcenter;
col_rf_center_contra = data_contra_primord.cRFcenter ;
pref_orimap_contra = data_contra_primord.OriPreferred ;
% ODI_contra = data_contra_primord.ODI;

rf_on_ipsi = rf_ipsi_primord.allCXrfON ;
rf_off_ipsi = rf_ipsi_primord.allCXrfOFF ;
row_rf_center_ipsi = data_ipsi_primord.rRFcenter;
col_rf_center_ipsi = data_ipsi_primord.cRFcenter;
pref_orimap_ipsi = data_ipsi_primord.OriPreferred;
% ODI_ipsi = data_ipsi_primord.ODI;

%%
%   clear rfAllPrefEye
od = zeros(size(RetONOFFsorted));
Ori_map = zeros(size(od));
CV_map = zeros(size(od));
SF_map = zeros(size(Ori_map));
SF50_map = zeros(size(Ori_map));
LPI_map = zeros(size(Ori_map));

% allsfcurve = zeros([size(od) length(sf_lSamp)]);
% alloricurve = zeros([size(od) length(sf_lSamp)]);
% allangles = zeros([size(od) length(sf_lSamp)]);
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
R_iso_retinotopy = zeros(size(od));
CxOnOffBal = zeros(size(od));
zero_polarity = zeros(size(od));

%% Main 
dilate_str = strel('disk',round(sdspread),0);
 
if strcmp(eyepref,'contra')
    ODMap = ODCrtxPlt == 1;
elseif strcmp(eyepref,'ipsi')
    ODMap = ODCrtxPlt == -1;
end
 
%%
% for supplementary figure 1 
% rowRange = 30 ; 
% colRange = 46 ; 
% debug = 1 ; 
for ii = rowRange 
    for jj = colRange 
        
        rf_power = 2; % adding power makes orientation map more similar to reference orientation map 
        % weights of affs are assigned based on dominant eye 
        if ODCrtxPlt(ii,jj) == 1
            rfON  = rf_on_contra{ii,jj} .^ rf_power;
            rfOFF = abs(rf_off_contra{ii,jj} .^ rf_power) * -1;
            rowRFcenter = round(row_rf_center_contra(ii,jj));
            colRFcenter = round(col_rf_center_contra(ii,jj));
            oriPref = pref_orimap_contra(ii,jj) ;            
        elseif ODCrtxPlt(ii,jj) == -1
            rfON  = rf_on_ipsi{ii,jj} .^ rf_power;
            rfOFF = abs(rf_off_ipsi{ii,jj} .^ rf_power) * -1;
            rowRFcenter = round(row_rf_center_ipsi(ii,jj));
            colRFcenter = round(col_rf_center_ipsi(ii,jj));
            oriPref = pref_orimap_ipsi(ii,jj) ;            
        end
        
        %% adjusting the weights of RFs of primodial map based on reference orientation map

        rfONOFF = rfON + rfOFF;
        % finding the rf new weights by rotating to match the desired orientation
        rfRotationAngle = -oriPref + reference_ori_map(ii,jj) ;
        
        if mod(size(rfONOFF,1),2) == 0
            rfONOFF2 = rfONOFF(1:end-1,1:end-1);
        else
            rfONOFF2 = rfONOFF;
        end
        rowCenterRotation = floor (size(rfONOFF2,1) / 2);
        colCenterRotation = floor (size(rfONOFF2,2) / 2);
        
        rowShift = rowCenterRotation - rowRFcenter;
        colShift = colCenterRotation - colRFcenter;
        
        rfONOFFshift = circshift(rfONOFF, [rowShift , colShift]);
        rfONOFFshiftRotated = imrotate(rfONOFFshift,rfRotationAngle,'crop'); % counter clockwise
        rfONOFFshiftback = circshift(rfONOFFshiftRotated,[ -rowShift , -colShift ]);
        weightsRFspaceON = rfONOFFshiftback  .* (rfONOFFshiftback>0);
        weightsRFspaceOFF = rfONOFFshiftback .* (rfONOFFshiftback<0);

        % using rfONOFF rather than ON/OFF separately as the weight for the new RFs, showed better resutls for both ori/sf
%         rfONshift = circshift(rfON, [rowShift , colShift]);
%         rfONshiftRotated = imrotate(rfONshift,rfRotationAngle,'crop'); % counter clockwise
%         rfONshiftback = circshift(rfONshiftRotated,[ -rowShift , -colShift ]);
%         rfOFFshift = circshift(rfOFF, [rowShift , colShift]);
%         rfOFFshiftRotated = imrotate(rfOFFshift,rfRotationAngle,'crop'); % counter clockwise
%         rfOFFshiftback = circshift(rfOFFshiftRotated,[ -rowShift , -colShift ]);
%         weightsRFspaceON = rfONshiftback;
%         weightsRFspaceOFF = rfOFFshiftback;
      
        %% applying the new weights to the available afferents 
        %SelectionCondition = SelectCellCortexSpread(sdspread,ODMap,ii,jj);
        [SelectionCondition,~] = SelectCellCortexSpread2(dilate_str,ODMap,ii,jj);
        [inda, indb] = find( SelectionCondition );

        aff_center_weight_plot = nan(length(inda),4); 
        
        allCXrfTempON = zeros(size(idxAffRfspace));
        allCXrfTempOFF = zeros(size(idxAffRfspace));
        for pp = 1:length(inda)
            indRf = RetONOFFsorted(inda(pp),indb(pp));
            addon = RetinaRF{indRf};
            [rCirc,cCirc] = find(indRf == idxAffRfspace);

            if ONOFFCrtxPlt(inda(pp),indb(pp)) == 1
                tmpweights = weightsRFspaceON(rCirc,cCirc) ;
                allCXrfTempON = allCXrfTempON +  tmpweights * addon; 
            elseif ONOFFCrtxPlt(inda(pp),indb(pp)) == -1
                tmpweights = abs(weightsRFspaceOFF(rCirc,cCirc)) ;
                allCXrfTempOFF =  allCXrfTempOFF + tmpweights * addon; %   figure,imagesc(allCXrfTempOFF),caxis([-1 1])
            end
            
            if debug
               
                aff_center_weight_plot(pp,1) = cCirc;
                aff_center_weight_plot(pp,2) = rCirc;
                aff_center_weight_plot(pp,3) = tmpweights;
                aff_center_weight_plot(pp,4) = ONOFFCrtxPlt(inda(pp),indb(pp));
            end
        end
        
        %% Modifying synaptic weights based on pw locations  
        max_resp_on = sum(allCXrfTempON(:));        % max(allCXrfTempON(:));
        max_resp_off = sum(abs(allCXrfTempOFF(:))); % max(abs(allCXrfTempOFF(:)));

        if max_resp_on ~= 0 && max_resp_off ~= 0
            weight_rf = weight_coef_RF(ii,jj); 
            if ONOFFCrtxPlt(ii,jj) == 1
                LWON =  1;
                LWOFF =  weight_rf * (max_resp_on/max_resp_off);
            elseif ONOFFCrtxPlt(ii,jj) == -1
                LWON  =  weight_rf * (max_resp_off/max_resp_on);
                LWOFF =  1;
            end
        else 
            LWON = 1; 
            LWOFF =  1;
            zero_polarity(ii,jj) = 1;  % to measure the number of cortical locations that non-dominant polarity is zero
        end 
             
        allCXrfTempON1 = allCXrfTempON   * LWON; 
        allCXrfTempOFF1 = allCXrfTempOFF * LWOFF; 
        
        %% tuning curve measurements 
        RFSpaceSim =   allCXrfTempON1 + allCXrfTempOFF1 ;
        SingleRFNorm = RFSpaceSim / max(abs(RFSpaceSim(:)));
        [ori_tuning_resp,sf_tuning_resp,sf_tuning_resp_dog,cv,ori_pref,sf_pref_deg,sf50_deg,LPI,angles_bin,num_cycles_bin_deg] = ... 
               fft_ori_sf_tuning6(SingleRFNorm,sf_lSamp,0,pix2deg,0); 
        
       if debug
           % debug_RF(SingleRFNorm,angles_bin,ori_tuning_resp,ori_pref,sf_tuning_resp_dog,sf_lSamp,sf_pref_deg,pix2deg)
%            allCXrfTempON_norm = allCXrfTempON / max(abs(RFSpaceSim(:)));
%            allCXrfTempOFF_norm = allCXrfTempOFF / max(abs(RFSpaceSim(:)));
           debug_RF2(allCXrfTempON,allCXrfTempOFF,SingleRFNorm,aff_center_weight_plot,angles_bin,ori_tuning_resp,ori_pref,num_cycles_bin_deg,sf_tuning_resp_dog,sf_pref_deg,sf50_deg, ONOFFCrtxPlt, SelectionCondition, sdspread, ii, jj)
           text = [eyepref ', row : ' num2str(ii) ', col : ' num2str(jj)] ;
           annotation('textbox',[0.03 0.87 0.98 0.08],'String',text,'EdgeColor','none','fontsize',12)
       end
            
        %%
        max_on_response(ii,jj) = max(allCXrfTempON1(:)); 
        max_off_response(ii,jj) = max(abs(allCXrfTempOFF1(:)));         
        CxOnOffBal(ii,jj) = ONOFF_balance(allCXrfTempON1,allCXrfTempOFF1); 
        
        % iso retinotopic (patch alpha = .5) % 50 of the maximum value is set to find the 
        radius_rf = measure_radius_rf(rowRFcenter,colRFcenter,allCXrfTempON1,allCXrfTempOFF1,alpha); 
        
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

OriSFdata.num_projected_aff = num_projected_aff; 
OriSFdata.max_on_response = max_on_response; 
OriSFdata.max_off_response = max_off_response;
OriSFdata.CxOnOffBal = CxOnOffBal;
OriSFdata.R_iso_retinotopy = R_iso_retinotopy; 
OriSFdata.rRFcenter = row_rf_center_contra; % the center does not change comparing to the reference RF
OriSFdata.cRFcenter = col_rf_center_contra;
OriSFdata.zero_polarity = zero_polarity;

% RF_data.rfAllPrefEye = rfAllPrefEye; 
RF_data.allCXrfON = rfAllPrefEyeON; 
RF_data.allCXrfOFF = rfAllPrefEyeOFF;

if debug
    % plot_receptive_field(OriSFdata,RF_data,electrode_position,reference_ori_map)
    plot_rf_eye_dominance(OriSFdata,RF_data,ODCrtxPlt_smooth,electrode_position,eyepref)
end

end


%% function

function [SelectionCondition,mask_spread] = SelectCellCortexSpread2(dilate_str,ODMap,rAffCrtx,cAffCrtx)
    % selecting neighboring cells in cortex based on spread value
    % the spread function is a Gaussian distribution
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dilate_str = strel('disk',round(sdspread),0);
    mask_spread = zeros(size(ODMap));
    mask_spread(rAffCrtx,cAffCrtx) = 1;
    mask_spread = imdilate(mask_spread,dilate_str); %   figure,imagesc(mask_spread)

    % Selecting those cells in retina that have similar organization to cortex
     SelectionCondition =  ODMap == 1 & mask_spread == 1 ;
end

function CxOnOffBal = ONOFF_balance(CXRFON,CXRFOFF)
    tmpCXRFON_norm = max(max(abs(CXRFON )));
    tmpCXRFOFF_norm = max(max(abs(CXRFOFF)));
    CxOnOffBal =  (tmpCXRFON_norm - tmpCXRFOFF_norm) /(tmpCXRFON_norm + tmpCXRFOFF_norm);
end

function r_fit_circle = measure_radius_rf(rowRFcenter,colRFcenter,allCXrfTempON,allCXrfTempOFF,alpha)
%  measure the radius from the center of RF to a percentage (alpha) of the
%  maximum response for iso-retinotopic patch measurement
   % rf_abs = abs(allCXrfTempON) + abs(allCXrfTempOFF);
   rf_abs = abs(allCXrfTempON + allCXrfTempOFF);
   rf_thresh = max(rf_abs(:)) * alpha ;
   rf_diff = abs(rf_abs - rf_thresh);% figure, imagesc(rf_diff), colorbar
   valid_data1 = rf_diff < .1; % the difference with the set alpha is less than 10%
    
   not_valid_area = rf_abs < .05 ; % removing the least 5% value of all the map 

   valid_data = valid_data1 .* (~not_valid_area); 
   [point_row,point_col] = find(valid_data); 
   dist_center_rf = sqrt((point_row(:)-rowRFcenter).^2 + (point_col(:)-colRFcenter).^2); 
   
   r_not_valid = measure_r_not_valid(allCXrfTempON,allCXrfTempOFF,rowRFcenter,colRFcenter); 
   dist_center_rf( dist_center_rf < r_not_valid) = []; 
   r_fit_circle = round(mean(dist_center_rf));
   
   if isnan(r_fit_circle)
       r_fit_circle = r_not_valid; 
   end 
       
   if 1 == 0
       figure,
       rf_norm = allCXrfTempON + allCXrfTempOFF; 
       rf_norm = rf_norm / max(abs(rf_norm(:))); 
       subplot(121),imagesc(rf_norm), colormap(gca,'jet'), axis square
       subplot(122),imagesc(rf_abs), hold on, plot(colRFcenter,rowRFcenter,'ko'), axis square
       hold on
       viscircles([colRFcenter,rowRFcenter],r_fit_circle)
   end
end 

function   r_remove = measure_r_not_valid(allCXrfTempON,allCXrfTempOFF,rowRFcenter,colRFcenter)
% removing the range for RF from the center of RF and center of each polarity 
% 
% r_remove : distance of center of RF and max center of ON or OFF RF separately
    
   if  max(abs(allCXrfTempON(:))) > 0    
       [r_on,c_on] = find(max(abs(allCXrfTempON(:))) == allCXrfTempON); 
       dist_on = norm([rowRFcenter-r_on,colRFcenter-c_on]); 
   else 
       dist_on = 0; 
   end 
   
   if  max(abs(allCXrfTempOFF(:))) > 0 
       [r_off,c_off] = find(max(abs(allCXrfTempOFF(:)))== allCXrfTempOFF); 
       dist_off = norm([rowRFcenter-r_off,colRFcenter-c_off]);
   else 
       dist_off = 0; 
   end 
   
   if dist_on>dist_off
       lim_dist = dist_on; 
   else 
       lim_dist = dist_off; 
   end 
   r_remove = round(lim_dist);
   
%    temp1 = zeros(size(rf_abs));
%    temp1(rowRFcenter,colRFcenter) = 1;
%    disk_filt = strel('disk',r_remove,0);
%    not_valid_area_in = imdilate(temp1,disk_filt); % figure, imagesc(not_valid_area)
%    
%    not_valid_area_out = rf_abs < .05 ; % removing the least 5% value of all the map 
%    not_valid_area = not_valid_area_in + not_valid_area_out ; 
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
    title(sprintf('cv:%.2f \n lhi:%.2f \n ',...
        CVRaw(electrodePos,ii),LHIRaw(electrodePos,ii)))
    
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

function debug_RF2(allCXrfTempON,allCXrfTempOFF,SingleRFNorm,aff_center_weight_plot,angles_bin,ori_tuning_resp,ori_pref,sf_bin,sf_tuning_curve,sf_pref_deg,sf50_deg, ONOFFCrtxPlt, SelectionCondition, sdspread, ii, jj)
    %debug
    ind_nan = aff_center_weight_plot(:,3) == 0; 
    aff_center_weight_plot(ind_nan,3) = nan; 
    
    data = aff_center_weight_plot;
    ind_on = data(:,4) == 1;
    rfONXcenter = data(ind_on,1);
    rfONYcenter = data(ind_on,2);
    weightON = data(ind_on,3);

    ind_off = ~ind_on;
    rfOFFXcenter = data(ind_off,1);
    rfOFFYcenter = data(ind_off,2);
    weightOFF = data(ind_off,3);
    %%
    figure,clf
    %RFSpaceSim = allCXrfTempON + allCXrfTempOFF ;
    %SingleRFNorm = RFSpaceSim / max(abs(RFSpaceSim(:)));
    caxisVal = 1;%0.6;
    rf_on_plot = (SingleRFNorm > 0) .* SingleRFNorm ; 
    rf_off_plot = (SingleRFNorm < 0) .* SingleRFNorm ; 
    subplot(241),imagesc(rf_on_plot),caxis([-caxisVal caxisVal]),axis square,colormap('jet'),colorbar
    subplot(242),imagesc(rf_off_plot),caxis([-caxisVal caxisVal]),axis square,colormap('jet'),colorbar
    subplot(243),imagesc(SingleRFNorm),caxis([-1 1]),axis square,colormap('jet'),colorbar
    
    subplot(245), cla, imagesc(ONOFFCrtxPlt),axis square,colormap('jet'), hold on, viscircles([jj, ii],sdspread)
    title('cortical space')
    subplot(246),imagesc(ONOFFCrtxPlt .* SelectionCondition),caxis([-1 1]),axis square,colormap('jet')
    title('cortical space')
    %title(sprintf('R aff covered:%.2f',r_covered_aff))
    
    %hold on, plot(rfcenONOFFcol,rfcenONOFFrow,'ko')
    %subplot(242),imagesc(SelectionCondition.*ONOFFCrtxPlt),axis square,colormap('jet'),colorbar
    subplot(244),cla, imagesc(zeros(size(SingleRFNorm))),set(gca,'ydir','reverse'),axis 'square'; colorbar
    hold on , scatter(rfONXcenter,rfONYcenter,weightON*1,'r','filled')
    hold on , scatter(rfOFFXcenter,rfOFFYcenter,weightOFF*1,'b','filled')
   % title(sprintf('dist pol:%.2f',dist_pol))
     title('visual space (aff weights)')
    
    %%
    subplot(247), polar1(angles_bin,ori_tuning_resp','k');title(sprintf('Ori = %.2f',ori_pref))
    subplot(248),plot(sf_bin,sf_tuning_curve,'lineWidth',2); axis square
    set(gca, 'XScale', 'log','box','off','Tickdir','out','YTickLabel',[],'XTickLabel',[])
    title(sprintf('sf:%.2f \n sf50:%.2f',sf_pref_deg,sf50_deg ))
end