function outputRF = measure_ori_primordial_map(RetinaRF,rowRange,colRange,ODCrtxPlt,ONOFFCrtxPlt,RetinaDist, ... 
    RetONOFFsorted,rAffSpread,i_aff_rf_space,synaptic_weight_factor,ODCrtxPlt_smooth,pix2deg,sf_lSamp,eyepref, ...
    electrode_position,debug,show_fig)
% Measure the primordial orientationm map based on the fft of receptive field in 
% cortex after afferent spread in cortical plate 
%
%
% Initialization Getting Input 
%   RetinaDist =  Average distance of cells in retina (5 pixels distance or 20 microns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization (pre assigning the matrices)
allCXrfnorm{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];
allCXrfON{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];
allCXrfOFF{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];
rRFcenter = zeros(size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2));
cRFcenter = zeros(size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2));
rRFcenterON = zeros(size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2));
cRFcenterON = zeros(size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2));
rRFcenterOFF = zeros(size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2));
cRFcenterOFF = zeros(size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2));

OriPreferred = zeros(size(ONOFFCrtxPlt));
SFPref = zeros(size(ONOFFCrtxPlt));
CxOnOffBal = zeros(size(ONOFFCrtxPlt));

OriHist{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];
OriBin{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];
SFHist{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];
    %% Main
    SFMapInitialnorm = synaptic_weight_factor / max(synaptic_weight_factor(:)); % 0 to 1 Cycle per degree

    for ii =  rowRange
        for jj = colRange
             
            if ODCrtxPlt(ii,jj) == 1 
                ODMap = ODCrtxPlt == 1;
            elseif ODCrtxPlt(ii,jj) == -1 
                ODMap = ODCrtxPlt == -1;
            end 
               
            rAffCrtx = ii;
            cAffCrtx = jj;

            SelectionCondition = SelectCellCortexSpread(rAffSpread,ODMap,rAffCrtx,cAffCrtx);% figure,imagesc(SelectionCondition .* ONOFFMap)
            [inda, indb] = find( SelectionCondition );
            
            % to have smooth retinotopy both eyes are considered to measure the average center of RF
            SelectionCellsRetinotopyON = SelectionCondition .* (ONOFFCrtxPlt == 1);% figure,imagesc(SelectionCellsRetinotopyON)
            [rfcenONrow,rfcenONcol] = CenterRF(SelectionCellsRetinotopyON==1,RetONOFFsorted,i_aff_rf_space);
            SelectionCellsRetinotopyOFF = SelectionCondition .* (ONOFFCrtxPlt == -1);% figure,imagesc(ONOFFSelectedRFspace)
            [rfcenOFFrow,rfcenOFFcol] = CenterRF(SelectionCellsRetinotopyOFF==1,RetONOFFsorted,i_aff_rf_space);
            rfcenONOFFrow = (rfcenONrow * LWON + rfcenOFFrow * LWOFF) / (LWON + LWOFF); 
            rfcenONOFFcol = (rfcenONcol * LWON + rfcenOFFcol * LWOFF) / (LWON + LWOFF); 
    
            if ONOFFCrtxPlt(ii,jj) == 1
                LWON =  1;
                LWOFF = SFMapInitialnorm(ii,jj) * .9; % .9 to make it lower than the dominant polarity
            elseif ONOFFCrtxPlt(ii,jj) == -1
                LWON = SFMapInitialnorm(ii,jj) * .9;
                LWOFF =  1;
            end

            weightsRFspaceON = spreadfunction(rfcenONrow,rfcenONcol,rAffSpread,LWON,RetinaDist,i_aff_rf_space);
            weightsRFspaceOFF = spreadfunction(rfcenOFFrow,rfcenOFFcol,rAffSpread,LWOFF,RetinaDist,i_aff_rf_space);

            allCXrfTempON = zeros(size(i_aff_rf_space));
            allCXrfTempOFF = zeros(size(i_aff_rf_space));
            for pp = 1:length(inda)
                indRf = RetONOFFsorted(inda(pp),indb(pp));
                addon = RetinaRF{indRf};%   figure,imagesc(addon),caxis([-1 1])
                [rCirc,cCirc] = find(indRf == i_aff_rf_space);
                if ONOFFCrtxPlt(inda(pp),indb(pp)) == 1
                    tmpweights = weightsRFspaceON(rCirc,cCirc);
                    allCXrfTempON = allCXrfTempON +  tmpweights * addon; %   figure,imagesc(allCXrfTempON),caxis([-1 1])
                elseif ONOFFCrtxPlt(inda(pp),indb(pp)) == -1
                    tmpweights = weightsRFspaceOFF(rCirc,cCirc);
                    allCXrfTempOFF =  allCXrfTempOFF + tmpweights * addon;
                end
            end
             
            CxOnOffBal(ii,jj) = ONOFF_balance(allCXrfTempON,allCXrfTempOFF);
            
            RFSpaceSim = allCXrfTempON + allCXrfTempOFF ;
            SingleRFNorm = RFSpaceSim / max(abs(RFSpaceSim(:)));
            %[ortmp, sfall, cv, angles_grad,snr,pref,rSFmax,angles,img_fft_abs,xr,yr,snrfft] = FftOriSFTuning2(SingleRFNorm);
            
            [ori_tuning_resp,sf_tuning_resp,~,ori_pref,sf_pref_deg,angles_bin] = fft_ori_sf_tuning(SingleRFNorm,sf_lSamp,pix2deg,debug); 

            OriPreferred(ii,jj) = ori_pref;
            SFPref(ii,jj) = sf_pref_deg;
            
            allCXrfnorm{ii,jj} = RFSpaceSim; 
            allCXrfON{ii,jj} = allCXrfTempON; 
            allCXrfOFF{ii,jj} = allCXrfTempOFF; 
            
            rRFcenter(ii,jj) = rfcenONOFFrow; 
            cRFcenter(ii,jj) = rfcenONOFFcol; 
%             rRFcenterON(ii,jj) = rfcenONrow; 
%             cRFcenterON(ii,jj) = rfcenONcol; 
%             rRFcenterOFF(ii,jj) = rfcenOFFrow; 
%             cRFcenterOFF(ii,jj) = rfcenOFFcol; 
            
            OriHist{ii,jj} = ori_tuning_resp;
            %OriBin{ii,jj} = angles_bin;
            SFHist{ii,jj} = sf_tuning_resp;
        end
    end

    outputRF.rRFcenter = rRFcenter;
    outputRF.cRFcenter = cRFcenter;
%     outputRF.rRFcenterON = rRFcenterON;
%     outputRF.cRFcenterON = cRFcenterON;
%     outputRF.rRFcenterOFF = rRFcenterOFF;
%     outputRF.cRFcenterOFF = cRFcenterOFF;
    outputRF.allCXrfnorm = allCXrfnorm;
    outputRF.allCXrfON = allCXrfON;
    outputRF.allCXrfOFF = allCXrfOFF;
    outputRF.CxOnOffBal = CxOnOffBal; 
    
    outputRF.OriPreferred = OriPreferred;
    outputRF.SFPref = SFPref ;

    outputRF.OriHist = OriHist;
    outputRF.OriBin = angles_bin;
    outputRF.SFHist = SFHist;
    
    if show_fig 
        show_onoff_balance(CxOnOffBal,ODCrtxPlt_smooth)
        
    end 
end

%%
function SelectionCondition = SelectCellCortexSpread(sdspread,ODMap,rAffCrtx,cAffCrtx)
    % selecting neighboring cells in cortex based on spread value
    % the spread function is a Gaussian distribution
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dilate_str = strel('disk',round(sdspread),0);
    mask_spread = zeros(size(ODMap));
    mask_spread(rAffCrtx,cAffCrtx) = 1;
    mask_spread = imdilate(mask_spread,dilate_str); %   figure,imagesc(mask_spread)

    % Selecting those cells in retina that have similar organization to cortex
     SelectionCondition =  ODMap == 1 & mask_spread == 1 ;
end


function weightsRFspace = spreadfunction(affrow,affcol,sdspread,weight,RetinaDist,idxAffRfspace)
    % giving weight to afferents in cortex based on the location in RF space
    %
    % affrow,affcol     :   position of RF center 
    % sdspread          :   spread radius
    % weight            :   different weight for ON and OFF
    % RetinaDist        :   average distance of cells in retina
    % idxAffRfspace     :   index afferents in Receptive field space
    % Retinotopy3mIndex :   Retinotopic map of afferents in cortex
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isnan(affrow)
        r2 = ceil(sdspread*RetinaDist);
        mu = [0 0]; %   Center of mask
        rx1 = r2-1;
        x1 = -rx1:rx1;
        x2 = -rx1:rx1;
        [X1,X2] = meshgrid(x1,-x2);
        covariance1 = [sdspread*RetinaDist 0;0 sdspread*RetinaDist];
        w = mvnpdf([X1(:) X2(:)],mu,covariance1);
        w =  reshape(w,length(x2),length(x1));
        w = w/max(w(:));
        WTemp = zeros(size(idxAffRfspace));
        WTemp(affrow,affcol) = 1;
        weightsRFspace = weight * filter2(w,WTemp);   %  figure,imagesc(weightsRFspace)
    else 
        weightsRFspace = zeros(size(i_aff_rf_space));
    end 
end

function [ref_cen_row,ref_cen_col] = CenterRF(SelectionCondition,RetONOFFsorted,IndRetinaRfspace)
%   Calculating the center of RF by averaging the AFF location of neighboring
%   Cortical cells in RF space 

    neighborcell = RetONOFFsorted(SelectionCondition);
    rrRetinotopy = 0;
    ccRetinotopy = 0;
    for kk = 1 : length(neighborcell)
        [rrRet,ccRet] = find( neighborcell(kk) == IndRetinaRfspace);
        rrRetinotopy = rrRetinotopy + rrRet;
        ccRetinotopy = ccRetinotopy + ccRet;
    end
    ref_cen_row = round(rrRetinotopy / length(neighborcell));
    ref_cen_col = round(ccRetinotopy / length(neighborcell));
end 

function CxOnOffBal = ONOFF_balance(CXRFON,CXRFOFF)
    tmpCXRFON_norm = max(max(abs(CXRFON )));
    tmpCXRFOFF_norm = max(max(abs(CXRFOFF)));
    CxOnOffBal =  (tmpCXRFON_norm - tmpCXRFOFF_norm) /(tmpCXRFON_norm + tmpCXRFOFF_norm);
end

%% plot function 
function show_onoff_balance(CxOnOffBal,ODCrtxPlt_smooth)
    % ONOFF balance
    figure
    set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
    od_contour_levels = 1 ; 
    imagesc(CxOnOffBal), colormap(gca,'jet'), axis('square'), colorbar, title('ONOFF balance','fontsize',20)
    hold on, contour(ODCrtxPlt_smooth, od_contour_levels, 'k', 'LineWidth', 5);
end 