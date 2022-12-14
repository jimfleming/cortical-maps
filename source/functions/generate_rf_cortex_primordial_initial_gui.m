function [output_rf,output_data] = generate_rf_cortex_primordial_initial3_2_TAA_uf(pref_eye,RetinaRF,r_covered_aff,rowRange,colRange,ODCrtxPlt,ONOFFCrtxPlt, ... 
    RetONOFFsorted,sdspread,i_aff_rf_space,synaptic_weight_factor,ODCrtxPlt_smooth,pix2deg,sf_lSamp,n_ori_smooth,sigma_LHI_2d,debug,show_fig)
% The spread Gaussian function is applied based on the distance in Cortrex not RF space 
% Measure receptive field in cortex after afferent spread in cortical plate 
% Cortical RFs : Gaussian weighted sum of afferents 
%
% the RF is measured both for one eye (ipsi or contra)
%   
% 
% r_covered_aff   % Average radius of area covered by an afferent in visual space (diameter ~= 1 deg)
%
% Sohrab Note       :     the old name was "make_rf_cortex_same_eye2"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug2 = 0; 
%% Initialization (pre assigning the matrices)
    allCXrfnorm{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];
    allCXrfON{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];
    allCXrfOFF{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];
    rRFcenter = zeros(size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2));
    cRFcenter = zeros(size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2));
%     rRFcenterON = zeros(size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2));
%     cRFcenterON = zeros(size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2));
%     rRFcenterOFF = zeros(size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2));
%     cRFcenterOFF = zeros(size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2));

    OriPreferred = zeros(size(ONOFFCrtxPlt));
    dist_center_on_off = zeros(size(ONOFFCrtxPlt));
    SFPref = zeros(size(ONOFFCrtxPlt));
    
    max_on_response = zeros(size(ONOFFCrtxPlt)); 
    max_off_response = zeros(size(ONOFFCrtxPlt));
    max_on_response_norm = zeros(size(ONOFFCrtxPlt));
    max_off_response_norm = zeros(size(ONOFFCrtxPlt));

    CxOnOffBal = zeros(size(ONOFFCrtxPlt));
    SF50_pref = zeros(size(ONOFFCrtxPlt));
    LPI = zeros(size(ONOFFCrtxPlt));
    CV = zeros(size(ONOFFCrtxPlt));
    ODI = zeros(size(ONOFFCrtxPlt));

    OriHist{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];
    %OriBin{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];
    sf_tuning_all{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];
    sf_tuning_dog_all{size(ONOFFCrtxPlt,1),size(ONOFFCrtxPlt,2)} = [];

    val_rf_dom_at_r_coverage = zeros(size(ONOFFCrtxPlt));
    rf_dom_pow2 = zeros(size(ONOFFCrtxPlt));

    TAA = initialize_TAA(ONOFFCrtxPlt, rowRange, colRange);
%% Main 
    mask_covered_aff  = strel('disk', r_covered_aff, 0);
    mask_covered_aff2 = strel('disk', r_covered_aff-1 , 0);
    %dilate_str = strel('disk',round(sdspread) ,0);

    weight_arbor_spread = make_cortical_spread_dist(sdspread); 
    weight_coef_RF = synaptic_weight_factor / max(synaptic_weight_factor(:)); % 0 to 1 Cycle per degree
    
    if strcmp(pref_eye,'contra')
        ODMap = ODCrtxPlt == 1;
    elseif strcmp(pref_eye,'ipsi')
        ODMap = ODCrtxPlt == -1;
    end
    
    [col_grid, row_grid] = meshgrid(1:size(ODMap, 1)); 
    width_rf_space = size(RetinaRF{1}, 1); 
    [col_grid_rf, row_grid_rf] = meshgrid(1:width_rf_space); 
    %%
    for ii =  rowRange
        for jj = colRange 
            
            WTemp = zeros(size(ODCrtxPlt));
            WTemp(ii,jj) = 1; 
            tmpweights = filter2(weight_arbor_spread,WTemp);   
               
            mask_spread = imdialte_grid(row_grid, col_grid, ii, jj, sdspread); 
%             dist_matrix = sqrt((col_grid - jj).^2 + (row_grid - ii).^2); 
%             mask_spread = dist_matrix <= sdspread; %   figure,imagesc(mask_spread)
            SelectionCondition = ODMap == 1 & mask_spread == 1;%   figure,imagesc(SelectionCondition)
            %[SelectionCondition,~] = SelectCellCortexSpread2(dilate_str, ODMap, ii, jj);
            [inda, indb] = find( SelectionCondition );           
            
            % [ODI(ii,jj),ODI_contra(ii,jj),ODI_ipsi(ii,jj)] = measure_ODI(SelectionCondition,tmpweights,mask_spread,ODCrtxPlt(ii,jj)); 
            ODI(ii,jj) = sum(sum(SelectionCondition.*tmpweights));
            
            %% dominant RF
            aff_center_weight_plot = nan(length(inda),4); 
            dominant_pol = ONOFFCrtxPlt(ii,jj);
            non_dominant_pol = -dominant_pol;
            [allCXrfTempDominant, aff_center_weight_plot, TAA] = gen_rf(RetinaRF, RetONOFFsorted, ONOFFCrtxPlt, i_aff_rf_space, inda, indb, tmpweights, dominant_pol, 1, aff_center_weight_plot, TAA, ii, jj, debug);
            
            % this method worked better to find the cneter of RFs this center in mature_rf code is used to rotate this RF and find new weight for afferents
            rf_norm = allCXrfTempDominant / max(abs(allCXrfTempDominant(:)));
            
            %  center of rfs = max resp of dominant polarity ( reason : rotation of non_dom rf around the max response  )
%             [center_rf_row,center_rf_col] = find(abs(rf_norm) == 1);
            
            % finding the center of rfs by averaging all the pixels 
            [center_rf_row, center_rf_col] = find_center_rf(rf_norm); 
            if isnan(center_rf_row)
                figure, 
                subplot(121), imagesc(ONOFFCrtxPlt), axis square, hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1)
                hold on, viscircles([jj, ii], sdspread, 'color', 'k')
                subplot(122), imagesc(SelectionCondition .* ONOFFCrtxPlt), axis square
                hold on, viscircles([jj, ii], sdspread, 'color', 'k')
                colormap jet 
                text = ['There is no Afferent reaching at this cortical location from ' pref_eye ' eye']; 
                sgtitle(text)
                error(text)
            end 
            % check
%             [rf_norm2,val_rf_dom_at_r_coverage(ii,jj),rf_dom_pow2(ii,jj)] = add_power_rf_dominant(rf_norm, mask_covered_aff, mask_covered_aff2, center_rf_row, center_rf_col, .1); % add explanation in method 
            [rf_norm2, val_rf_dom_at_r_coverage(ii,jj), rf_dom_pow2(ii,jj)] = add_power_rf_dominant2(rf_norm, row_grid_rf, col_grid_rf, r_covered_aff, center_rf_row, center_rf_col, .1); % add explanation in method 
            
            % another methods that I used 
            % in case there are separate subregions 
%             Feature = abs(rf_norm) > .1;
%             Feature_props = regionprops(Feature,'centroid','Area');
%             [~,n_features] = bwlabel(Feature);
%             areas = zeros(n_features,1);
%             for qq = 1 : n_features
%                 areas(qq) = Feature_props(qq).Area ;
%             end
%             [~,ind] = max(areas);
%             center_rf_col = round(Feature_props(ind).Centroid(1));
%             center_rf_row = round(Feature_props(ind).Centroid(2));
 
            % finding the center of RF by calculating the average of weighted row/col
%             Feature = abs(rf_norm) > .1;
%             size_rf = size(rf_norm,1); 
%             [col_grid,row_grid] = meshgrid(1:size_rf); 
%             row_weighted = row_grid(Feature(:)) .* rf_norm(Feature(:)); 
%             row_mean = sum(row_weighted(:)) / sum(rf_norm(Feature(:))) ; 
%             center_rf_row = round(row_mean); 
%             col_weighted = col_grid(Feature(:)) .* rf_norm(Feature(:)); 
%             col_mean = sum(col_weighted(:)) / sum(rf_norm(Feature(:))) ; 
%             center_rf_col = round(col_mean); 

            %% non_dominant RF 
%             temp_mask1 = zeros(size(rf_norm2));
%             temp_mask1(center_rf_row, center_rf_col) = 1;
%             temp_mask1 = imdilate(temp_mask1, mask_covered_aff);
            temp_mask1 = imdialte_grid(row_grid_rf, col_grid_rf, center_rf_row, center_rf_col, r_covered_aff); 
            
            weight_other_pol = (1 - abs(rf_norm2 .^ 1)) .* temp_mask1;
            [allCXrfTempNonDominant, aff_center_weight_plot, TAA] = gen_rf(RetinaRF, RetONOFFsorted, ONOFFCrtxPlt, i_aff_rf_space, inda, indb, weight_other_pol, non_dominant_pol, 0, aff_center_weight_plot, TAA, ii, jj, debug);

            if dominant_pol == 1
                allCXrfTempON = allCXrfTempDominant;
                allCXrfTempOFF = allCXrfTempNonDominant;
            elseif dominant_pol == -1
                allCXrfTempOFF = allCXrfTempDominant;
                allCXrfTempON = allCXrfTempNonDominant;
            end
            
            [row_center_on,col_center_on] =  find(max(allCXrfTempON(:)) == allCXrfTempON);
            [row_center_off,col_center_off] =  find(max(abs(allCXrfTempOFF(:))) == abs(allCXrfTempOFF));
             
            % in case there are two separate subregions 
            try
                dist_pol = norm(row_center_on-row_center_off, col_center_on-col_center_off);
                dist_center_on_off(ii,jj) = dist_pol;
            catch
                dist_center_on_off(ii,jj) = nan; 
            end 
            
            %% Modifying the weight and tuning curve measurement 
            [LWON, LWOFF ] = modify_synaptic_weight(allCXrfTempON, allCXrfTempOFF, weight_coef_RF(ii,jj), ONOFFCrtxPlt(ii,jj)); 

            %% 
            CXrfON = allCXrfTempON * LWON;
            CXrfOFF = allCXrfTempOFF * LWOFF;
            CxOnOffBal(ii,jj) = ONOFF_balance(CXrfON,CXrfOFF);
            
            RFSpaceSim = CXrfON + CXrfOFF ;
            SingleRFNorm = RFSpaceSim / max(abs(RFSpaceSim(:)));
            
            max_on_response(ii,jj) = max(allCXrfTempON(:));
            max_off_response(ii,jj) = max(abs(allCXrfTempOFF(:)));
            
            [Max_on, Max_off] = find_max_resp_norm(SingleRFNorm);
            max_on_response_norm(ii,jj) = Max_on;
            max_off_response_norm(ii,jj) = Max_off;
        
%             [ori_tuning_resp, sf_tuning_resp, sf_tuning_resp_dog, cv, ori_pref, sf_pref_deg, sf50_deg, lpi, angles_bin, sf_bin_deg] = ... 
%                fft_ori_sf_tuning6(SingleRFNorm,sf_lSamp,0,pix2deg,0); 
            [ori_tuning_resp, sf_tuning_resp, sf_tuning_resp_dog, cv, ori_pref, sf_pref_deg, sf50_deg, lpi, angles_bin, sf_bin_deg] = ... 
               fft_ori_tuning_uf(SingleRFNorm,sf_lSamp); 
           
            if debug2 == 1 
                % debug_RF(SingleRFNorm,angles_bin,ori_tuning_resp,ori_pref,sf_tuning_resp_dog,sf_lSamp,sf_pref_deg,pix2deg)
                 debug_RF_primord(allCXrfTempON,allCXrfTempOFF,SingleRFNorm,aff_center_weight_plot,angles_bin,ori_tuning_resp, ori_pref, sf_bin_deg, sf_tuning_resp_dog, sf_pref_deg, sf50_deg, dist_pol, weight_other_pol, r_covered_aff)
            end 

            % using the center of dominant polarity (in making mature RFs, it might not work because after rotation, there might not be enough aff available for the other polarity)
            rfcenONOFFrow = center_rf_row(1);
            rfcenONOFFcol = center_rf_col(1); 
            % using the center of abs(rf)
            % [rfcenONOFFrow,rfcenONOFFcol] = measure_center_rf(allCXrfTempON,allCXrfTempOFF,SingleRFNorm); 
            
            %%
            OriPreferred(ii,jj) = ori_pref;
            SFPref(ii,jj) = sf_pref_deg;
            SF50_pref(ii,jj) = sf50_deg;
            LPI(ii,jj) = lpi;
            CV(ii,jj) = cv;
            
            allCXrfnorm{ii,jj} = RFSpaceSim; 
            allCXrfON{ii,jj} = CXrfON; 
            allCXrfOFF{ii,jj} = CXrfOFF; 
            
            rRFcenter(ii,jj) = rfcenONOFFrow; 
            cRFcenter(ii,jj) = rfcenONOFFcol; 
            
            OriHist{ii,jj} = ori_tuning_resp;
            % OriBin{ii,jj} = angles_bin;
            % SFHist{ii,jj} = sf_tuning_resp;
            sf_tuning_all{ii,jj} = sf_tuning_resp;
            sf_tuning_dog_all{ii,jj} = sf_tuning_resp_dog;
            
        end
    end
    
    %% ori smooth  
    [ori_map_smooth,ori_map_interpolated] = smooth_ori_euler(OriPreferred,ODCrtxPlt_smooth,n_ori_smooth,0,0); 

    % LHI measurement 
    LHI_map = measure_LHI(OriPreferred,sigma_LHI_2d); 
    LHI_map_smooth = measure_LHI(ori_map_smooth,sigma_LHI_2d); 
    %%  output
    %output_rf.allCXrfnorm = allCXrfnorm;
    output_rf.allCXrfON = allCXrfON;
    output_rf.allCXrfOFF = allCXrfOFF;
    
    output_data.rRFcenter = rRFcenter;
    output_data.cRFcenter = cRFcenter;
    output_data.CxOnOffBal = CxOnOffBal; 
    
    output_data.max_on_response = max_on_response; 
    output_data.max_off_response = max_off_response; 
    output_data.max_on_response_norm = max_on_response_norm; 
    output_data.max_off_response_norm = max_off_response_norm; 
    
    
    output_data.OriPreferred = OriPreferred;
    output_data.ori_map_smooth = ori_map_smooth;
    output_data.ori_map_interpolated = ori_map_interpolated;
    output_data.SFPref = SFPref ;
    output_data.SF50Pref = SF50_pref;
    output_data.LPI = LPI ;
    output_data.CV = CV;
    output_data.LHI_map = LHI_map;
    output_data.LHI_map_smooth = LHI_map_smooth;
    
    ODI = ODI / max(ODI(:)); 
    output_data.ODI = ODI; 
  
    output_data.OriHist = OriHist;
    output_data.OriBin = angles_bin;
    output_data.sf_bin_deg = sf_bin_deg; 
    output_data.SFHist = sf_tuning_all;
    output_data.SFHist_dog = sf_tuning_dog_all; 
    
    output_data.dist_center_on_off = dist_center_on_off;
    
    output_data.val_rf_dom_at_r_coverage = val_rf_dom_at_r_coverage; 
    output_data.rf_dom_pow2 = rf_dom_pow2; 
    
    output_data.TAA = TAA;
%%  debug  
    if debug 
        % show_onoff_balance(CxOnOffBal,ODCrtxPlt_smooth)
        plot_receptive_field(output_rf,output_data,ii)
        % plot_rf_linear_lhi_test(output_data,ODCrtxPlt_smooth,ii,180,'reference')
    end 
end


%% functions 
function w = make_cortical_spread_dist(sdspread)
    r2 = ceil(sdspread*3);
    mu = [0 0]; %   Center of mask
    rx1 = r2-1;
    x1 = -rx1:rx1;
    x2 = -rx1:rx1;
    [X1,X2] = meshgrid(x1,-x2);
    covariance1 = [sdspread 0;0 sdspread];
    w = mvnpdf([X1(:) X2(:)],mu,covariance1);
    w =  reshape(w,length(x2),length(x1));
    w = w/max(w(:));
    w = w .^ 1;%3
end

function [allCXrf, aff_center_weight_plot, TAA] = gen_rf(RetinaRF, RetONOFFsorted, ONOFFCrtxPlt, i_aff_rf_space, row_ind, col_ind, weights, pol, ...
                                                    is_dominant, aff_center_weight_plot, TAA, ii, jj, debug)
    allCXrf = zeros(size(i_aff_rf_space));
    for pp = 1:length(row_ind)
        indRf = RetONOFFsorted(row_ind(pp), col_ind(pp));
        addon = RetinaRF{indRf};

        % using onoff_smooth instead of onoff may result in some off aff in rf_on and vice versa so just to avoid this condition the following line is added
        sgn_rf_retina = sign(sum(addon(:))); 

        if ONOFFCrtxPlt(row_ind(pp), col_ind(pp)) == pol && sgn_rf_retina == pol
            if is_dominant
                tmp_w = weights(row_ind(pp), col_ind(pp));
                allCXrf = allCXrf + tmp_w .* addon;
            else
                [rCirc, cCirc] = find(indRf == i_aff_rf_space);
                tmp_w = weights(rCirc, cCirc);
                allCXrf = allCXrf + tmp_w * addon;   % check!
            end
            
            % TAA{row_ind(pp), col_ind(pp)}(ii, jj) = tmp_w; 
            
            if debug == 1
                [rCirc, cCirc] = find(indRf == i_aff_rf_space);
                aff_center_weight_plot(pp,1) = cCirc;
                aff_center_weight_plot(pp,2) = rCirc;
                aff_center_weight_plot(pp,3) = tmp_w;
                aff_center_weight_plot(pp,4) = ONOFFCrtxPlt(row_ind(pp), col_ind(pp));
            end
        end
        
        if ONOFFCrtxPlt(row_ind(pp), col_ind(pp)) == pol
        % using onoff_smooth instead of onoff may result in some off aff in rf_on and vice versa so just to avoid this condition the following line is added
        % there are conditions that ONOFFCrtxPlt(row_ind(pp), col_ind(pp)
        % == pol is met but sgn_rf_retina == pol is not true, because
        % smoothing makes some cortical locations to switch the polarity 
            if is_dominant
                tmp_w = weights(row_ind(pp), col_ind(pp));
            else
                [rCirc, cCirc] = find(indRf == i_aff_rf_space);
                tmp_w = weights(rCirc, cCirc);
            end
            TAA{row_ind(pp), col_ind(pp)}(ii, jj) = tmp_w; 
        end 
    end
end



function [ODI,ODI_contra,ODI_ipsi] = measure_ODI(SelectionCondition,tmpweights,mask_spread,eye_pref)
    ODI = sum(sum(SelectionCondition.*tmpweights));
    aff_dominant_eye_weighted =    (mask_spread .* SelectionCondition) .* tmpweights;
    aff_nondominant_eye_weighted  = (mask_spread .* ~SelectionCondition) .* tmpweights;

    if eye_pref == 1
        ODI_contra = sum(aff_dominant_eye_weighted(:));
        ODI_ipsi = sum(aff_nondominant_eye_weighted(:));
    elseif eye_pref == -1
        ODI_contra = sum(aff_nondominant_eye_weighted(:));
        ODI_ipsi = sum(aff_dominant_eye_weighted(:));
    end
end 
            
function CxOnOffBal = ONOFF_balance(CXRFON,CXRFOFF)
    tmpCXRFON_norm = max(max(abs(CXRFON )));
    tmpCXRFOFF_norm = max(max(abs(CXRFOFF)));
    CxOnOffBal =  (tmpCXRFON_norm - tmpCXRFOFF_norm) /(tmpCXRFON_norm + tmpCXRFOFF_norm);
end

function [rfcenONOFFrow,rfcenONOFFcol] = measure_center_rf(allCXrfTempON,allCXrfTempOFF,SingleRFNorm)
%     max_on = max(allCXrfTempON(:));
%     [row_on,col_on] = find( max_on == allCXrfTempON);
% 
%     allCXrfTempOFF = abs(allCXrfTempOFF);
%     max_off = max(allCXrfTempOFF(:));
%     [row_off,col_off] = find( max_off == allCXrfTempOFF);
% 
%     rfcenONOFFrow = (row_on * max_on + row_off * max_off) / (max_on + max_off);
%     rfcenONOFFcol = (col_on * max_on + col_off * max_off) / (max_on + max_off);

    rf_abs = abs(SingleRFNorm); 
    [rr,cc] = find(rf_abs > 0.01);
    rrRetinotopy = 0;
    ccRetinotopy = 0;
    Weight = 0; 
    for ii = 1 : length(rr)
        weight = rf_abs(rr(ii),cc(ii)); 
        rrRetinotopy = rrRetinotopy + rr(ii) * weight;
        ccRetinotopy = ccRetinotopy + cc(ii) * weight;
        Weight = Weight + weight; 
    end 
    rfcenONOFFrow = round(rrRetinotopy / Weight);
    rfcenONOFFcol = round(ccRetinotopy / Weight);
end

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

%% debug
function debug_RF_primord(allCXrfTempON,allCXrfTempOFF,SingleRFNorm,aff_center_weight_plot,angles_bin,ori_tuning_resp,ori_pref,sf_bin,sf_tuning_curve,sf_pref_deg,sf50_deg,dist_pol,weight_other_pol,r_covered_aff)
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
    caxisVal = 5;%0.6;
    subplot(241),imagesc(allCXrfTempON),axis square,colormap('jet'),colorbar%,caxis([-caxisVal caxisVal])
    subplot(242),imagesc(allCXrfTempOFF),axis square,colormap('jet'),colorbar%,caxis([-caxisVal caxisVal])
    subplot(243),imagesc(SingleRFNorm),caxis([-1 1]),axis square,colormap('jet'),colorbar
    
    subplot(245),imagesc(weight_other_pol),caxis([-1 1]),axis square,colormap('jet'),colorbar
    title(sprintf('R aff covered:%.2f',r_covered_aff))
    
    %hold on, plot(rfcenONOFFcol,rfcenONOFFrow,'ko')
    %subplot(242),imagesc(SelectionCondition.*ONOFFCrtxPlt),axis square,colormap('jet'),colorbar
    subplot(244),imagesc(zeros(size(SingleRFNorm))),set(gca,'ydir','reverse'),axis 'square'; colorbar
    hold on , scatter(rfONXcenter,rfONYcenter,weightON*50,'r','filled')
    hold on , scatter(rfOFFXcenter,rfOFFYcenter,weightOFF*50,'b','filled')
    title(sprintf('dist pol:%.2f',dist_pol))
    
    
    %%
    subplot(247), polar1(angles_bin,ori_tuning_resp','k');title(sprintf('Ori = %.2f',ori_pref))
    subplot(248),plot(sf_bin,sf_tuning_curve,'lineWidth',2); axis square
    set(gca, 'XScale', 'log','box','off','Tickdir','out','YTickLabel',[],'XTickLabel',[])
    title(sprintf('sf:%.2f \n sf50:%.2f',sf_pref_deg,sf50_deg ))
end

%%
function plot_receptive_field(output_rf,output_data,electrodePos)

    allCXrfONReference = output_rf.allCXrfON ;
    allCXrfOFFReference = output_rf.allCXrfOFF ;

    OriHistRFReference = output_data.OriHist;
    OriBinRFReference = output_data.OriBin;

    sf_bin = output_data.sf_bin_deg; 
    SFHistRFReference = output_data.SFHist_dog;
    sfpref = output_data.SFPref; 
    sfpref50 = output_data.SF50Pref; 
    %sf_bin_deg = OutputRFCrtxReference.sf_bin_deg; 
    
    dist_center_on_off = output_data.dist_center_on_off; 
    
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

        axes('Position',[xs(ii/2+1),0.25,width,height])
        ori_tuning = OriHistRFReference{electrodePos,ii};
        %angles = OriBinRFReference{electrodePos,ii};%{electrodePos,ii};
        angles = OriBinRFReference; 
        ph = polar1(angles,ori_tuning','k');

    %     r = sfpref(electrodePos,ii);
    %     CyclePix = size(RFreONOFF,1) / r; % Cycle in pixel
    %     CycleDeg = CyclePix / deg2pix;
    %     CyclePerDeg = 1/CycleDeg;

        axes('Position',[xs(ii/2+1),0.00,width,height])
        SFPlot = SFHistRFReference{electrodePos,ii};
        %SFPlot = (SFHistRFReference(electrodePos,ii,:));
        plot(sf_bin,SFPlot,'lineWidth',2);
        axis square
        set(gca, 'XScale', 'log');
        set(gca,'box','off','Tickdir','out','YTickLabel',[],'XTickLabel',[])
        % title(sprintf('%.2f \n %.2f',sfpref(electrodePos,ii),sfpref50(electrodePos,ii)))
        title(sprintf('%2.0f \n %.2f \n %.2f',dist_center_on_off(electrodePos,ii),sfpref(electrodePos,ii),sfpref50(electrodePos,ii)))
        
        annotation('textbox',[0.01 0.88 0.1 0.1],'String',['row = ' num2str(electrodePos)],'EdgeColor','none','fontsize',14)
    end

    colormap('jet')
end 


function [rf_norm2,val_rf_dom_at_r_coverage,rf_dom_pow2] = add_power_rf_dominant(rf_norm, mask_covered_aff, mask_covered_aff2, center_rf_row, center_rf_col, thresh_rf)
% The function calculates the average value of rf_norm at the r_covered_aff
% it adds power function to the rf_norm if it is above the thresh_rf
% to bring the average 
%
% Inputs 
% thresh_rf = .1;   the desired value of the sum of RF at r afferent coeverage value below thresh_rf
% mask_covered_aff  = strel('disk',r_covered_aff,0);
% mask_covered_aff2 = strel('disk',r_covered_aff-1 ,0);
%
% Outputs 
% rf_norm : The new rf_nom, if the value is under the threshold, this does not change 
% val_rf_dom_at_r_coverage : the average value of rf_norm at the r_covered_aff
% rf_dom_pow2              : the power used to calculate the new rf_norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    temp_mask1 = zeros(size(rf_norm));
    temp_mask1(center_rf_row,center_rf_col) = 1;
    mask_coverage_circle = imdilate(temp_mask1, mask_covered_aff) - imdilate(temp_mask1, mask_covered_aff2);
    rf_dom_at_r_coverage = (mask_coverage_circle .* rf_norm) ;
    val_rf_dom_at_r_coverage = abs(sum(rf_dom_at_r_coverage(:)) / sum(mask_coverage_circle(:)));

    rf_dom_pow2 = 0; 
    if val_rf_dom_at_r_coverage > thresh_rf
        rf_dom_pow2 = log(thresh_rf)/log(val_rf_dom_at_r_coverage); % val_rf_dom_at_r_coverage .^ rf_dom_pow2 = thresh_rf ==> rf_dom_pow2 = LOGval_rf_dom_at_r_coverage(thresh_rf)
        rf_norm2 = abs(rf_norm .^ rf_dom_pow2) ;%* sgn_rf_cortex;
        %   allCXrfTempDominant = rf_norm .* max_dominant_rf_value;
        rf_dom_at_r_coverage = (mask_coverage_circle .* rf_norm2) ;
        val_rf_dom_at_r_coverage = sum(rf_dom_at_r_coverage(:)) / sum(mask_coverage_circle(:));
    else 
        rf_norm2 = rf_norm; 
    end
    
end 

function [rf_norm2,val_rf_dom_at_r_coverage,rf_dom_pow2] = add_power_rf_dominant2(rf_norm, row_grid_rf, col_grid_rf, r_covered_aff, center_rf_row, center_rf_col, thresh_rf)
% The function calculates the average value of rf_norm at the r_covered_aff
% it adds power function to the rf_norm if it is above the thresh_rf
% to bring the average 
%
% Inputs 
% thresh_rf = .1;   the desired value of the sum of RF at r afferent coeverage value below thresh_rf
% mask_covered_aff  = strel('disk',r_covered_aff,0);
% mask_covered_aff2 = strel('disk',r_covered_aff-1 ,0);
%
% Outputs 
% rf_norm : The new rf_nom, if the value is under the threshold, this does not change 
% val_rf_dom_at_r_coverage : the average value of rf_norm at the r_covered_aff
% rf_dom_pow2              : the power used to calculate the new rf_norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     temp_mask1 = zeros(size(rf_norm));
%     temp_mask1(center_rf_row,center_rf_col) = 1;
    % mask_coverage_circle = imdilate(temp_mask1, mask_covered_aff) - imdilate(temp_mask1, mask_covered_aff2);
    mask_spread1 = imdialte_grid(row_grid_rf, col_grid_rf, center_rf_row, center_rf_col, r_covered_aff-1); 
    mask_spread2 = imdialte_grid(row_grid_rf, col_grid_rf, center_rf_row, center_rf_col, r_covered_aff); 
    mask_coverage_circle = mask_spread2 - mask_spread1; 
    rf_dom_at_r_coverage = (mask_coverage_circle .* rf_norm) ;
    val_rf_dom_at_r_coverage = abs(sum(rf_dom_at_r_coverage(:)) / sum(mask_coverage_circle(:)));

    rf_dom_pow2 = 0; 
    if val_rf_dom_at_r_coverage > thresh_rf
        rf_dom_pow2 = log(thresh_rf)/log(val_rf_dom_at_r_coverage); % val_rf_dom_at_r_coverage .^ rf_dom_pow2 = thresh_rf ==> rf_dom_pow2 = LOGval_rf_dom_at_r_coverage(thresh_rf)
        rf_norm2 = abs(rf_norm .^ rf_dom_pow2) ;%* sgn_rf_cortex;
        %   allCXrfTempDominant = rf_norm .* max_dominant_rf_value;
        rf_dom_at_r_coverage = (mask_coverage_circle .* rf_norm2) ;
        val_rf_dom_at_r_coverage = sum(rf_dom_at_r_coverage(:)) / sum(mask_coverage_circle(:));
    else 
        rf_norm2 = rf_norm; 
    end
    
end 

function [rf_center_row, rf_center_col] = find_center_rf(rf_norm)
% finding the center of RFs by averaging the rf value for all the pixels 
        rf_abs = abs(rf_norm); 
        rf_size = size( rf_abs , 1) ;
        [ XX , YY ] = meshgrid( 1:rf_size , 1:rf_size ); 
        % measure rf center by averaging the rf_abs
        sum_rf = sum( rf_abs , [1,2]);
        rf_xx = sum(rf_abs .* XX , [1,2]) / sum_rf ; 
        rf_yy = sum(rf_abs .* YY , [1,2]) / sum_rf ;
        
        rf_center_row = round(rf_yy); 
        rf_center_col = round(rf_xx); 
end 

function TAA = initialize_TAA(ONOFFCrtxPlt, rowRange, colRange)
    TAA{size(ONOFFCrtxPlt,1), size(ONOFFCrtxPlt,2)} = [];
    for ii =  rowRange
        for jj = colRange 
            TAA{ii, jj} = zeros(size(ONOFFCrtxPlt,1), size(ONOFFCrtxPlt,2)); 
        end 
    end
end 

