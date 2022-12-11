function fig1_supp_rf_eyes(eyepref,rf_input,rowRange,colRange,idxAffRfspace,sdspread,RetONOFFsorted,ONOFFCrtxPlt,ODCrtxPlt,sf_lSamp,pix2deg,weight_coef_RF,debug)

%%
rf_on = rf_input.allCXrfON ;
rf_off = rf_input.allCXrfOFF ;

% row_rf_center = data_eye.rRFcenter;
% col_rf_center = data_eye.cRFcenter ;

%% Main 
dilate_str = strel('disk',round(sdspread),0);
 
if strcmp(eyepref,'contra')
    ODMap = ODCrtxPlt == 1;
elseif strcmp(eyepref,'ipsi')
    ODMap = ODCrtxPlt == -1;
end
 
%%
for ii = rowRange 
    for jj = colRange 
        
        allCXrfTempON = rf_on{ii, jj}; 
        allCXrfTempOFF = rf_off{ii, jj}; 
        
        weightsRFspaceON = allCXrfTempON  .* (allCXrfTempON>0);
        weightsRFspaceOFF = allCXrfTempOFF .* (allCXrfTempOFF<0);
        
        %% applying the new weights to the available afferents 
        %SelectionCondition = SelectCellCortexSpread(sdspread,ODMap,ii,jj);
        [SelectionCondition,~] = SelectCellCortexSpread2(dilate_str,ODMap,ii,jj);
        [inda, indb] = find( SelectionCondition );

        aff_center_weight_plot = nan(length(inda),4); 
        
        for pp = 1:length(inda)
            indRf = RetONOFFsorted(inda(pp),indb(pp));
            [rCirc,cCirc] = find(indRf == idxAffRfspace);

            if ONOFFCrtxPlt(inda(pp),indb(pp)) == 1
                tmpweights = weightsRFspaceON(rCirc,cCirc) ;
            elseif ONOFFCrtxPlt(inda(pp),indb(pp)) == -1
                tmpweights = abs(weightsRFspaceOFF(rCirc,cCirc)) ;
            end
            
            if debug
                aff_center_weight_plot(pp,1) = cCirc;
                aff_center_weight_plot(pp,2) = rCirc;
                aff_center_weight_plot(pp,3) = tmpweights;
                aff_center_weight_plot(pp,4) = ONOFFCrtxPlt(inda(pp),indb(pp));
            end
        end

        %% tuning curve measurements 
        
        RFSpaceSim =   allCXrfTempON + allCXrfTempOFF ;
        SingleRFNorm = RFSpaceSim / max(abs(RFSpaceSim(:)));
        [ori_tuning_resp,sf_tuning_resp,sf_tuning_resp_dog,cv,ori_pref,sf_pref_deg,sf50_deg,LPI,angles_bin,num_cycles_bin_deg] = ... 
               fft_ori_sf_tuning6(SingleRFNorm,sf_lSamp,0,pix2deg,0); 
        
       if debug
           debug_RF2(allCXrfTempON,allCXrfTempOFF,SingleRFNorm,aff_center_weight_plot,angles_bin,ori_tuning_resp,ori_pref,num_cycles_bin_deg,sf_tuning_resp_dog,sf_pref_deg,sf50_deg, ONOFFCrtxPlt, SelectionCondition, sdspread, ii, jj)
           text = [eyepref ', row : ' num2str(ii) ', col : ' num2str(jj)] ;
           annotation('textbox',[0.03 0.87 0.98 0.08],'String',text,'EdgeColor','none','fontsize',12)
       end
            
    end
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

function [LWON, LWOFF, zero_polarity ] = modify_synaptic_weight(allCXrfTempON, allCXrfTempOFF, weight_coef_RF, ONOFFCrtxPlt)
    % Modifying synaptic weights based on pw locations
    max_resp_on = sum(allCXrfTempON(:));        % max(allCXrfTempON(:));
    max_resp_off = sum(abs(allCXrfTempOFF(:))); % max(abs(allCXrfTempOFF(:)));

    if max_resp_on ~= 0 && max_resp_off ~= 0
        weight_rf = weight_coef_RF;
        if ONOFFCrtxPlt == 1
            LWON =  1;
            LWOFF =  weight_rf * (max_resp_on/max_resp_off);
        elseif ONOFFCrtxPlt == -1
            LWON  =  weight_rf * (max_resp_off/max_resp_on);
            LWOFF =  1;
        end
    else
        LWON = 1;
        LWOFF =  1;
        zero_polarity = 1;  % to measure the number of cortical locations that non-dominant polarity is zero
    end
end 


%% debug
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
    hold on , scatter(rfONXcenter,rfONYcenter,weightON*.5,'r','filled')
    hold on , scatter(rfOFFXcenter,rfOFFYcenter,weightOFF*.5,'b','filled')
   % title(sprintf('dist pol:%.2f',dist_pol))
     title('visual space (aff weights)')
    
    %%
    subplot(247), polar1(angles_bin,ori_tuning_resp','k');title(sprintf('Ori = %.2f',ori_pref))
    subplot(248),plot(sf_bin,sf_tuning_curve,'lineWidth',2); axis square
    set(gca, 'XScale', 'log','box','off','Tickdir','out','YTickLabel',[],'XTickLabel',[])
    title(sprintf('sf:%.2f \n sf50:%.2f',sf_pref_deg,sf50_deg ))
end
