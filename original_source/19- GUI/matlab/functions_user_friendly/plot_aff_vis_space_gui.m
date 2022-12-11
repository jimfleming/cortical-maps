function plot_aff_vis_space_gui(appdata, RFreONOFFnorm, crtx_row, crtx_col, ax_visapce, ax_cortex, ax_color, eye_input)
    
% plot RF mature, afferents with rfs (shown with circle) in visual space,
% selected afferents in cortical space, arbor spread in cortical space,
% orientation map with OD ONOFF contours, tuning curve

ii = crtx_row; 
jj = crtx_col; 

% cortex_RF = rf_contra_mature; 

sdspread      = appdata.rAffSpread;
idxAffRfspace = appdata.IndRetinaRfspace;
ret_rf_rad    = appdata.RFONCenter; % rf_sd_RGC; 

ODCrtxPlt    = appdata.ODCrtxPlt;
ONOFFCrtxPlt = appdata.ONOFFCrtxPlt;
ODCrtxPlt_smooth = appdata.ODCrtxPlt_smooth;
RetONOFFsorted = appdata.RetONOFFsorted;

%%
% ODMap = ODCrtxPlt == ODCrtxPlt(ii, jj);
if strcmp(eye_input,'Contra')
    ODMap = ODCrtxPlt == 1;
elseif strcmp(eye_input,'Ipsi')
    ODMap = ODCrtxPlt == -1;
end
SelectionCondition = find_selected_afferents(ODMap, ii, jj, sdspread); 
% RFreONOFFnorm = find_receptive_field_fig5S(cortex_RF, ii, jj); 

aff_center_weight_plot = find_afferents_weights(SelectionCondition, RetONOFFsorted, ONOFFCrtxPlt, idxAffRfspace, RFreONOFFnorm);

%% plot 

plot_Arbor_spread(ODCrtxPlt, ONOFFCrtxPlt, SelectionCondition, sdspread, crtx_row, crtx_col, ODCrtxPlt_smooth, ax_cortex, ax_color, eye_input)
plot_afferents_vis_space(aff_center_weight_plot, RFreONOFFnorm, ret_rf_rad, ax_visapce, ax_color)

end
%% functions 
function SelectionCondition = find_selected_afferents(ODMap, ii, jj, sdspread)
    [col_grid, row_grid] = meshgrid(1:size(ODMap, 1));
    mask_spread = imdialte_grid(row_grid, col_grid, ii, jj, sdspread);
    SelectionCondition = ODMap == 1 & mask_spread == 1;
end
   
function mask_spread = imdialte_grid(row_grid, col_grid, ii, jj, sdspread)
% replacing imdilate function of Matlab (to make the code faster)
% the function uses the meshgrid with the number of rows and cols of input image 

%%%%%%%%%%%%%%%%%%%%%%%
    dist_matrix = sqrt((col_grid - jj).^2 + (row_grid - ii).^2);
    mask_spread = dist_matrix <= sdspread; 

end
%%         - Figure, Receptive Field Reference
function RFreONOFFnorm = find_receptive_field_fig5S(cortex_RF, ii, jj)

    RFrefON = cortex_RF.allCXrfON{ii, jj} ;
    RFrefOFF = cortex_RF.allCXrfOFF{ii, jj} ;

    RFreONOFF = RFrefON + RFrefOFF;
    MaxValrf = max(max(abs(RFrefON(:))),max(abs(RFrefOFF(:))));
    RFreONOFFnorm = RFreONOFF / MaxValrf;
    
end

function aff_center_weight_plot = find_afferents_weights(SelectionCondition, RetONOFFsorted, ONOFFCrtxPlt, idxAffRfspace, RFreONOFFnorm)

[inda, indb] = find( SelectionCondition ); 
aff_center_weight_plot = nan(length(inda),4);

weightsRFspaceON  = RFreONOFFnorm  .* (RFreONOFFnorm>0);
weightsRFspaceOFF = RFreONOFFnorm  .* (RFreONOFFnorm<0);

for pp = 1:length(inda)
    indRf = RetONOFFsorted(inda(pp),indb(pp));
   % addon = RetinaRF{indRf};
    [rCirc, cCirc] = find(indRf == idxAffRfspace);
    
    if ONOFFCrtxPlt(inda(pp), indb(pp)) == 1
        tmpweights = weightsRFspaceON(rCirc, cCirc) ;
%         allCXrfTempON = allCXrfTempON +  tmpweights * addon;
    elseif ONOFFCrtxPlt(inda(pp), indb(pp)) == -1
        tmpweights = abs(weightsRFspaceOFF(rCirc,cCirc)) ;
%         allCXrfTempOFF =  allCXrfTempOFF + tmpweights * addon; %   figure,imagesc(allCXrfTempOFF),caxis([-1 1])
    end
    
    aff_center_weight_plot(pp,1) = cCirc;
    aff_center_weight_plot(pp,2) = rCirc;
    aff_center_weight_plot(pp,3) = tmpweights;
    aff_center_weight_plot(pp,4) = ONOFFCrtxPlt(inda(pp),indb(pp)); 
end
        
end 


function plot_afferents_vis_space(aff_center_weight_plot, SingleRFNorm, ret_rf_rad, ax, axis_color)

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
    

    %% Cortical RF with afferents RF in visual space 
    % Afferents Weights in Visual Space 
    % nexttile, 
    imagesc(ax, SingleRFNorm),
    caxis(ax, [-1 1]),
    axis(ax, 'square'),
    colormap(ax, 'jet'),
    % axis(ax, 'off')
    % title(ax, 'RF (Visual Space)')
    
    ax.XColor = axis_color;
    ax.YColor = axis_color;
    
    % ON afferents in Visual Space 
    hold(ax, 'on')
    radiCirc = ret_rf_rad;
    num_on = 0;
    MeanON = mean(weightON(~isnan(weightON))); 
    for qq1 = 1 : length(rfONXcenter)
        if ~isnan(weightON(qq1))
            if weightON(qq1) > MeanON  %* .5
                %if temp_mask1(round(rfONYcenter(qq1)), round(rfONXcenter(qq1)))
                    % viscircles(ax, [rfONXcenter(qq1), rfONYcenter(qq1)], radiCirc, 'Color', 'r', 'Linewidth', weightON(qq1)*1);
                    appviscircles(ax, [rfONXcenter(qq1), rfONYcenter(qq1)], radiCirc, 'Color', 'r', 'Linewidth', weightON(qq1)*1);

                %end 
                num_on = num_on + 1;
            end 
        end
    end

    % OFF afferents in Visual Space 
    hold(ax, 'on')
    num_off = 0;
    MeanOFF = mean(weightOFF(~isnan(weightOFF))); 
    for qq2 = 1 : length(rfOFFXcenter)
        if ~isnan(weightOFF(qq2))
            if weightOFF(qq2) >  MeanOFF  %* .5
                %if temp_mask1(round(rfOFFYcenter(qq2)), round(rfOFFXcenter(qq2)))
                 % viscircles(ax, [rfOFFXcenter(qq2), rfOFFYcenter(qq2)], radiCirc, 'Color', 'b', 'Linewidth', weightOFF(qq2)*1);
               appviscircles(ax, [rfOFFXcenter(qq2), rfOFFYcenter(qq2)], radiCirc, 'Color', 'b', 'Linewidth', weightOFF(qq2)*1);
                 %end 
                num_off = num_off + 1;
            end 
        end
    end

end 

function plot_orimap_figS5(OriSFdata, ODCrtxPlt_smooth, ONOFF_smoothed, crtx_row, crtx_col, sdspread)
    contour_level_show =1; 
     nexttile, imagesc(OriSFdata.ori_map_smooth), colormap(gca, 'hsv') , axis('square'), colorbar, axis off
     hold on, contour(ODCrtxPlt_smooth, contour_level_show, 'k', 'LineWidth', 2);
     hold on, contour(ONOFF_smoothed, contour_level_show, 'k', 'LineWidth', 2);
     hold on, viscircles([crtx_col, crtx_row], sdspread,'Color','y');
     hold on, plot(crtx_col, crtx_row, 'yo')
     title('Orientation map')
end

function plot_Arbor_spread(ODCrtxPlt, ONOFFCrtxPlt, SelectionCondition, sdspread, crtx_row, crtx_col, ODCrtxPlt_smooth, ax, axis_color, eye_input)
%     nexttile, %imagesc( ONOFFCrtxPlt ), axis square, colormap(gca, 'jet'),
%     map1 = [255 100 100; 255 0 0; 0 175 255; 0 0 255]./255;
    z1 = (double((ODCrtxPlt>0))+1)*1 ;
    z2 = (double(ONOFFCrtxPlt==-1)+0)*2 ;
%     imagesc(z1+z2)
%     colormap(gca, map1); axis square;
%     caxis([1 4])
%     axis off
%     hold on, contour(ODCrtxPlt_smooth, 1,'Color','k', 'linewidth', 2);
%     hold on, viscircles([crtx_col, crtx_row], sdspread,'Color','y');
%     hold on, plot(crtx_col, crtx_row, 'yo')
%     title('Arbor Spread (Cortical Space)')

    %% Selected Afferents in Cortex
    % nexttile, % imagesc(SelectionCondition .* ONOFFCrtxPlt), axis square, colormap(gca, 'jet'), axis off
    if strcmp(eye_input, 'Contra') %  ODCrtxPlt(crtx_row, crtx_col) > 0
        map2 = [131 255 123; 255 0 0 ; 0 0 255 ]./255;
    elseif strcmp(eye_input, 'Ipsi')  % ODCrtxPlt(crtx_row, crtx_col) < 0
        map2 = [131 255 123; 255 100 100; 0 175 255]./255;
    end
    imagesc(ax, SelectionCondition .* (z1+z2)), axis(ax, 'square'), colormap(ax, map2),
    % axis(ax, 'on')
    hold(ax, 'on'), plot(ax, crtx_col, crtx_row, 'ko')
    % title(ax, sprintf('Selected Aff (Cortical Space) \n Total number of Aff : %.0f', sum(SelectionCondition(:))) )
    
     
    ax.XColor = axis_color;
    ax.YColor = axis_color;
end



    


 
