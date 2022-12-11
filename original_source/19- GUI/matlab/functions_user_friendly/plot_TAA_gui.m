function plot_TAA_gui(appdata, crtx_row, crtx_col, ax_arbor, ax_TAA, ax_color, eye_input)
    
ii = crtx_row; 
jj = crtx_col; 
ONOFF_smoothed = appdata.ONOFF_smoothed;

n_interp = 2; 
if strcmp(eye_input, 'Contra')
    TAA = appdata.data_primord_contra.TAA{ii, jj};
elseif strcmp(eye_input, 'Ipsi')
    TAA = appdata.data_primord_ipsi.TAA{ii, jj};
end 
TAA_resized = interp2(TAA, n_interp); 
TAA_resized = imgaussfilt(TAA_resized,2);

%% to plot spread 
sdspread     = appdata.rAffSpread;
ODCrtxPlt    = appdata.ODCrtxPlt;
ONOFFCrtxPlt = appdata.ONOFFCrtxPlt;
ODCrtxPlt_smooth = appdata.ODCrtxPlt_smooth;

% ODMap = ODCrtxPlt == ODCrtxPlt(ii, jj);
if strcmp(eye_input,'Contra')
    ODMap = ODCrtxPlt == 1;
elseif strcmp(eye_input,'Ipsi')
    ODMap = ODCrtxPlt == -1;
end
    
SelectionCondition = find_selected_afferents(ODMap, ii, jj, sdspread); 
plot_Arbor_spread(ODCrtxPlt, ONOFFCrtxPlt, SelectionCondition, sdspread, crtx_row, crtx_col, ODCrtxPlt_smooth, ax_arbor, ax_color, eye_input)
%% 
% pol_map = (ONOFF_smoothed<.5) * -1 + (ONOFF_smoothed>=.5) * 1; 
% pol = pol_map(ii, jj); 
pol = ONOFFCrtxPlt(ii, jj);
contour_level = 1; 
imagesc(ax_TAA, TAA_resized * pol),
axis(ax_TAA, 'square'),
colormap(ax_TAA, 'jet'), 
% caxis(ax_TAA, [-1 1])
caxis(ax_TAA, [-1.25 1.25])

hold(ax_TAA, 'on'),
contour(ax_TAA, TAA_resized, contour_level, 'k', 'LineWidth', 1);

ax_TAA.XColor = ax_color;
ax_TAA.YColor = ax_color;
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



    


 
