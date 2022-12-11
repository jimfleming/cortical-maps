function [appdata] = show_maps_segregation(app, appdata)


ONOFF_smoothed = appdata.ONOFF_smoothed;
ODCrtxPlt_smooth = appdata.ODCrtxPlt_smooth;
% ONOFF_interpolated = appdata.ONOFF_interpolated;
% ODCrtxPlt_interpolated = appdata.ODCrtxPlt_interpolated;

% map1 = [133 94 82;255 153 153;33 31 133;153 153 255]./255;
onoff_od_colormap = [255 100 100 ; 255 0 0    ; 0 175 255 ; 0 0 255]./255;
appdata.map1 = onoff_od_colormap;
% appdata.onoff_od_colormap = onoff_od_colormap;

width_cortex = size(ODCrtxPlt_smooth, 1); 
%% figure 
axes_od = app.UIAxes_Segregation_OD ; 
axes_onoff = app.UIAxes_Segregation_ONOFF ; 
axes_od_onoff = app.UIAxes_OD_ONOFF ; 

% axes_od.Visible = 'on';
cla(axes_od), imagesc(ODCrtxPlt_smooth>.5, 'Parent', axes_od)
colormap(axes_od, 'gray'), axis(axes_od, 'square')
title(axes_od, 'Segregation OD', 'fontsize', 14)
xlim(axes_od, [1 width_cortex]), 
ylim(axes_od, [1 width_cortex])
% axis(axes_od, 'off')
setAxes(axes_od)

axes_onoff.Visible = 'on';
cla(axes_onoff), imagesc(ONOFF_smoothed>.5, 'Parent', axes_onoff)
colormap(axes_onoff, 'jet'), axis(axes_onoff, 'square')
title(axes_onoff, 'Segregation ONOFF', 'fontsize', 14)
xlim(axes_onoff, [1 width_cortex]), ylim(axes_onoff, [1 width_cortex])
%axis(axes_onoff, 'off')
setAxes(axes_onoff)

axes_od_onoff.Visible = 'on';
z3 = (double((ODCrtxPlt_smooth>0.5))+1)*1 ;  
z4 = (double(ONOFF_smoothed<.5)+0)*2 ;
imagesc(z3+z4, 'Parent', axes_od_onoff)
colormap(axes_od_onoff, onoff_od_colormap);  axis(axes_od_onoff, 'square'), 
hold(axes_od_onoff, 'on')
title(axes_od_onoff, 'OD ONOFF','fontsize',14)
hold(axes_od_onoff, 'on'), contour(axes_od_onoff, ODCrtxPlt_smooth, app.od_contour_levels, 'k', 'LineWidth', 1);
xlim(axes_od_onoff, [1 width_cortex]), ylim(axes_od, [1 width_cortex])
%axis(axes_od_onoff, 'off')
setAxes(axes_od_onoff)
        
end

