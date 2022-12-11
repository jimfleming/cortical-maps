function [appdata] = smooth_map_app_sn(app, appdata)

ODCrtxPlt = appdata.ODCrtxPlt;
ONOFFCrtxPlt = appdata.ONOFFCrtxPlt;
n_interp = appdata.n_interpol;

% OD/ONOFF smoothing


kern_std = 2;   % The std of kernel for smoothing the OD map
%n_interp = 1;   % The number of interpolation for OD map (OD and ori maps must have same size, use 0 for non-smoothed ori map and use 2 for smoothed ori map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kern = fspecial('gaussian', kern_std*6, kern_std);
ODCrtxPlt_smooth = conv2(ODCrtxPlt, kern, 'same');
ODCrtxPlt_smooth = (ODCrtxPlt_smooth - min(min(ODCrtxPlt_smooth))) ./ (max(max(ODCrtxPlt_smooth)) - min(min(ODCrtxPlt_smooth)));
%ODCrtxPlt_plot = interp2(ODCrtxPlt_smooth, n_interp);
ODCrtxPlt_interpolated = imresize(ODCrtxPlt_smooth, n_interp);

ONOFF_smoothed = conv2(ONOFFCrtxPlt, kern, 'same');
ONOFF_smoothed = (ONOFF_smoothed - min(min(ONOFF_smoothed))) ./ (max(max(ONOFF_smoothed)) - min(min(ONOFF_smoothed)));
%     ONOFF_plot = interp2(ONOFF_smoothed, n_interp);
ONOFF_interpolated = imresize(ONOFF_smoothed, n_interp);

appdata.ONOFF_smoothed = ONOFF_smoothed;
appdata.ODCrtxPlt_smooth = ODCrtxPlt_smooth;
appdata.ONOFF_interpolated = ONOFF_interpolated;
appdata.ODCrtxPlt_interpolated = ODCrtxPlt_interpolated;

map1 = [133 94 82;255 153 153;33 31 133;153 153 255]./255;
appdata.map1 = map1;

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
colormap(axes_od_onoff, map1);  axis(axes_od_onoff, 'square'), 
hold(axes_od_onoff, 'on')
title(axes_od_onoff, 'OD ONOFF','fontsize',14)
hold(axes_od_onoff, 'on'), contour(axes_od_onoff, ODCrtxPlt_smooth, app.od_contour_levels, 'k', 'LineWidth', 1);
xlim(axes_od_onoff, [1 width_cortex]), ylim(axes_od, [1 width_cortex])
%axis(axes_od_onoff, 'off')
setAxes(axes_od_onoff)
        
end

function setAxes(axes)
    axes.Visible = 'on';
    set(axes, 'tickdir', 'out'); set(axes, 'linewidth', 2); set(axes, 'Box', 'on');
    set(axes, 'xtick', []); set(axes, 'ytick', []); %set(axes, 'square');
    pause(0.001)
end