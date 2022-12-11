function [appdata] = smooth_map_app(appdata)

ODCrtxPlt = appdata.ODCrtxPlt;
ONOFFCrtxPlt = appdata.ONOFFCrtxPlt;
n_interp = appdata.n_interpol;

% OD/ONOFF smoothing


kern_std = 2;   % The std of kernel for smoothing the OD map
%n_interp = 1;   % The number of interpolation for OD map (OD and ori maps must have same size, use 0 for non-smoothed ori map and use 2 for smoothed ori map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kern = fspecial('gaussian', kern_std*6, kern_std);

if appdata.aff_sampling_density < 1 
    kern = fspecial('gaussian', 1, 1);
end 

ODCrtxPlt_smooth = conv2(ODCrtxPlt, kern, 'same');
ODCrtxPlt_smooth = (ODCrtxPlt_smooth - min(min(ODCrtxPlt_smooth))) ./ (max(max(ODCrtxPlt_smooth)) - min(min(ODCrtxPlt_smooth)));
%   ODCrtxPlt_plot = interp2(ODCrtxPlt_smooth, n_interp);
ODCrtxPlt_interpolated = imresize(ODCrtxPlt_smooth, n_interp);

ONOFF_smoothed = conv2(ONOFFCrtxPlt, kern, 'same');
ONOFF_smoothed = (ONOFF_smoothed - min(min(ONOFF_smoothed))) ./ (max(max(ONOFF_smoothed)) - min(min(ONOFF_smoothed)));
%     ONOFF_plot = interp2(ONOFF_smoothed, n_interp);
ONOFF_interpolated = imresize(ONOFF_smoothed, n_interp);

appdata.ONOFF_smoothed = ONOFF_smoothed;
appdata.ODCrtxPlt_smooth = ODCrtxPlt_smooth;
appdata.ONOFF_interpolated = ONOFF_interpolated;
appdata.ODCrtxPlt_interpolated = ODCrtxPlt_interpolated;

% appdata.map1 = [133 94 82;255 153 153;33 31 133;153 153 255]./255;


end