function [appdata] = smooth_map_app(appdata)
    % Self: Blurs
    ODCrtxPlt = appdata.ODCrtxPlt;
    ONOFFCrtxPlt = appdata.ONOFFCrtxPlt;
    n_interp = appdata.n_interpol;

    % OD/ONOFF smoothing

    kern_std = 2;   % The std of kernel for smoothing the OD map

    kern = fspecial('gaussian', kern_std*6, kern_std);

    if appdata.aff_sampling_density < 1 
        kern = fspecial('gaussian', 1, 1);
    end 

    ODCrtxPlt_smooth = conv2(ODCrtxPlt, kern, 'same');
    ODCrtxPlt_smooth = (ODCrtxPlt_smooth - min(min(ODCrtxPlt_smooth))) ./ (max(max(ODCrtxPlt_smooth)) - min(min(ODCrtxPlt_smooth)));
    ODCrtxPlt_interpolated = imresize(ODCrtxPlt_smooth, n_interp);

    ONOFF_smoothed = conv2(ONOFFCrtxPlt, kern, 'same');
    ONOFF_smoothed = (ONOFF_smoothed - min(min(ONOFF_smoothed))) ./ (max(max(ONOFF_smoothed)) - min(min(ONOFF_smoothed)));
    ONOFF_interpolated = imresize(ONOFF_smoothed, n_interp);

    appdata.ONOFF_smoothed = ONOFF_smoothed;
    appdata.ODCrtxPlt_smooth = ODCrtxPlt_smooth;
    appdata.ONOFF_interpolated = ONOFF_interpolated;
    appdata.ODCrtxPlt_interpolated = ODCrtxPlt_interpolated;
end
