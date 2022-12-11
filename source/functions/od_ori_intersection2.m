function [intersections, angles] = od_ori_intersection2(OD_raw,ori)

    kern_std = 2;   % The std of kernel for smoothing the OD map
    n_interp = 2;   % The number of interpolation for OD map (OD and ori maps must have same size, use 0 for non-smoothed ori map and use 2 for smoothed ori map)
    smooth_OD_border = 1;   % if 1, the border of the OD will be smoothed (unrounded matrix)
    ori_contour_levels = 10;%[0, 60, 120, 180];   % contour levels of orientation. If only one level (e.g. k) is needed, use [k k]. If a number is given, e.g. 10, then 10 contour levels will be used.
    label_contours = 0;   % if 1, the contours will be labeled

%     OD_raw = load(OD_path).ODCrtxPlt;
%     ori = load(ori_path).sumX;   % or .OriMapSmoothed for real value ori map, either complex or real ori map is fine

    [OD, ori] = preview_OD_ori(OD_raw, ori, kern_std, n_interp, smooth_OD_border, ori_contour_levels, label_contours);
    OD = imresize(OD,[size(ori,1) size(ori,2)]); 
    [intersections, angles] = compute_intersect_angle(ori, OD, ori_contour_levels, [0.5 0.5], label_contours);

end 