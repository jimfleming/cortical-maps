function output_ori_analysis = analysis_ori_map(reference_ori_map,ODCrtxPlt_smooth, ...
                                ODCrtxPlt_interpolated,pix2mm,n_interpol,debug,show_fig)
    % Orientation map geometrical analysis 
    [ori_map_smooth,ori_map_interpolated]  = smooth_ori_euler(reference_ori_map,ODCrtxPlt_smooth,n_interpol,debug,show_fig);

    % Ori OD contour intersection 
    [intersections, intersection_angles] = od_ori_intersection2(ODCrtxPlt_smooth,ori_map_interpolated);

    % Pinwheel location and charge detection (Kai Lun)
    [pw_locations, pw_charges] = get_pinwheel_location_and_charge(ori_map_interpolated,1);

    % pw distance to OD borders 
    pix2mm_interpol = pix2mm / n_interpol; 
    dist_pw_od = calculate_dist_pw_od(pw_locations,ODCrtxPlt_interpolated,pix2mm_interpol,show_fig); 

    output_ori_analysis.ori_map_smooth = ori_map_smooth;
    output_ori_analysis.ori_map_interpolated = ori_map_interpolated;
    output_ori_analysis.intersections = intersections;
    output_ori_analysis.intersection_angles = intersection_angles;
    output_ori_analysis.pw_locations = pw_locations;
    output_ori_analysis.pw_charges = pw_charges;
    output_ori_analysis.dist_pw_od = dist_pw_od;
end 
