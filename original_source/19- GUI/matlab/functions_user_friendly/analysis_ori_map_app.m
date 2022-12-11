function analysis_ori_map_app(app, appdata)
          

    reference_ori_map = appdata.data_contra_mature.ori_map_smooth; 
    ODCrtxPlt_smooth = appdata.ODCrtxPlt_smooth; 
    ODCrtxPlt_interpolated = appdata.ODCrtxPlt_interpolated; 
%     pix2mm
     n_interpol = appdata.n_interpol; 
%     debug
%     show_fig
      
    %% 
    % Orientation map geometrical analysis 
    [ori_map_smooth,ori_map_interpolated]  = smooth_ori_euler(reference_ori_map, ODCrtxPlt_smooth, n_interpol, 0, 0);

    % Ori OD contour intersection 
    %[intersections, intersection_angles] = od_ori_intersection2(ODCrtxPlt_smooth, ori_map_interpolated);
    [~, angles_ori_od] = compute_intersect_angle_app( ori_map_interpolated, ODCrtxPlt_interpolated, 10, 1, 0, 3, hsv, 'Ori vs OD', 1, app); 

    %%
    axes_orimap_od = app.UIAxes_orimap_contour ; 
    axes_orimap_od.Visible = 'on'; 
    imagesc(axes_orimap_od, ori_map_interpolated ), 
    axis(axes_orimap_od, 'square'), colormap(axes_orimap_od,'hsv'), 
    title(axes_orimap_od, 'Orientation'), axis(axes_orimap_od,'off')
    hold(axes_orimap_od, 'on'), contour(axes_orimap_od, ODCrtxPlt_interpolated, 1, 'k', 'LineWidth', 1);
    hold(axes_orimap_od, 'on'), contour(axes_orimap_od, ori_map_interpolated, 10, 'k', 'LineWidth', 1);
    
    % Pinwheel location and charge detection (Kai Lun)
    % [pw_locations, pw_charges] = get_pinwheel_location_and_charge(ori_map_interpolated,1);

    % pw distance to OD borders 
%     pix2mm_interpol = pix2mm / n_interpol; 
%     dist_pw_od = calculate_dist_pw_od(pw_locations,ODCrtxPlt_interpolated,pix2mm_interpol, 1); 


end 

function [intersections, angles] = compute_intersect_angle_app(map1, map2, contour1_levels, contour2_levels, label_contours, least_polygon_side_num, color_map, figure_title, show_fig, app)
    %
    % map1: main map, all other parameters are referring to map1 if not stated otherwise
    % map2: another map with same shape as map1, only the contour of this map will be plotted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 5
        label_contours = 0;
    end
    if nargin < 6
        least_polygon_side_num = 3;
    end
    if nargin < 7
        show_fig = 1;
    end
    if any(size(map1) ~= size(map2))
        error(["The shape of map1 (" num2str(size(map1)) ") is different from the shape of map2 (" num2str(size(map2)) ")!"]);
    end
    
    c1 = contourc(map1, contour1_levels);
    c2 = contourc(map2, contour2_levels);
    %c2 = contourc(map2, [0 0]);
    
    [contour_coord1, coord_nums1] = get_contour_coord(c1, least_polygon_side_num);
    [contour_coord2, coord_nums2] = get_contour_coord(c2, least_polygon_side_num);
    num_coord1 = sum(coord_nums1); num_coord2 = sum(coord_nums2);
%     disp(["The map1 has " num2str(num_coord1) " contour points and map2 has " num2str(num_coord2) " contour points."]);
    
    intersections_float = [];
    for i = 1:length(contour_coord2)
        polygon2 = contour_coord2{i};
        poly2_distance = sqrt(sum((polygon2(1:end-1, :)-polygon2(2:end, :)).^2, 2));
        for j = 1:length(contour_coord1)
            polygon1 = contour_coord1{j};
            poly1_distance = sqrt(sum((polygon1(1:end-1, :)-polygon1(2:end, :)).^2, 2));
            thres_distance = 2 * max(max(poly1_distance), max(poly2_distance));
            cross_tmp = get_intersection(polygon1, polygon2, thres_distance);
            if size(cross_tmp)
                intersections_float = [intersections_float; cross_tmp];
            end
        end
    end
    intersections = round(intersections_float);
    [~, idx, ~] = unique(intersections, 'rows');
    intersections = intersections(sort(idx), :);
    angles = get_intersect_angle(map1, map2, intersections);
%     disp(["Number of intersections: " num2str(length(intersections))]);
%     if intersections_float
%         %scatter(intersections(:,1), intersections(:,2), 25, 'w', 'filled');
%         scatter(intersections_float(:,1), intersections_float(:,2), 25, 'w', 'filled');
%     end
    %%
    if show_fig == 1 
        xtick = [0 45 90]; 
        ytick = [0 20]; 
        
%         figure 
%         subplot(3,4,[1,2,5,6]), imagesc(map1);  colormap(color_map); title(figure_title,'fontsize',12), axis off; axis square;
%         hold on, contour(map2, contour2_levels, 'k', 'LineWidth', 1);
% %         caxis([-1 1])
%          caxis( [min(map1(:)) max(map1(:))] )
%         
%         subplot(3,4,[3,4,7,8]),imagesc(map1); hold on; colormap(gca, color_map); axis('square'); axis(gca,'off')
%         contour(map1, contour1_levels, 'k', 'LineWidth', 1);
%         contour(map2, contour2_levels, 'w', 'LineWidth', 1);
%         
% %         contour(map1, contour1_levels, 'k', 'LineWidth', 1,'ShowText','on')
% %         contour(map2, contour2_levels, 'w', 'LineWidth', 1,'ShowText','on')
%         
%         if intersections_float
%             %scatter(intersections(:,1), intersections(:,2), 25, 'w', 'filled');
%             scatter(intersections_float(:,1), intersections_float(:,2), 25, 'w', 'filled');
%         end
% %         caxis([-1 1])
%          caxis( [min(map1(:)) max(map1(:))] )
         
%% 
        axes_ori_od_hist = app.UIAxes_ori_od_hist ; 
        axes_ori_od_hist.Visible = 'on'; 
        axes_ori_od_polar_hist = app.UIAxes_ori_od_polar_hist ;  
        axes_ori_od_polar_hist.Visible = 'on'; 
        
        %subplot(3,4,9),
        cla(axes_ori_od_hist),
        [freqs x_ticks] = hist( angles, 5:10:85);
        total_num = sum(freqs); percents = freqs / total_num * 100;
        bar(axes_ori_od_hist, x_ticks+5, percents, 'w');
        ylabel(axes_ori_od_hist, "Percent"); xlabel(axes_ori_od_hist, 'Angle (deg)');
        axis(axes_ori_od_hist, 'square')
        set(axes_ori_od_hist,'XTick',xtick,'YTick',ytick,'XTickLabel',xtick,'YTickLabel',ytick,'Box','off','TickDir','out'); 
        
        
        %subplot(3,4,10)
        %[freqs2, ~] = hist(angles, 2.5:5:87.5);
        freqs_polar = make_freq_polar_format(freqs);
        ori_bin_deg =  5 : 10 : 365;
        %ori_bin_deg2 =  2.5 : 5 : 357.5;
        ori_bin_rad = ori_bin_deg * (pi/180);
        polar1(axes_ori_od_polar_hist, ori_bin_rad, freqs_polar', 'k');
    end 
end

function [contour_coord_arr, element_nums] = get_contour_coord(contour_matrix, least_polygon_side_num)
    done = 0;
    cur_idx = 1;
    level_idx = 1;
    element_nums = [];
    while done == 0
        num_elements = contour_matrix(2, cur_idx);
        element_nums = [element_nums num_elements];
        if num_elements > least_polygon_side_num
            contour_coord_arr{level_idx} = contour_matrix(:, cur_idx+1:cur_idx+num_elements)';
            level_idx = level_idx + 1;
        end
        cur_idx = cur_idx + num_elements+1;
        if cur_idx > length(contour_matrix)
            done = 1;
        end
    end
end

function [crossed, m, n] = line_cross(l1_p1, l1_p2, l2_p1, l2_p2)
    ax = l1_p1(1); ay = l1_p1(2); bx = l1_p2(1); by = l1_p2(2);
    cx = l2_p1(1); cy = l2_p1(2); dx = l2_p2(1); dy = l2_p2(2);
    crossed = false; m = Inf; n = Inf;
    if ((cx-dx)*(by-ay)-(cy-dy)*(bx-ax))
        m=((cx-ax)*(by-ay)-(cy-ay)*(bx-ax))/((cx-dx)*(by-ay)-(cy-dy)*(bx-ax));
        if ((bx-ax)*(cy-dy)-(by-ay)*(cx-dx))
            n=((cx-ax)*(cy-dy)-(cy-ay)*(cx-dx))/((bx-ax)*(cy-dy)-(by-ay)*(cx-dx));
            if ((0 <= m) && (m <= 1) && (0 <= n) && (n <= 1))
                crossed = true;
            end
        end
    end
end

function [crossed, intersect_coord] = is_crossed(coord, contour2_coor, thres_distance)
    %
    % coord: [prev_point, current_point, next_point, next_2_point] of the first contour
    % contour2_coor: the coordinates of the second contour
    % thres_distance: The threshold of distance in pixel that the point will be considered
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    prev_coor = coord(1,:); coor = coord(2,:); next_coor = coord(3,:); next2_coor = coord(4,:);
    distances = sqrt(sum((coor-contour2_coor).^2, 2));
    close_idx = find(distances<=thres_distance);
    crossed = false; intersect_coord = 0;
    for i = 1:length(close_idx)
        idx = close_idx(i);
        if idx < length(contour2_coor)
            c2_prev_coor = [0 0]; c2_next2_coor = [0 0];
            if idx > 1
                c2_prev_coor = contour2_coor(idx-1,:);
            end
            c2_coor = contour2_coor(idx,:);
            c2_next_coor = contour2_coor(idx+1,:);
            if idx < length(contour2_coor)-2
                c2_next2_coor = contour2_coor(idx+2,:);
            end
            [crossed, m, ~] = line_cross(coor, next_coor, c2_coor, c2_next_coor);
            if crossed
                backcrosses = [0 0 0 0 0];
                if c2_next2_coor
                    [backcrosses(1), ~, ~] = line_cross(coor, next_coor, c2_next_coor, c2_next2_coor);
                end
                if next2_coor
                    [backcrosses(2), ~, ~] = line_cross(next_coor, next2_coor, c2_coor, c2_next_coor);
                    if c2_next2_coor
                        [backcrosses(3), ~, ~] = line_cross(next_coor, next2_coor, c2_next_coor, c2_next2_coor);
                    end
                end
                if prev_coor
                    [backcrosses(4), ~, ~] = line_cross(prev_coor, coor, c2_coor, c2_next_coor);
                    if c2_prev_coor
                        [backcrosses(5), ~, ~] = line_cross(prev_coor, coor, c2_prev_coor, c2_coor);
                    end
                end
                if any(backcrosses)
                    crossed = false;
                else
                    crossed = true;
                    x = c2_coor(1) + m * (c2_next_coor(1) - c2_coor(1));
                    y = c2_coor(2) + m * (c2_next_coor(2) - c2_coor(2));
                    intersect_coord = [x y];
                    break
                end
            end
        end
    end
end

function cross = get_intersection(contour1_coor, contour2_coor, thres_distance)
    %
    % contour1_coor: Contour coordinates of OD columns
    % contour2_coor: Contour coordinates of another map (spatial frequency, orientations, etc.)
    % thres_distance: The threshold of distance in pixel that the point will be considered
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cross = [];
    for i = 1:length(contour1_coor)-1
        coor = contour1_coor(i, :);
        prev_coor = [0 0]; next2_coor = [0 0];
        if i > 1
            prev_coor = contour1_coor(i-1, :);
        end
        next_coor = contour1_coor(i+1, :);
        if i < length(contour1_coor)-2
            next2_coor = contour1_coor(i+2, :);
        end
        [crossed, intersection] = is_crossed([prev_coor; coor; next_coor; next2_coor], contour2_coor, thres_distance);
        if crossed
            cross = [cross; intersection];
        end
    end
end

function angle_degs = get_intersect_angle(map1, map2, intersections)
    %
    % intersections: the xy-indices of intersection in map1 and map2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    angle_degs = [];
    [map1_field_x, map1_field_y] = gradient(map1);
    [map2_field_x, map2_field_y] = gradient(map2);
    for i = 1:length(intersections)
        indx = intersections(i, 1); indy = intersections(i, 2);
        vx = map1_field_x(indy, indx);
        vy = map1_field_y(indy, indx);
        ux = map2_field_x(indy, indx);
        uy = map2_field_y(indy, indx);
        m = [vx vy];
        n = [ux uy];
        cos_angle = sum(m.*n) / (sqrt(sum(m.^2)) * sqrt(sum(n.^2)));
        angle_deg = rad2deg(acos(cos_angle));
        if angle_deg > 90
            angle_deg = 180-angle_deg;
        end
        angle_degs = [angle_degs angle_deg];
    end
end

  
function freqs_out = make_freq_polar_format(freqs_input)
    % polar plot
    freq2 = [freqs_input fliplr(freqs_input)];
    freq3 = [freq2 freq2];
    %freq3 = smooth(freq3);
    %freqs_out = freq3 / max(freq3(:)) ;
    %freqs_out(end+1) = freqs_out(1);
    freq4= [freq3 freq3 freq3]; 
    freq5 = smooth(freq4); 
    freq6 = freq5( length(freq3)+1 : 2*length(freq3) ); 
    freqs_out = freq6 / max(freq6(:)); 
    freqs_out(end+1) = freqs_out(1);
end 