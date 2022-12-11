function [pw_locations] = get_pinwheel_location_and_charge_coverage(ori_map,show_fig, border_size )
    %
    % ori_map: the 2D orientation map
    % border_size: the size (in pixel) of the borders that should have no pinwheels
    %
    % Return:
    % pw_locations: the xy-coordinates of the pinwheel centers; size = [num_pinwheel, 2]
    % pw_charges: the charges of the pinwheels. 1: clockwise; -1: anticlockwise.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nargin < 3
        border_size = 0;
    end
    
    z = degree_to_complexNum(ori_map);
    
    if show_fig == 1
        figure;
        imagesc(ori_map); axis square; colormap(hsv); colorbar; hold on;
        contour(real(z), [0 0], 'r', 'LineWidth', 2);
        contour(imag(z), [0 0], 'b', 'LineWidth', 2);
    end
    
    cReal = contourc(real(z), [0 0]);
    cImag = contourc(imag(z), [0 0]);

    if isempty(cImag) || isempty(cReal)
        pw_locations = [];
    else 
        % get the contour coordinates from contours with at least four points
        [cReal_coord, ~] = get_contour_coord(cReal, 3);
        [cImag_coord, ~] = get_contour_coord(cImag, 3);
        
        pw_locations = [];
        for i = 1:length(cImag_coord)
            polygonImag = cImag_coord{i};
            polyIm_distance = sqrt(sum((polygonImag(1:end-1, :)-polygonImag(2:end, :)).^2, 2));
            for j = 1:length(cReal_coord)
                polygonReal = cReal_coord{j};
                polyRe_distance = sqrt(sum((polygonReal(1:end-1, :)-polygonReal(2:end, :)).^2, 2));
                thres_distance = 2 * max(max(polyRe_distance), max(polyIm_distance));
                cross_tmp = get_intersection(polygonReal, polygonImag, thres_distance);
                if size(cross_tmp)
                    pw_locations = [pw_locations; cross_tmp];
                end
            end
        end
        pw_locations_int = round(pw_locations);
        [~, idx, ~] = unique(pw_locations_int, 'rows');
        pw_locations_int = pw_locations_int(sort(idx), :);
        pw_locations = pw_locations(sort(idx), :);
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
        cur_idx = cur_idx + num_elements + 1;
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
    % coord: [prev_point, current_point, next_point, next_2_point] of the first contour
    % contour2_coor: the coordinates of the second contour
    % thres_distance: The threshold of distance in pixel that the point will be considered
    
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
    % contour1_coor: Contour coordinates of OD columns
    % contour2_coor: Contour coordinates of another map (spatial frequency, orientations, etc.)
    % thres_distance: The threshold of distance in pixel that the point will be considered
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
