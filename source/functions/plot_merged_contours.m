function [] = plot_merged_contours(map, contour_levels, buffer_factor, color, linewidth)
    %
    % to plot merged contours for the map with discontinuos surface
    % buffer_factor: [0, 1], higher buffer_factor squeezes the contours less and thus having smoother lines
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    map_mx = max(max(map));
    new_mx = map_mx * 1000;
    contour_matrix = contourc(map, contour_levels);
    clevels = get_contour_levels(contour_matrix);
    clevels = sort(clevels);
    map_tmp = map;
    buffer_value = (clevels(end) - clevels(end-1)) * buffer_factor;
    map_tmp(map_tmp > clevels(end)+buffer_value) = new_mx;
    for cc = 1:length(clevels)
        contour(map_tmp, [clevels(cc) clevels(cc)], color, 'LineWidth', linewidth);
    end
end

function [clevels] = get_contour_levels(contour_matrix)
    done = 0;
    cur_idx = 1;
    clevels = [];
    while done == 0
        clevel = contour_matrix(1, cur_idx);
        num_elements = contour_matrix(2, cur_idx);
        if find(clevels == clevel)
            % do nothing
        else
            clevels = [clevels clevel];
        end
        cur_idx = cur_idx + num_elements + 1;
        if cur_idx > length(contour_matrix)
            done = 1;
        end
    end
end