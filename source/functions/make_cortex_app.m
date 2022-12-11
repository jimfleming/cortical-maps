function [appdata] = make_cortex_app(appdata)


CrtxLength = appdata.CrtxLength;
CrtxWidth = appdata.CrtxWidth;
onoff_rf_space = appdata.onoff_rf_space;
od_rf_space= appdata.od_rf_space;
RetinotopyRFspace_plot = appdata.RetinotopyRFspace_plot;
IndRetinaRfspace = appdata.IndRetinaRfspace;
cmap = appdata.cmap;
rng_trial = appdata.rng_trial;
show_fig = 0;

CrtxONOFF3m = zeros(CrtxLength, CrtxWidth);
CrtxOD3m = zeros(CrtxLength, CrtxWidth);
Retinotopy3mIndPlot = zeros(CrtxLength, CrtxWidth); % index based on color map
Retinotopy3mIndex = zeros(CrtxLength, CrtxWidth);   % index based on color map
connection_ind = map_retinotopically(onoff_rf_space ~= 0, size(CrtxOD3m), 1, 1, rng_trial);
CrtxONOFF3m(:) = onoff_rf_space(connection_ind);
CrtxOD3m(:) = od_rf_space(connection_ind);
Retinotopy3mIndPlot(:) = RetinotopyRFspace_plot(connection_ind);
Retinotopy3mIndex(:) = IndRetinaRfspace(connection_ind);


appdata.CrtxONOFF3m = CrtxONOFF3m;
appdata.CrtxOD3m = CrtxOD3m;
appdata.Retinotopy3mIndPlot = Retinotopy3mIndPlot;
appdata.Retinotopy3mIndex = Retinotopy3mIndex;

end


function [connection_ind] = map_retinotopically(input_cell_bin_arr, target_shape, preconnect, termination_thresh, rng_trial)
%
% input_cell_bin_arr: a 2D binary array representing the position of input cells
% target_shape: the shape of the target, [num_row, num_col]; num_row * num_col = length(find(input_cell_bin_arr))
% preconnect: {0, 1}, if 1, the input cells will be connected to the closest postnaptic targets first before optimization
% termination_thresh: the number of iterations without improvement for termination of optimization process, max = num_targets
%
% Return
% connection_ind: the indices for connections, i.e. target(:) = input_cell_bin_arr(connection_ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(rng_trial)

connection_ind = zeros(1, target_shape(1) * target_shape(2));
% normalize the x- and y-distances
[rr, cc] = find(input_cell_bin_arr);
rr_norm = (rr - min(rr)) / (max(rr) - min(rr));
cc_norm = (cc - min(cc)) / (max(cc) - min(cc));

if logical(preconnect)
    target_rr_norm = (0:target_shape(1)-1) / (target_shape(1)-1);   % min = 1
    target_cc_norm = (0:target_shape(2)-1) / (target_shape(2)-1);
    rr_norm_tmp = rr_norm;
    cc_norm_tmp = cc_norm;
    input2target_conn_ind = zeros(1, length(connection_ind));
    for c = 1:target_shape(2)   % column index
        for r = 1:target_shape(1)   % row index
            distances = vecnorm([target_rr_norm(r) target_cc_norm(c)] - [rr_norm_tmp cc_norm_tmp], 2, 2);
            [~, sorted_ind] = sort(distances);
            jj = find(rr_norm == rr_norm_tmp(sorted_ind(1)));
            ii = find(cc_norm == cc_norm_tmp(sorted_ind(1)));
            [i, ~] = find(ii==jj');
            connection_ind(sub2ind(target_shape, r, c)) = sub2ind(size(input_cell_bin_arr), rr(ii(i)), cc(ii(i)));
            target_positions(sub2ind(target_shape, r, c), :) = [r c];
            input2target_conn_ind(sub2ind(target_shape, r, c)) = ii(i);
            rr_norm_tmp(sorted_ind(1)) = [];
            cc_norm_tmp(sorted_ind(1)) = [];
        end
    end
else   % use random permutation for pre-connection
    [X, Y] = meshgrid(1:target_shape(1), 1:target_shape(2));
    C = cat(2, X', Y');
    target_positions = reshape(C, [], 2);
    input2target_conn_ind = randperm(length(rr));
end

target_rr = target_positions(:, 1);
target_cc = target_positions(:, 2);
target_rr_norm = (target_rr - min(target_rr)) / (max(target_rr) - min(target_rr));
target_cc_norm = (target_cc - min(target_cc)) / (max(target_cc) - min(target_cc));

input2target_conn_ind = get_minimal_wiring([rr_norm cc_norm], [target_rr_norm target_cc_norm], input2target_conn_ind, termination_thresh, rng_trial);
connection_ind = sub2ind(size(input_cell_bin_arr), rr(input2target_conn_ind), cc(input2target_conn_ind));
end


function [connection_ind] = get_minimal_wiring(input_positions, target_positions, connection_ind, termination_thresh, rng_trial)
%
% input_positions: the xy-coordinates of the inputs, shape: (num_inputs, 2)
% target_positions: the xy-coordinates of the targets, shape: (num_targets, 2), where num_targets = num_inputs.
% connection_ind: the indices for connections, i.e. target_positions receives inputs from input_positions(connection_ind)
% termination_thresh: the number of iterations without improvement for termination, max = num_targets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(rng_trial)

count = 0;
while count < termination_thresh
    distances = vecnorm(target_positions - input_positions(connection_ind, :), 2, 2);
    [~, sorted_ind] = sort(distances);
    cur_ind = length(sorted_ind) - count;
    cur_longest_conn = distances(sorted_ind(cur_ind));
    for ii = 1:length(connection_ind) - 1
        test_ind = randi(length(sorted_ind));
        cur_test_conn = distances(sorted_ind(test_ind));
        cur_dist_sum = cur_longest_conn + cur_test_conn;
        new_input_pos = [input_positions(connection_ind(sorted_ind(test_ind)), :); input_positions(connection_ind(sorted_ind(cur_ind)), :)];
        new_dist = vecnorm([target_positions(sorted_ind(cur_ind), :); target_positions(sorted_ind(test_ind), :)] - new_input_pos, 2, 2);
        if sum(new_dist) < cur_dist_sum
            if max(new_dist) < cur_longest_conn
                connection_ind([sorted_ind(cur_ind) sorted_ind(test_ind)]) = connection_ind([sorted_ind(test_ind) sorted_ind(cur_ind)]);
                count = 0;
                break
            end
        end
        if ii == length(connection_ind) - 1
            count = count + 1;
        end
    end
end
end