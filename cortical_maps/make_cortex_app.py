import numpy as np


def make_cortex_app(appdata):
    CrtxLength = appdata.CrtxLength
    CrtxWidth = appdata.CrtxWidth
    onoff_rf_space = appdata.onoff_rf_space
    od_rf_space = appdata.od_rf_space
    RetinotopyRFspace_plot = appdata.RetinotopyRFspace_plot
    IndRetinaRfspace = appdata.IndRetinaRfspace
    rng_trial = appdata.rng_trial

    CrtxONOFF3m = np.zeros(shape=(CrtxLength, CrtxWidth))
    CrtxOD3m = np.zeros(shape=(CrtxLength, CrtxWidth))
    Retinotopy3mIndPlot = np.zeros(
        shape=(CrtxLength, CrtxWidth)
    )  # index based on color map
    Retinotopy3mIndex = np.zeros(
        shape=(CrtxLength, CrtxWidth)
    )  # index based on color map
    connection_ind = map_retinotopically(
        onoff_rf_space != 0, CrtxOD3m.shape, False, 1, rng_trial
    )

    CrtxONOFF3m.flat[:] = onoff_rf_space.flat[connection_ind]
    CrtxOD3m.flat[:] = od_rf_space.flat[connection_ind]
    Retinotopy3mIndPlot.flat[:] = RetinotopyRFspace_plot.flat[connection_ind]
    Retinotopy3mIndex.flat[:] = IndRetinaRfspace.flat[connection_ind]

    appdata.CrtxONOFF3m = CrtxONOFF3m
    appdata.CrtxOD3m = CrtxOD3m
    appdata.Retinotopy3mIndPlot = Retinotopy3mIndPlot
    appdata.Retinotopy3mIndex = Retinotopy3mIndex


def map_retinotopically(
    input_cell_bin_arr, target_shape, preconnect, termination_thresh, rng_trial
):
    # input_cell_bin_arr: a 2D binary array representing the position of input cells
    # target_shape: the shape of the target, [num_row, num_col]; num_row * num_col = length(find(input_cell_bin_arr))
    # preconnect: {0, 1}, if 1, the input cells will be connected to the closest postnaptic targets first before optimization
    # termination_thresh: the number of iterations without improvement for termination of optimization process, max = num_targets
    #
    # Return
    # connection_ind: the indices for connections, i.e. target(:) = input_cell_bin_arr(connection_ind)
    rng = np.random.default_rng(seed=rng_trial)
    connection_ind = np.zeros(shape=(target_shape[0] * target_shape[1]), dtype=np.int64)

    # normalize the x- and y-distances
    rr, cc = np.nonzero(input_cell_bin_arr)
    rr_norm = (rr - rr.min()) / (rr.max() - rr.min())
    cc_norm = (cc - cc.min()) / (cc.max() - cc.min())

    if preconnect:
        target_rr_norm = np.arange(target_shape[0]) / (target_shape[0] - 1)
        target_cc_norm = np.arange(target_shape[1]) / (target_shape[1] - 1)
        rr_norm_tmp = rr_norm.copy()
        cc_norm_tmp = cc_norm.copy()
        input2target_conn_ind = np.zeros_like(connection_ind)
        target_positions = np.zeros(shape=(len(connection_ind), 2), dtype=np.int64)

        for c in range(target_shape[1]):  # column index
            for r in range(target_shape[0]):  # row index
                distances = np.linalg.norm(
                    np.stack([target_rr_norm[r], target_cc_norm[c]], axis=-1)
                    - np.stack([rr_norm_tmp, cc_norm_tmp], axis=-1),
                    axis=-1,
                )
                sorted_ind = np.argsort(distances)

                (jj,) = np.nonzero(rr_norm == rr_norm_tmp[sorted_ind[0]])
                (ii,) = np.nonzero(cc_norm == cc_norm_tmp[sorted_ind[0]])
                _, (i,) = np.nonzero(ii[None] == jj[None].T)

                connection_ind[
                    np.ravel_multi_index((r, c), target_shape)
                ] = np.ravel_multi_index(
                    (rr[ii[i]], cc[ii[i]]), input_cell_bin_arr.shape
                )
                target_positions[np.ravel_multi_index((r, c), target_shape)] = [r, c]
                input2target_conn_ind[np.ravel_multi_index((r, c), target_shape)] = ii[
                    i
                ]

                rr_norm_tmp = np.delete(rr_norm_tmp, sorted_ind[0])
                cc_norm_tmp = np.delete(cc_norm_tmp, sorted_ind[0])
    else:
        X, Y = np.meshgrid(np.arange(target_shape[0]), np.arange(target_shape[1]))
        C = np.stack([X, Y], axis=-1)
        target_positions = C.reshape(-1, 2)
        input2target_conn_ind = rng.permutation(len(rr))

    target_rr = target_positions[:, 0]
    target_cc = target_positions[:, 1]
    target_rr_norm = (target_rr - min(target_rr)) / (max(target_rr) - min(target_rr))
    target_cc_norm = (target_cc - min(target_cc)) / (max(target_cc) - min(target_cc))

    input2target_conn_ind = get_minimal_wiring(
        np.stack([rr_norm, cc_norm], axis=-1),
        np.stack([target_rr_norm, target_cc_norm], axis=-1),
        input2target_conn_ind,
        termination_thresh,
        rng_trial,
    )
    connection_ind = np.ravel_multi_index(
        (rr[input2target_conn_ind], cc[input2target_conn_ind]), input_cell_bin_arr.shape
    )
    return connection_ind


def get_minimal_wiring(
    input_positions, target_positions, connection_ind, termination_thresh, rng_trial
):
    # input_positions: the xy-coordinates of the inputs, shape: (num_inputs, 2)
    # target_positions: the xy-coordinates of the targets, shape: (num_targets, 2), where num_targets = num_inputs.
    # connection_ind: the indices for connections, i.e. target_positions receives inputs from input_positions(connection_ind)
    # termination_thresh: the number of iterations without improvement for termination, max = num_targets
    rng = np.random.default_rng(seed=rng_trial)
    count = 0
    while count < termination_thresh:
        distances = np.linalg.norm(
            target_positions - input_positions[connection_ind], axis=-1
        )
        sorted_ind = np.argsort(distances)
        cur_ind = len(sorted_ind) - count - 1
        cur_longest_conn = distances[sorted_ind[cur_ind]]
        for ii in range(len(connection_ind)):
            test_ind = rng.integers(len(sorted_ind))
            cur_test_conn = distances[sorted_ind[test_ind]]
            cur_dist_sum = cur_longest_conn + cur_test_conn
            new_input_pos = np.stack(
                [
                    input_positions[connection_ind[sorted_ind[test_ind]]],
                    input_positions[connection_ind[sorted_ind[cur_ind]]],
                ],
                axis=0,
            )
            new_dist = np.linalg.norm(
                np.stack(
                    [
                        target_positions[sorted_ind[cur_ind]],
                        target_positions[sorted_ind[test_ind]],
                    ],
                    axis=0,
                )
                - new_input_pos,
                axis=-1,
            )
            if sum(new_dist) < cur_dist_sum:
                if max(new_dist) < cur_longest_conn:
                    connection_ind[
                        [sorted_ind[cur_ind], sorted_ind[test_ind]]
                    ] = connection_ind[[sorted_ind[test_ind], sorted_ind[cur_ind]]]
                    count = 0
                    break
            if ii == len(connection_ind) - 1:
                count = count + 1

    return connection_ind
