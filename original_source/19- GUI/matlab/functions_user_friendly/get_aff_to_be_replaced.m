function [Aff_to_be_replaced, Conv_arr, Conv_Result, srch_r, srch_c] = get_aff_to_be_replaced(pos_arr, pos_arr_mask, OD_sorted, ONOFF_sorted)
    [i1, j1] = find(pos_arr);
    mask_search = filter2(ones(3), pos_arr) - pos_arr;
    if nargin < 4
        arr_interested = OD_sorted;
        [srch_r, srch_c] = find(mask_search & arr_interested ~= arr_interested(i1, j1));
    else
        arr_interested = ONOFF_sorted;
        [srch_r, srch_c] = find(mask_search & arr_interested ~= arr_interested(i1, j1) & OD_sorted == OD_sorted(i1, j1));
    end
    
    Conv_CC = (arr_interested > 0) .* pos_arr_mask .* (pos_arr_mask > 0);   % Convolution Center Contra
    Conv_CI = (arr_interested < 0) .* pos_arr_mask .* (pos_arr_mask > 0);   % Convolution Center Ipsi
    Conv_SC = (arr_interested > 0) .* pos_arr_mask .* (pos_arr_mask < 0);   % Convolution Surround Contra
    Conv_SI = (arr_interested < 0) .* pos_arr_mask .* (pos_arr_mask < 0);   % Convolution Surround Ipsi
    Sum_CC = sum(abs(Conv_CC(:)));
    Sum_CI = sum(abs(Conv_CI(:)));
    Sum_SC = sum(abs(Conv_SC(:)));
    Sum_SI = sum(abs(Conv_SI(:)));
    if arr_interested(i1, j1) == 1
        Conv_Result = (Sum_CC - Sum_CI) - (Sum_SC - Sum_SI);
    elseif arr_interested(i1, j1) == -1
        Conv_Result = (Sum_CI - Sum_CC) - (Sum_SI - Sum_SC);
    end
    Conv_arr = zeros(size(srch_r, 1), 1);
    Aff_to_be_replaced = zeros(size(srch_r, 1), 1);
    
    for k1 = 1 : size(srch_r, 1)
        temp_result = arr_interested;
        temp_result(i1, j1) = arr_interested(srch_r(k1), srch_c(k1));
        temp_result(srch_r(k1), srch_c(k1)) = 1;
        temp_mask_Contra = circshift(pos_arr_mask, [srch_r(k1)-i1 srch_c(k1)-j1]);

        Conv_CC_temp = sum(sum((temp_result > 0) .* temp_mask_Contra .* (temp_mask_Contra > 0)));   % Convolution Center Contra
        Conv_CI_temp = sum(sum((temp_result < 0) .* temp_mask_Contra .* (temp_mask_Contra > 0)));   % Convolution Center Ipsi
        Conv_SC_temp = sum(sum((temp_result > 0) .* temp_mask_Contra .* (temp_mask_Contra < 0)));   % Convolution Surround Contra
        Conv_SI_temp = sum(sum((temp_result < 0) .* temp_mask_Contra .* (temp_mask_Contra < 0)));   % Convolution Surround Ipsi
        if arr_interested(i1, j1) == 1
            Conv_arr(k1) = (Conv_CC_temp - Conv_CI_temp) - (Conv_SC_temp - Conv_SI_temp);
        elseif arr_interested(i1, j1) == -1
            Conv_arr(k1) = (Conv_CI_temp - Conv_CC_temp) - (Conv_SI_temp - Conv_SC_temp);
        end
        Aff_to_be_replaced(k1) = arr_interested(srch_r(k1), srch_c(k1));
    end
end