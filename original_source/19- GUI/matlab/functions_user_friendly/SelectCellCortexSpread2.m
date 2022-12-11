function [SelectionCondition, mask_spread] = SelectCellCortexSpread2(dilate_str, ODMap, rAffCrtx, cAffCrtx)
    % selecting neighboring cells in cortex based on spread value
    % the spread function is a Gaussian distribution
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dilate_str = strel('disk',round(sdspread),0);
    mask_spread = zeros(size(ODMap));
    mask_spread(rAffCrtx, cAffCrtx) = 1;
    mask_spread = imdilate(mask_spread, dilate_str); %   figure,imagesc(mask_spread)

    % Selecting those cells in retina that have similar organization to cortex
    SelectionCondition = ODMap == 1 & mask_spread == 1;
end