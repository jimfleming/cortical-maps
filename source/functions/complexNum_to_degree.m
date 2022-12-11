function [ori_map] = complexNum_to_degree(complex_map)
    ori_map = angle(complex_map) / 2 * (180/pi) + 90;
end