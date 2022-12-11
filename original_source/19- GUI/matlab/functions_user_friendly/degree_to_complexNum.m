function [complex_map] = degree_to_complexNum(ori_map)
    map_tmp = (ori_map - 90) * (pi/180) * 2;
    real = cos(map_tmp);
    imag = sin(map_tmp);
    complex_map = real + imag * 1i;
end