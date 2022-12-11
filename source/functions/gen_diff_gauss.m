function [kern] = gen_diff_gauss(xlen, ylen, A, B, beta, d_1, d_2, shift)
    if nargin < 8
        shift = 0;
    end
    x = -floor(xlen/2):floor(xlen/2);
    y = -floor(ylen/2):floor(ylen/2);
    [X, Y] = meshgrid(x, y);
    if shift == 0
        kern = A * exp(-(beta * X.^2 + Y.^2) / d_1) - B * exp(-(X.^2 + Y.^2) / d_2);
    else
        kern = A * exp(-(beta * X.^2 + (Y-shift).^2) / d_1) - B * exp(-(beta * X.^2 + (Y+shift).^2) / d_2);
    end
end