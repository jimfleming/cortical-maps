function mask_spread = imdialte_grid(row_grid, col_grid, ii, jj, sdspread)
% replacing imdilate function of Matlab (to make the code faster)
% the function uses the meshgrid with the number of rows and cols of input image 

%%%%%%%%%%%%%%%%%%%%%%%
    dist_matrix = sqrt((col_grid - jj).^2 + (row_grid - ii).^2);
    mask_spread = dist_matrix <= sdspread; 

end