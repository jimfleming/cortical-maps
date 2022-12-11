function [rf_center_row, rf_center_col, sigma_rfs] = find_center_sigma_rfs_S4( rf_input, row_range, col_range, minsigma, maxsigma, sigma_smooth, sd_factor, debug )
% find center of RF by measuring the average position of x and y of abs(rfon+rfoff)
% sigma of gaussian fit to the abs(rfon+rfoff)
% rf_input   :   rf_mature_contra or rf_mature_ipsi 
% Searching limit sigmas to find the best Gaussian fit to the abs(rfon + rfoff)
%  minsigma = 10 ;   
%  maxsigma = 30 ;
% sigma_smooth = 5,  to make abs(rf_onoff) uniform for Gaussian fit 
% sd_factor      :   standard deviation factor : R estimate of Cortical Receptive field = estimated_sigma * sd_factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rf_on_all = rf_input.allCXrfON ;
    rf_off_all = rf_input.allCXrfOFF ;
    rf_size = size( rf_on_all{1} , 1) ;
    [ XX , YY ] = meshgrid( 1:rf_size , 1:rf_size ); 

    sigma_rfs = zeros( length(row_range), length(col_range)); 
    rf_center_row = zeros( length(row_range), length(col_range)); 
    rf_center_col = zeros( length(row_range), length(col_range)); 
    
    for ii = row_range  
        for jj = col_range
            rf_on = rf_on_all{ ii , jj };
            rf_off = rf_off_all{ ii , jj };
            rf_onoff = rf_on + rf_off;
            max_rf = max(max(abs(rf_on(:))),max(abs(rf_off(:))));
            rf_norm = rf_onoff / max_rf; 
            rf_abs = abs(rf_norm); 

            % measure rf center by averaging the rf_abs
            sum_rf = sum( rf_abs , [1,2]);
            rf_xx = sum(rf_abs .* XX , [1,2]) / sum_rf ; 
            rf_yy = sum(rf_abs .* YY , [1,2]) / sum_rf ;
            rf_center_row(ii, jj) = rf_yy; 
            rf_center_col(ii, jj) = rf_xx; 

%             rf_abs = imgaussfilt(rf_abs,sigma_smooth); % to remove onoff boundary 
%             [optsigma] = autoGaussianSurf( XX , YY , rf_abs , minsigma , maxsigma , rf_xx , rf_yy ) ; 
%             R_rf = optsigma * sd_factor ; 
%             sigma_rfs(ii, jj) = R_rf; 
            
            % debug 
            if debug
                figure(80), clf
                imagesc( rf_norm ), axis('square'), %colorbar 
                hold on, plot( rf_xx , rf_yy , 'ko')
                hold on, viscircles( [ rf_xx , rf_yy ] , R_rf ,'color', 'k', 'linewidth', 5), axis off 
                colormap jet 
                %title( num2str(R_rf) , 'fontsize', 16)
                pause
            end 

        end
    end

end 

function [optsigma] = autoGaussianSurf(xi , yi , zi , minsigma , maxsigma , x0 , y0 )

    xi = xi(:);
    yi = yi(:);
    
    sigmas = minsigma : maxsigma ; 

    res = zeros(length(sigmas),7);
    
    % Run through all the different values for sigma
    for ii = 1:length(sigmas)
        % Determine the residual error for the optimal x, y for this sigma
        G = exp(-((xi-x0).^2+(yi-y0).^2)/2/sigmas(ii)^2);
        X = [G,ones(length(G),1)];
        ps = X\zi(:);
        res(ii,:) = [sum((zi(:) - X*ps).^2),ps(:)',x0,y0,sigmas(ii),sigmas(ii)];
    end
    
    % Find sigma with corresponding least error
    [~,optsigma_ind] = min(res(:,1));
    optsigma = sigmas(optsigma_ind); 
end 
