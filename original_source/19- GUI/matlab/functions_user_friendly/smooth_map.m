function [ONOFF_smoothed,ODCrtxPlt_smooth,ONOFF_interpolated,ODCrtxPlt_interpolated] = ...
                                        smooth_map(ODCrtxPlt,ONOFFCrtxPlt,n_interp,show_fig)
    % OD/ONOFF smoothing 
    
    kern_std = 2; %2  % The std of kernel for smoothing the OD map
    %n_interp = 1;   % The number of interpolation for OD map (OD and ori maps must have same size, use 0 for non-smoothed ori map and use 2 for smoothed ori map)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    kern = fspecial('gaussian', kern_std*6, kern_std);
    ODCrtxPlt_smooth = conv2(ODCrtxPlt, kern, 'same');
    ODCrtxPlt_smooth = (ODCrtxPlt_smooth - min(min(ODCrtxPlt_smooth))) ./ (max(max(ODCrtxPlt_smooth)) - min(min(ODCrtxPlt_smooth)));
    %ODCrtxPlt_plot = interp2(ODCrtxPlt_smooth, n_interp);
    ODCrtxPlt_interpolated = imresize(ODCrtxPlt_smooth, n_interp);

    ONOFF_smoothed = conv2(ONOFFCrtxPlt, kern, 'same');
    ONOFF_smoothed = (ONOFF_smoothed - min(min(ONOFF_smoothed))) ./ (max(max(ONOFF_smoothed)) - min(min(ONOFF_smoothed)));
%     ONOFF_plot = interp2(ONOFF_smoothed, n_interp);
    ONOFF_interpolated = imresize(ONOFF_smoothed, n_interp);
    
    % - Plot ON OFF OD on top 
    if show_fig == 1
        figure,
        od_contour_levels = 1;
        map1 = [133 94 82;255 153 153;33 31 133;153 153 255]./255;
        set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
        z1 = (double((ODCrtxPlt>0))+1)*1 ;
        z2 = (double(ONOFFCrtxPlt==-1)+0)*2 ;        
        subplot(131), imagesc(z1+z2)
        colormap(gca, map1); axis square; title('Raw Map','fontsize',20)
        %   hold on, contour(ODCrtxPlt_smooth, od_contour_levels, 'k', 'LineWidth', 5);
        % hold on , plot(2:2:60,ones(1,30)*25,'k.','MarkerSize',10)
        % hold on , plot(37,30,'wo','MarkerSize',10,'MarkerFaceColor','w'),viscircles([37 30],10,'color','w')
        
        z3 = (double((ODCrtxPlt_smooth>0.5))+1)*1 ;  % figure, imagesc(z1)
        z4 = (double(ONOFF_smoothed<.5)+0)*2 ;
        subplot(132), imagesc(z3+z4)
        colormap(gca, map1); axis square; hold on, 
        title('Smooth Map','fontsize',20)
        hold on, contour(ODCrtxPlt_smooth, od_contour_levels, 'k', 'LineWidth', 5);
        %   xlim([.5 8.5]), ylim([.5 8.5])
        
        z5 = (double((ODCrtxPlt_interpolated>0.5))+1)*1 ;  % figure, imagesc(z1)
        z6 = (double(ONOFF_interpolated<.5)+0)*2 ;
        subplot(133), imagesc(z5+z6)
        colormap(gca, map1); axis square; hold on, contour(ODCrtxPlt_interpolated, od_contour_levels, 'k', 'LineWidth', 5);
        title('Smooth Map (interpolated)','fontsize',20)
        %   xlim([1 33]), ylim([1 33])
    end 

end 