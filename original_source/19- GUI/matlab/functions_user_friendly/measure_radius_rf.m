function r_fit_circle = measure_radius_rf(rowRFcenter,colRFcenter,allCXrfTempON,allCXrfTempOFF,alpha)
%  measure the radius from the center of RF to a percentage (alpha) of the
%  maximum response for iso-retinotopic patch measurement
   %rf_abs = abs(allCXrfTempON) + abs(allCXrfTempOFF);
   rf_abs = abs(allCXrfTempON + allCXrfTempOFF);
   rf_thresh = max(rf_abs(:)) * alpha ;
   rf_diff = abs(rf_abs - rf_thresh);% figure, imagesc(rf_diff), colorbar
   valid_data1 = rf_diff < .1; % the difference with the set alpha is less than 10%
    
   not_valid_area = rf_abs < .05 ; % removing the least 5% value of all the map 

   valid_data = valid_data1 .* (~not_valid_area); 
   [point_row,point_col] = find(valid_data); 
   dist_center_rf = sqrt((point_row(:)-rowRFcenter).^2 + (point_col(:)-colRFcenter).^2); 
   
   r_not_valid = measure_r_not_valid(allCXrfTempON,allCXrfTempOFF,rowRFcenter,colRFcenter); 
   dist_center_rf( dist_center_rf < r_not_valid) = []; 
   r_fit_circle = round(mean(dist_center_rf));
   
   if isnan(r_fit_circle)
       r_fit_circle = r_not_valid; 
   end 
       
   if 1 == 0
       figure,
       rf_norm = allCXrfTempON + allCXrfTempOFF; 
       rf_norm = rf_norm / max(abs(rf_norm(:))); 
       subplot(121),imagesc(rf_norm), colormap(gca,'jet'), axis square
       subplot(122),imagesc(rf_abs), hold on, plot(colRFcenter,rowRFcenter,'ko'), axis square
       hold on
       viscircles([colRFcenter,rowRFcenter],r_fit_circle)
   end
end

function  r_remove = measure_r_not_valid(allCXrfTempON,allCXrfTempOFF,rowRFcenter,colRFcenter)
% removing the range for RF from the center of RF and center of each polarity 
% 
% r_remove : distance of center of RF and max center of ON or OFF RF separately
    
   if  max(abs(allCXrfTempON(:))) > 0    
       [r_on,c_on] = find(max(abs(allCXrfTempON(:))) == allCXrfTempON); 
       dist_on = norm([rowRFcenter-r_on,colRFcenter-c_on]); 
   else 
       dist_on = 0; 
   end 
   
   if  max(abs(allCXrfTempOFF(:))) > 0 
       [r_off,c_off] = find(max(abs(allCXrfTempOFF(:)))== allCXrfTempOFF); 
       dist_off = norm([rowRFcenter-r_off,colRFcenter-c_off]);
   else 
       dist_off = 0; 
   end 
   
   if dist_on>dist_off
       lim_dist = dist_on; 
   else 
       lim_dist = dist_off; 
   end 
   r_remove = round(lim_dist);
   
%    temp1 = zeros(size(rf_abs));
%    temp1(rowRFcenter,colRFcenter) = 1;
%    disk_filt = strel('disk',r_remove,0);
%    not_valid_area_in = imdilate(temp1,disk_filt); % figure, imagesc(not_valid_area)
%    
%    not_valid_area_out = rf_abs < .05 ; % removing the least 5% value of all the map 
%    not_valid_area = not_valid_area_in + not_valid_area_out ; 
end