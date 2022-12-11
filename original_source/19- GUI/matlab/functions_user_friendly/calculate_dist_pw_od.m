function dist_pw_od_micron = calculate_dist_pw_od(pw_locations,ODCrtxPlt,pix2mm,show_fig)
% pw distance to OD borders 

[od_y,od_x] = find(edge(ODCrtxPlt > .5)); 
%[pw_locations, ~] = get_pinwheel_location_and_charge(ori_map,1);
dist_pw_od = zeros(size(pw_locations,1),1); 
for ii = 1 : size(pw_locations,1)
    pw_x = pw_locations(ii,1); 
    pw_y = pw_locations(ii,2); 
    
    dist_pw_od(ii) = min(sqrt((pw_x - od_x(:)).^2  + (pw_y - od_y(:)).^2)); 
end 

dist_pw_od_micron = dist_pw_od * pix2mm; 

mean_dist = mean(dist_pw_od_micron); 
std_dist = std(dist_pw_od_micron); 

%% 
if show_fig
    figure
    histogram(dist_pw_od_micron,10)
    xlabel('distance pw OD','fontsize',20)
    ylabel('frequency','fontsize',20)
    title(sprintf(' pinwheels OD distance (Microns) \n mean: %.2f \n std : %.2f',mean_dist,std_dist),'fontsize',20)
    set(gca,'tickdir','out','box','off')
end 