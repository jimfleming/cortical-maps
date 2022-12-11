function [New_ON_x,New_ON_y,New_OFF_x,New_OFF_y] = ONOFF_retina_grid(temp_xON,temp_yON,temp_xOFF,temp_yOFF,RfsizeRet,ON_j,OFF_j,min_onoff,max_onoff)

ON_x = zeros(size(temp_xON));
ON_y = zeros(size(temp_yON));
OFF_x = zeros(size(temp_xOFF));
OFF_y = zeros(size(temp_yOFF));

delta_y = zeros(length(temp_xON),1);
delta_x = zeros(length(temp_xON),1);
dist_onoff_old = zeros(length(temp_xON),1);
angle = zeros(length(temp_xON),1);
closest_off_ind  = zeros(length(temp_xON),1);
new_dist = zeros(length(temp_xON),1);
move_dist = zeros(length(temp_xON),1);
New_ON_x = zeros(length(temp_xON),1);
New_ON_y = zeros(length(temp_xON),1);

for ii = 1 : length(temp_xON)
    %ON_x(ii) = temp_xON(ii) + (rand(1)-0.5) * 2 * (ON_j*RfsizeRet);
    %ON_y(ii) = temp_yON(ii) + (rand(1)-0.5) * 2 * (ON_j*RfsizeRet);
    ON_x(ii) = temp_xON(ii) + (rand(1)-0.5) * 2 * (ON_j);
    ON_y(ii) = temp_yON(ii) + (rand(1)-0.5) * 2 * (ON_j);
end

for ii = 1 : length(temp_xOFF)
    %OFF_x(ii) = temp_xOFF(ii) + (rand(1)-0.5) * 2 * (OFF_j*RfsizeRet);
    %OFF_y(ii) = temp_yOFF(ii) + (rand(1)-0.5) * 2 * (OFF_j*RfsizeRet);
    OFF_x(ii) = temp_xOFF(ii) + (rand(1)-0.5) * 2 * (OFF_j);
    OFF_y(ii) = temp_yOFF(ii) + (rand(1)-0.5) * 2 * (OFF_j);
end

for on_n = 1 : length(temp_xON)
    disonofftemp = sqrt((OFF_x - ON_x(on_n)).^2 + (OFF_y - ON_y(on_n)).^2);
    [closes_dist_onoff,ind_closest_off] =  min(disonofftemp(:));
    delta_y(on_n) = ( -OFF_y(ind_closest_off) + ON_y(on_n) );
    delta_x(on_n) = ( -OFF_x(ind_closest_off) + ON_x(on_n) );
    dist_onoff_old(on_n) = closes_dist_onoff;
    angle(on_n) = atan( delta_y(on_n) / delta_x(on_n) ); % angle(on_n) * (180/pi)  
    closest_off_ind(on_n) = ind_closest_off;
end

max_dist = max(dist_onoff_old);
min_dist = min(dist_onoff_old);
m = (max_onoff - min_onoff)/(max_dist - min_dist);

New_OFF_x = OFF_x;
New_OFF_y = OFF_y;
for on_n = 1 : length(temp_xON)
    new_dist(on_n) = min_onoff + m *( dist_onoff_old(on_n) - min_dist );
    move_dist(on_n) = ( (new_dist(on_n) - dist_onoff_old(on_n)  )/2 ); % moving both ON and OFF 
    if delta_y(on_n) > 0 && delta_x(on_n) > 0
        project_angle_sin = abs( sin(angle(on_n)) );
        project_angle_cos =  abs( cos(angle(on_n)) );
    elseif delta_y(on_n) > 0 && delta_x(on_n) < 0
        project_angle_sin = abs( sin(angle(on_n)) );
        project_angle_cos =  -abs( cos(angle(on_n)) );
    elseif delta_y(on_n) < 0 && delta_x(on_n) < 0
        project_angle_sin = -abs( sin(angle(on_n)) );
        project_angle_cos =  -abs( cos(angle(on_n)) );
    elseif delta_y(on_n) < 0 && delta_x(on_n) > 0
        project_angle_sin = -abs( sin(angle(on_n)) );
        project_angle_cos =  abs( cos(angle(on_n)) );
    end
    New_OFF_x(closest_off_ind(on_n)) = OFF_x(closest_off_ind(on_n)) - (move_dist(on_n) * project_angle_cos) ;
    New_OFF_y(closest_off_ind(on_n)) = OFF_y(closest_off_ind(on_n)) - (move_dist(on_n) * project_angle_sin) ;
    New_ON_x(on_n) = ON_x(on_n) + (move_dist(on_n) * project_angle_cos) ;
    New_ON_y(on_n) = ON_y(on_n) + (move_dist(on_n) * project_angle_sin) ;
end

