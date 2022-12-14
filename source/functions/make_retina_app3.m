function [appdata] = make_retina_app3(appdata, show_fig)
    jitter_seed = appdata.rng_trial; % Jitter seed number
    
    %%
    ADF = appdata.retina.aff_density_factor; % afferent density factor 
    
    %%
    distAPON = appdata.retina.distAPON;
    distMLON = appdata.retina.distMLON;
    distAPOFF = appdata.retina.distAPOFF;
    distMLOFF = appdata.retina.distMLOFF; 
   
    dist_on_x  = distAPON;  % def val 5; % Anterior Posterior
    dist_on_y  = distMLON;  %5;         % Medio Lateral
    dist_off_x = distAPOFF; %5
    dist_off_y = distMLOFF; %5
    dist_min = min([dist_on_x, dist_on_y, dist_off_x, dist_off_y]);
    appdata.retina.dist_min_RGC = dist_min; 

    jitterONx = appdata.retina.jitterONx;
    jitterONy = appdata.retina.jitterONy;
    jitterOFFx = appdata.retina.jitterOFFx;
    jitterOFFy = appdata.retina.jitterOFFy; 
   
    ON_jx  = (dist_on_x) * jitterONx;
    ON_jy  = (dist_on_y) * jitterONy;
    OFF_jx = (dist_off_x) * jitterOFFx;
    OFF_jy = (dist_off_y) * jitterOFFy;
   
    %%
    CrtxLength = appdata.CrtxLength;
    MinDistONOFF = appdata.retina.MinDistONOFF;
    MaxDistONOFF = appdata.retina.MaxDistONOFF;
    pix2mic = appdata.pix2mic;
 
    %% Measure width of RF space 
    num_cell = CrtxLength^2 ;
    Width_rf_space_onx = ceil(sqrt( (num_cell/4) * (dist_on_x^2) )); % 4 because there are two eyes and two polarities
    Width_rf_space_ony = ceil(sqrt( (num_cell/4) * (dist_on_y^2) )); % 4 because there are two eyes and two polarities
    Width_rf_space_offx = ceil(sqrt( (num_cell/4) * (dist_off_x^2) )); % 4 because there are two eyes and two polarities
    Width_rf_space_offy = ceil(sqrt( (num_cell/4) * (dist_off_y^2) )); % 4 because there are two eyes and two polarities

    %% Number of Cells each eye 
    OFF_ON_ratio = 1;         % Number of OFF to ON cell ratio
    Contra_ipsi_ratio = 1;    % Number of Contra to Ipsi cell ratio

    num_cell_contra = (num_cell / (Contra_ipsi_ratio + 1)) * Contra_ipsi_ratio ;
    num_cell_contra = ceil(num_cell_contra);
    num_cell_ipsi = num_cell - num_cell_contra;

    %% Contra
    num_OFF_contra = (num_cell_contra / (OFF_ON_ratio + 1)) * OFF_ON_ratio ; 
    num_OFF_contra = ceil(num_OFF_contra); 
    num_ON_contra = num_cell_contra - num_OFF_contra; 

    % on 
    num_ONx_contra = sqrt(num_ON_contra); 
    num_ONy_contra = num_ONx_contra; 
    max_rangeX_ON_contra = round(Width_rf_space_onx); 
    max_rangeY_ON_contra = Width_rf_space_ony; 
    range_ONx_contra = linspace(dist_on_x, max_rangeX_ON_contra, num_ONx_contra);
    range_ONy_contra = linspace(dist_on_y, max_rangeY_ON_contra, num_ONy_contra);

    % off 
    num_OFFx_contra = sqrt(num_OFF_contra); 
    num_OFFy_contra = num_OFFx_contra; 
    max_rangeX_OFF_contra = round(Width_rf_space_offx); 
    max_rangeY_OFF_contra = Width_rf_space_offy; 
    range_OFFx_contra = linspace(dist_off_x, max_rangeX_OFF_contra, num_OFFx_contra);
    range_OFFy_contra = linspace(dist_off_y, max_rangeY_OFF_contra, num_OFFy_contra);

    %% ipsi 
    num_OFF_ipsi = (num_cell_ipsi / (OFF_ON_ratio + 1)) * OFF_ON_ratio ; 
    num_OFF_ipsi = ceil(num_OFF_ipsi); 
    num_ON_ipsi = num_cell_ipsi - num_OFF_ipsi; 

    % on 
    num_ONx_ipsi = sqrt(num_ON_ipsi); 
    num_ONy_ipsi = num_ONx_ipsi; 
    max_rangeX_ON_ipsi = round(Width_rf_space_onx); 
    max_rangeY_ON_ipsi = Width_rf_space_ony; 
    range_ONx_ipsi = linspace(dist_on_x, max_rangeX_ON_ipsi, num_ONx_ipsi);
    range_ONy_ipsi = linspace(dist_on_y, max_rangeY_ON_ipsi, num_ONy_ipsi);

    % off 
    num_OFFx_ipsi = sqrt(num_OFF_ipsi); 
    num_OFFy_ipsi = num_OFFx_ipsi; 
    max_rangeX_OFF_ipsi = round(Width_rf_space_offx); 
    max_rangeY_OFF_ipsi = Width_rf_space_offy; 
    range_OFFx_ipsi = linspace(dist_off_x, max_rangeX_OFF_ipsi, num_OFFx_ipsi);
    range_OFFy_ipsi = linspace(dist_off_y, max_rangeY_OFF_ipsi, num_OFFy_ipsi);

    %% 
    [ONx_Grid_contra , ONy_Grid_contra] = meshgrid( range_ONx_contra , range_ONy_contra );
    [OFFx_Grid_contra , OFFy_Grid_contra] = meshgrid( range_OFFx_contra , range_OFFy_contra );
    [RetONxContra, RetONyContra, RetOFFxContra, RetOFFyContra] = ONOFF_retina_grid3(ONx_Grid_contra(:), ONy_Grid_contra(:), OFFx_Grid_contra(:), OFFy_Grid_contra(:), [], ON_jx, ON_jy, OFF_jx, OFF_jy, MinDistONOFF, MaxDistONOFF, jitter_seed);

    [ONx_Grid_ipsi , ONy_Grid_ipsi] = meshgrid( range_ONx_ipsi , range_ONy_ipsi );
    [OFFx_Grid_ipsi , OFFy_Grid_ipsi] = meshgrid( range_OFFx_ipsi , range_OFFy_ipsi );
    [RetONxIpsi,RetONyIpsi,RetOFFxIpsi,RetOFFyIpsi] = ONOFF_retina_grid3(ONx_Grid_ipsi(:), ONy_Grid_ipsi(:), OFFx_Grid_ipsi(:), OFFy_Grid_ipsi(:), [], ON_jx, ON_jy, OFF_jx, OFF_jy, MinDistONOFF, MaxDistONOFF, jitter_seed+2); % rng_trial + 1 to be different with the other eye

    %%
    TotalRetinaCells = length(RetONxContra) + length(RetOFFxContra) + length(RetONxIpsi) + length(RetOFFxIpsi); 
    RetinaAllX = cat(1,RetONxContra,RetOFFxContra,RetONxIpsi,RetOFFxIpsi); 
    RetinaAllY = cat(1,RetONyContra,RetOFFyContra,RetONyIpsi,RetOFFyIpsi); 
    RetinaAllOD = cat(1, ones(length(RetONxContra),1), ones(length(RetOFFxContra),1), -1*ones(length(RetONxIpsi),1), -1*ones(length(RetOFFxIpsi), 1));
    RetinaAllONOFF = cat(1, ones(length(RetONxContra),1), -1*ones(length(RetOFFxContra),1), ones(length(RetONxIpsi),1), -1*ones(length(RetOFFxIpsi), 1));

    BoundXL = floor(min(RetinaAllX(:))) - 2;
    BoundXR = ceil(max(RetinaAllX(:))) + 2;
    BoundYL = floor(min(RetinaAllY(:))) - 2;
    BoundYR = ceil(max(RetinaAllY(:))) + 2;
    
    BoundL = min(BoundXL,BoundYL);
    BoundU = max(BoundXR,BoundYR);
    
    % make rf space 
    [rf_space_x, rf_space_y] = meshgrid(BoundL:BoundU, BoundL:BoundU);
    
    appdata.RetONxContra  = RetONxContra;
    appdata.RetONyContra  = RetONyContra;
    appdata.RetOFFxContra = RetOFFxContra;
    appdata.RetOFFyContra = RetOFFyContra;
    appdata.rf_space_x = rf_space_x;

    %% afferent sorting filter parameter based on RGC distribution
    XY_ratio = dist_on_x / dist_on_y ; 

    CenterRadiusX = (5 / dist_min) * 5;  % 5 is the base value we used in the 

    if CenterRadiusX < 3
        % 3 arbitrary number, I could add more options for RGC distance and the
        % aff_sampling_density would be lower than 1 automatically,
        % but the running time increases for larger initial RGC
        % distance because it increases the visual space size
        CenterRadiusX = .5 ;
    end

    CenterRadiusX = CenterRadiusX * ADF; 
    CenterRadiusY = CenterRadiusX; 
    large_surround_filter_size = CenterRadiusX*2; 
    if XY_ratio == 1 
        SurroundRadiusX = large_surround_filter_size; 
        SurroundRadiusY = large_surround_filter_size;
    elseif XY_ratio > 1 
        SurroundRadiusX = CenterRadiusX; 
        SurroundRadiusY = large_surround_filter_size;
    elseif XY_ratio < 1 
        SurroundRadiusX = large_surround_filter_size; 
        SurroundRadiusY = CenterRadiusX;
    end

    %% 
    appdata.RetinaAllX = RetinaAllX;
    appdata.RetinaAllY = RetinaAllY;
    appdata.RetinaAllOD    = RetinaAllOD;
    appdata.RetinaAllONOFF = RetinaAllONOFF;
    appdata.TotalRetinaCells = TotalRetinaCells;
    appdata.rf_space_x = rf_space_x;
    appdata.rf_space_y = rf_space_y;
    
    appdata.CenterRadiusX = CenterRadiusX;
    appdata.CenterRadiusY = CenterRadiusY;
    appdata.SurroundRadiusX = SurroundRadiusX;
    appdata.SurroundRadiusY = SurroundRadiusY;
end 

function [New_ON_x, New_ON_y, New_OFF_x, New_OFF_y] = ONOFF_retina_grid3(temp_xON, temp_yON, temp_xOFF, temp_yOFF, RfsizeRet, ON_jx, ON_jy, OFF_jx, OFF_jy, ...
    min_onoff, max_onoff, rng_trial)
    % rng_trial = 1; 
    rng(rng_trial)
    
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

    % adding jitter 
    for ii = 1 : length(temp_xON)
        ON_x(ii) = temp_xON(ii) + (rand(1)-0.5) * 2 * (ON_jx);
        ON_y(ii) = temp_yON(ii) + (rand(1)-0.5) * 2 * (ON_jy);
    end

    for ii = 1 : length(temp_xOFF)
        OFF_x(ii) = temp_xOFF(ii) + (rand(1)-0.5) * 2 * (OFF_jx);
        OFF_y(ii) = temp_yOFF(ii) + (rand(1)-0.5) * 2 * (OFF_jy);
    end

    % measuring the distance of the closest RGCs with different polairty 
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
end 
