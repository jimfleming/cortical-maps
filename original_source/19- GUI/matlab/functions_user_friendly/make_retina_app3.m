function [appdata] = make_retina_app3(appdata, show_fig)
    
    jitter_seed = appdata.rng_trial; % Jitter seed number
    
    %%
    ADF = appdata.retina.aff_density_factor; % afferent density factor 
    % rf_sd_RGC = ADF * 5; % default val : 5; 
    
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
%    dist_max = max([dist_on_x, dist_on_y, dist_off_x, dist_off_y]);

   %%
%    jitterON = appdata.retina.jitterON;
%    jitterOFF = appdata.retina.jitterOFF;

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
    %Width_rf_space = sqrt( (num_cell/4) * (dist_on_x^2) ); % 4 because there are two eyes and two polarities
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
    
%     if show_fig == 1 
%         num_bin = 10;
%         measure_dist(RetONxContra, RetONyContra, RetOFFxContra, RetOFFyContra, rf_space_x, num_bin,...
%             jitterON, jitterOFF, MinDistONOFF, MaxDistONOFF, pix2mic, 'Contra', 1);
% 
%     end 
    
%     show_RGC_distribution_app(RetONxContra, RetONyContra, RetOFFxContra, RetOFFyContra, rf_space_x, app.UIAxes_RGC_contra)
%     show_RGC_distribution_app(RetONxIpsi, RetONyIpsi, RetOFFxIpsi, RetOFFyIpsi, rf_space_x, app.UIAxes_RGC_ipsi)

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

 

function [New_ON_x,New_ON_y,New_OFF_x,New_OFF_y] = ONOFF_retina_grid2(temp_xON,temp_yON,temp_xOFF,temp_yOFF,RfsizeRet,ON_j,OFF_j,min_onoff,max_onoff,rng_trial)
    %rng_trial = 1; 
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
        %ON_x(ii) = temp_xON(ii) + (rand(1)-0.5) * 2 * (ON_j*RfsizeRet);
        %ON_y(ii) = temp_yON(ii) + (rand(1)-0.5) * 2 * (ON_j*RfsizeRet);
        ON_x(ii) = temp_xON(ii) + (rand(1)-0.5) * 2 * (ON_jx);
        ON_y(ii) = temp_yON(ii) + (rand(1)-0.5) * 2 * (ON_jy);
    end

    for ii = 1 : length(temp_xOFF)
        %OFF_x(ii) = temp_xOFF(ii) + (rand(1)-0.5) * 2 * (OFF_j*RfsizeRet);
        %OFF_y(ii) = temp_yOFF(ii) + (rand(1)-0.5) * 2 * (OFF_j*RfsizeRet);
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

function measure_dist(RetONx,RetONy,RetOFFx,RetOFFy,rf_space_x,num_bin,...
                jitterON,jitterOFF,MinDistONOFF,MaxDistONOFF,pix2mic,eye,show_fig)
    %% Measuring RGCs Distance
    on_center_x = RetONx(:);
    on_center_y = RetONy(:);
    off_center_x = RetOFFx(:);
    off_center_y = RetOFFy(:);

    OnOnDist = zeros(length(on_center_x),1);
    OffOffDist = zeros(length(off_center_x),1);
    OnOffDist = zeros(length(on_center_x),1);

    for aa = 1 : length(on_center_x)
        temp_dis_on_on = sqrt((on_center_x(aa) - on_center_x).^2+(on_center_y(aa) - on_center_y).^2);
        SortedDistONON =  sort(temp_dis_on_on);
        OnOnDist(aa) = SortedDistONON(2);%  The first index is the point itself
    end

    for bb = 1 : length(off_center_x)
        temp_dis_off_off = sqrt((off_center_x(bb) - off_center_x).^2+(off_center_y(bb) - off_center_y).^2);
        SortedDistOffOff =  sort(temp_dis_off_off);
        OffOffDist(bb) = SortedDistOffOff(2);%  The first index is the point itself
    end

    for cc = 1 : length(on_center_x)
        temp_dis_on_off = sqrt((on_center_x(cc) - off_center_x).^2+(on_center_y(cc) - off_center_y).^2);
        SortedDistOnOff =  sort(temp_dis_on_off);
        OnOffDist(cc) = SortedDistOnOff(1);
    end
    
    if show_fig == 1 
        figure,
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        show_RGC_distribution(RetONx,RetONy,RetOFFx,RetOFFy,rf_space_x)
        ONONDist_mic = OnOnDist * pix2mic ; 
        OffOffDist_mic = OffOffDist * pix2mic ; 
        OnOffDist_mic = OnOffDist * pix2mic ; 
        show_RGC_distance(ONONDist_mic,OffOffDist_mic,OnOffDist_mic,num_bin)
        
        title_text = sprintf([eye '\n jitter ON = %.2f \n jitter OFF = %.2f \n min distance ONOFF = %.2f pix (%.2f micron) \n max distance ONOFF = %.2f pix (%.2f micron)'], ...
            jitterON,jitterOFF,MinDistONOFF,MinDistONOFF* pix2mic,MaxDistONOFF,MaxDistONOFF* pix2mic); 
        annotation('textbox',[0.1 0.88 0.1 0.1],'String',title_text,'EdgeColor','none','fontsize',20)

        %sgtitle(sprintf(eye,'fontsize',20)
    end 
    
end 

function show_RGC_distribution(RetONx,RetONy,RetOFFx,RetOFFy,rf_space_x)
    % Plotting ON OFF RGCs on Retinotopy map
        axis_lim = size(rf_space_x,1);
        MarkerSize = 20;

        subplot(1,6,[1,2]),imshow(ones(size(rf_space_x,1)))
        hold on
        plot( RetONx(:) , RetONy(:) ,'r.' ,'MarkerSize', MarkerSize)
        plot( RetOFFx(:), RetOFFy(:),'b.' ,'MarkerSize', MarkerSize)
        %title(eye,'fontsize',20)
        xlim([0 axis_lim]),ylim([0 axis_lim])
        axis on 
        
        subplot(1,6,[3,4]),imshow(ones(size(rf_space_x,1)))
        hold on
        MarkerSize2 = 75;
        plot( RetONx(:) , RetONy(:) ,'r.' ,'MarkerSize', MarkerSize2)
        plot( RetOFFx(:), RetOFFy(:),'b.' ,'MarkerSize', MarkerSize2)
        %title(eye,'fontsize',20)
        xlim([2.5 32.5]),ylim([2.5 32.5])
        axis on 
end 

function show_RGC_distance(onon_dist,offoff_dist,onoff_dist,num_bin)
    xx = 0:0.1:150; 

    subplot(3,6,6), cla
    %     [h,edge] = hist(onon_dist,num_bin);
    %     plot(edge,h,'r'), 
    %hist_onon = histogram(onon_dist,num_bin,'Normalization','pdf');
    hist_onon = histogram(onon_dist,'BinEdges',[0 : 10 : 150],'Normalization','pdf');
    hist_onon.FaceColor = 'r';
    hold on, 
    mu_onon = mean(onon_dist);
    sig_onon = std(onon_dist);
    f_onon = exp(-(xx-mu_onon).^2./(2*sig_onon^2))./(sig_onon*sqrt(2*pi));
    plot(xx,f_onon,'k','linewidth',1.5)
    title(sprintf('M=%.2f, SD=%.2f',mu_onon,sig_onon),'fontsize',20)
    ylabel('ON-ON  ','FontSize',20), set(gca,'box','off','Tickdir','out'), axis square
    %xlim([min(onon_dist(:)) max(onon_dist(:))])
    ylim([0 0.05])
    xlim([min(xx(:)) max(xx(:))])
    
    subplot(3,6,12), cla
%     [h,edge] = hist(off_dist,num_bin);
%     plot(edge,h,'b'),
    % hist_offoff = histogram(offoff_dist,num_bin,'Normalization','pdf');
    hist_offoff = histogram(offoff_dist,'BinEdges',[0 : 10 : 150],'Normalization','pdf');
    hist_offoff.FaceColor = 'b';
    hold on, 
    mu_offoff = mean(offoff_dist);
    sig_offoff = std(offoff_dist);
    f_offoff = exp(-(xx-mu_offoff).^2./(2*sig_offoff^2))./(sig_offoff*sqrt(2*pi));
    plot(xx,f_offoff,'k','linewidth',1.5)
    title(sprintf('M=%.2f, SD=%.2f',mu_offoff,sig_offoff),'fontsize',20)
    ylabel('OFF-OFF  ','FontSize',20), set(gca,'box','off','Tickdir','out'), axis square
%     xlim([min(offoff_dist(:)) max(offoff_dist(:))])
    ylim([0 0.05])
    xlim([min(xx(:)) max(xx(:))])
    
    subplot(3,6,18), cla
%     [h,edge] = hist(on_off_dist,num_bin);
%     plot(edge,h,'k'), 
    %hist_onoff = histogram(onoff_dist,num_bin,'Normalization','pdf');
    hist_onoff = histogram(onoff_dist,'BinEdges',[0 : 10 : 150],'Normalization','pdf');
    hist_onoff.FaceColor = 'w';
    hold on, 
    mu_onoff = mean(onoff_dist);
    sig_onoff = std(onoff_dist);
    f_onoff = exp(-(xx-mu_onoff).^2./(2*sig_onoff^2))./(sig_onoff*sqrt(2*pi));
    plot(xx,f_onoff,'k','linewidth',1.5)
    title(sprintf('M=%.2f, SD=%.2f',mean(onoff_dist),std(onoff_dist)),'fontsize',20)
    ylabel('ON-OFF  ','FontSize',20), set(gca,'box','off','Tickdir','out'), axis square
%     xlim([min(onoff_dist(:)) max(onoff_dist(:))])
    ylim([0 0.05])
    xlim([min(xx(:)) max(xx(:))])
end


  
function overlap_hist = measure_overlap(dist_vec,R)
    %   Measuring the overlap by formula 
    overlap_hist = (2*R^2*acos(dist_vec./(2*R))-dist_vec.*sqrt(R^2-dist_vec.^2/4))/(pi*R^2)*100;
end 

function show_RGC_distribution_app(RetONx,RetONy,RetOFFx,RetOFFy,rf_space_x, ax)
    % Plotting ON OFF RGCs on Retinotopy map
%         axis_lim = size(rf_space_x,1);
%         MarkerSize = 20;

%         subplot(1,6,[1,2]),imshow(ones(size(rf_space_x,1)))
%         hold on
%         plot( RetONx(:) , RetONy(:) ,'r.' ,'MarkerSize', MarkerSize)
%         plot( RetOFFx(:), RetOFFy(:),'b.' ,'MarkerSize', MarkerSize)
%         %title(eye,'fontsize',20)
%         xlim([0 axis_lim]),ylim([0 axis_lim])
%         axis on 
        %ax = app.UIAxes_RGC_contra
        ax.Visible = 'on'; 
        cla(ax)
%         imshow( ones(size(rf_space_x,1)), 'Parent', ax)
%         hold(ax, 'on')
        MarkerSize2 = 20;
        plot( RetONx(:) , RetONy(:) ,'r.' ,'MarkerSize', MarkerSize2, 'Parent', ax)
        hold(ax, 'on')
        plot(ax, RetOFFx(:), RetOFFy(:),'b.' ,'MarkerSize', MarkerSize2, 'Parent', ax)
        %title(eye,'fontsize',20)
        xlim(ax, [2.5 32.5]),ylim(ax, [2.5 32.5])
        axis(ax, 'off') 
        axis(ax, 'square')
       % box(ax, 'on')
        set(ax,'ydir','reverse')
        pause(0.001)
end 