function [appdata] = make_retina_app2(app, appdata, show_fig)
    
   distAPON = appdata.retina.distAPON;
   distMLON = appdata.retina.distMLON;
   distAPOFF = appdata.retina.distAPOFF;
   distMLOFF = appdata.retina.distMLOFF; % -------------------------??
   jitterON = appdata.retina.jitterON;
   jitterOFF = appdata.retina.jitterOFF;
   CrtxLength = appdata.CrtxLength;
   MinDistONOFF = appdata.retina.MinDistONOFF;
   MaxDistONOFF = appdata.retina.MaxDistONOFF;
   pix2mic = appdata.pix2mic;
   rng_trial = appdata.rng_trial;

   % the only difference with make_retina5 is adding rng_trial 
    
    %ML_AP_rate = 1 ;           % ratio of the number of cells in medio-lateral (Y) to anterio-posterior (X) axis (distMLON = distAPON * ML_AP_rate)
    %distMLON = distAPON * ML_AP_rate; 
    ML_AP_rate = distMLON / distAPON; 
    %distAPOFF = distAPON; 
    
    num_cell = CrtxLength ^ 2; 
    Width = sqrt((num_cell/4) * (distAPON^2) * ML_AP_rate) ; % 4 because there are two eyes and two polarities 
    Width = ceil(Width); 
    
%     VisSpaceBoundX = floor((distAPON*CrtxLength)/4); % devided by 4(2*2), 2 for boundary is starting from [-VisSpaceBoundX to VisSpaceBoundX]
%     VisSpaceBoundY = floor((distMLON*CrtxLength)/4); %                    2 because there are two retinas 
%     GridONXstep = distAPON  ;
%     GridONYstep = distMLON  ;
    
%     range_x = GridONXstep : GridONXstep : (GridONXstep * Width); 
%     range_y = GridONXstep : GridONYstep : (GridONXstep * Width); 
    num_cell_x = ceil(Width / distAPON) ;
    num_cell_y = ceil(Width / distMLON) ;
    
    % to make the assymetry for macaque retina 
    if strcmp(app.retina_switch.Value, 'symmetrical')
        range_x = linspace(distAPON,Width * (1),num_cell_x); 
    elseif strcmp(app.retina_switch.Value, 'asymmetrical')
        range_x = linspace(distAPON,Width * (3/4),num_cell_x); 
    end 
    
    % range_y = linspace(distAPON,Width*(3/4),num_cell_y); 
    range_y = linspace(distAPON,Width * (1),num_cell_y); 
    
    [ONxGridRetina , ONyGridRetina] = meshgrid( range_x , range_y );
    OFFxGridRetina = ONxGridRetina; 
    OFFyGridRetina = ONyGridRetina; 
    

    ON_j = (distAPON) * jitterON; 
    OFF_j = (distAPOFF) * jitterOFF;

    [RetONxContra, RetONyContra, RetOFFxContra, RetOFFyContra] = ONOFF_retina_grid2(ONxGridRetina(:), ONyGridRetina(:), OFFxGridRetina(:), OFFyGridRetina(:), [], ON_j, OFF_j, MinDistONOFF, MaxDistONOFF, rng_trial);
    [RetONxIpsi, RetONyIpsi, RetOFFxIpsi, RetOFFyIpsi] = ONOFF_retina_grid2(ONxGridRetina(:), ONyGridRetina(:), OFFxGridRetina(:), OFFyGridRetina(:), [], ON_j, OFF_j, MinDistONOFF, MaxDistONOFF, rng_trial+2); % rng_trial + 1 to be different with the other eye

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
    
    if show_fig == 1 
        num_bin = 10;
        measure_dist(RetONxContra, RetONyContra, RetOFFxContra, RetOFFyContra, rf_space_x, num_bin,...
            jitterON, jitterOFF, MinDistONOFF, MaxDistONOFF, pix2mic, 'Contra', 1);

    end 
    
%     show_RGC_distribution_app(RetONxContra, RetONyContra, RetOFFxContra, RetOFFyContra, rf_space_x, app.UIAxes_RGC_contra)
%     show_RGC_distribution_app(RetONxIpsi, RetONyIpsi, RetOFFxIpsi, RetOFFyIpsi, rf_space_x, app.UIAxes_RGC_ipsi)

    appdata.RetinaAllX = RetinaAllX;
    appdata.RetinaAllY = RetinaAllY;
    appdata.RetinaAllOD = RetinaAllOD;
    appdata.RetinaAllONOFF = RetinaAllONOFF;
    appdata.TotalRetinaCells = TotalRetinaCells;
    appdata.rf_space_x = rf_space_x;
    appdata.rf_space_y = rf_space_y;
    
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