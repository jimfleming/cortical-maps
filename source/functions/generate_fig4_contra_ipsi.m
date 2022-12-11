function generate_fig4_contra_ipsi(data_contra, data_ipsi, ONOFF_smoothed, boundary, NumRandPoint, electrode_row1, electrode_row2, num_rows_c, num_rows_d, show_fig)
% generate_fig4_linear_corr_jin(data_contra_mature,.1,30,100,200,20,25,[],1)
% the number of points in panel a and b is 30
% panel c : 600
% panel d : 25 lines 
% in data there are 35 points in fig4a, and 47 points in fig4b

    
    CV2 = data_contra.cv_map_interpolated;
    bound2 = round(size(CV2,1) * boundary);
    %range2 = bound2 : (size(CV2,1)-bound2);
    range2 = round(linspace(bound2, (size(CV2,1)-bound2), NumRandPoint)); 
    
    onoff_map = (ONOFF_smoothed > .5) * 1 + (ONOFF_smoothed <= .5) * -1; 
    
    %%
    figure;
    set(gcf , 'units', 'normalized','OuterPosition',[0 0 1.2609 1.2609]); %% [0 0 1.2609 1.2609]
    generate_panel_a_b(data_contra,data_ipsi,onoff_map,electrode_row1,range2,1,show_fig); % panel a
    generate_panel_a_b(data_contra,data_ipsi,onoff_map,electrode_row2,range2,2,show_fig);% panel b
    
    electrode_rows_c = round(linspace(bound2,(size(CV2,1)-bound2),num_rows_c));
    generate_panel_c(data_contra,data_ipsi,onoff_map,electrode_rows_c,range2,3,show_fig) ;
    
    electrode_rows = round(linspace(bound2,(size(CV2,1)-bound2),num_rows_d));
    generate_panel_d(data_contra,data_ipsi,onoff_map,electrode_rows,range2);
    
%    generate_panel_correlations(size_rect2,NumRandPoint,cvcrop2,lpicrop2,sf50crop2,lhicrop2,1,show_fig) % panel a
%    generate_panel_correlations(size_rect2,NumRandPoint,cvcrop2,lpicrop2,sf50crop2,lhicrop2,2,show_fig) % panel b 
%    generate_panel_c(size_rect2,num_points_panel_c,cvcrop2,lpicrop2,sf50crop2,lhicrop2,3,show_fig) % panel c
%    generate_panel_d(50,size_rect2,NumRandPoint,cvcrop2,lpicrop2,sf50crop2,lhicrop2);
    
    annotation('textbox',[0.1 0.82 0.08 0.08],'String','a','EdgeColor','none','fontsize',14)
    annotation('textbox',[0.52 0.82 0.08 0.08],'String','b','EdgeColor','none','fontsize',14)
    annotation('textbox',[0.1 0.35 0.08 0.08],'String','c','EdgeColor','none','fontsize',14)
    annotation('textbox',[0.52 0.35 0.08 0.08],'String','d','EdgeColor','none','fontsize',14)

    % text = sprintf()
    text = ['Dominant eye' ', Num random points : ' num2str(NumRandPoint) ', row a : ' num2str(electrode_row1) ...
             ', row b : ' num2str(electrode_row2) ', num rows c : ' num2str(num_rows_c) ', num rows d : ' num2str(num_rows_d)];
    annotation('textbox',[0.03 0.87 0.98 0.08],'String',text,'EdgeColor','none','fontsize',12)
end 


%% functions 


function [rr,pp] = generate_panel_a_b(data_contra,data_ipsi,onoff_map,electrode_row,range2,subplot_count,show_fig)
    
    CV2_contra = data_contra.cv_map_interpolated;
    LPI2_contra = data_contra.LPI_intepolated;
    SF50_2_contra = data_contra.sf50_map_interpolated;
    LHI2_contra = data_contra.LHI_interpolated;
    
    CV2_ipsi = data_ipsi.cv_map_interpolated;
    LPI2_ipsi = data_ipsi.LPI_intepolated;
    SF50_2_ipsi = data_ipsi.sf50_map_interpolated;
    LHI2_ipsi = data_ipsi.LHI_interpolated;
    
    range_contra = range2(onoff_map(electrode_row,range2) == 1);
    cv2_contra = CV2_contra(electrode_row,range_contra);
    lpi2_contra = LPI2_contra(electrode_row,range_contra);
    sf50_2_contra = SF50_2_contra(electrode_row,range_contra);
    lhi2_contra = LHI2_contra(electrode_row,range_contra);
     
    range_ipsi = range2(onoff_map(electrode_row,range2) == -1) ;
    cv2_ipsi = CV2_ipsi(electrode_row,range_ipsi);
    lpi2_ipsi = LPI2_ipsi(electrode_row,range_ipsi);
    sf50_2_ipsi = SF50_2_ipsi(electrode_row,range_ipsi);
    lhi2_ipsi = LHI2_ipsi(electrode_row,range_ipsi);
    
    cv2 = [cv2_contra cv2_ipsi]; 
    lpi2 = [lpi2_contra lpi2_ipsi]; 
    sf50_2 = [sf50_2_contra sf50_2_ipsi];
    lhi2 = [lhi2_contra lhi2_ipsi]; 
    
    rr = zeros(5,1);
    pp = zeros(5,1);
    [rr(1),pp(1)]=corr(cv2', sf50_2');
    [rr(2),pp(2)]=corr(cv2', lpi2');
    [rr(3),pp(3)]=corr(lhi2',sf50_2');
    [rr(4),pp(4)]=corr(lhi2', lpi2');
    [rr(5),pp(5)]=corr(lhi2', cv2');
    
    %%
    if show_fig == 1
        subplot_count = (subplot_count - 1) * 5;
        marker_size = 4;
        %   ****** Correlation for Interpolated version ****************
        subplot(2,10, 1 + subplot_count);
%         plot(cv2, sf50_2,'ko','markersize',marker_size); lsline;
        plot(cv2, log(sf50_2),'ko','markersize',marker_size); lsline;
        title(sprintf('  R = %.4f \nP = %.4f',rr(1),pp(1)));
        xlabel('CV'), ylabel('SF50');   xlim([0 1]); ylim([log(.25) log(2.5)]); axis square;
        ax = gca;
        set(ax, 'YTick', [log(0.5) log(1) log(2)]);
        ax.YTickLabel=({'0.5','1','2'});
        set(gca,'Tickdir','OUT','Box','OFF')
        %         ylim([ 0.25 , 2])
        
        subplot(2,10, 2 + subplot_count);
        plot(cv2, lpi2,'ko','markersize',marker_size); lsline;
        title(sprintf('  R = %.4f \nP = %.4f',rr(2),pp(2)));
        xlabel('CV'), ylabel('LPI');   xlim([0 1]); ylim([-0.05 1.05]); axis square;
        set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,10, 3 + subplot_count);
        plot(lhi2, log(sf50_2),'ko','markersize',marker_size); lsline;
        title(sprintf('  R = %.4f \nP = %.4f',rr(3),pp(3)));
        xlabel('LHI'), ylabel('SF50');   xlim([0 1]); ylim([log(.25) log(2.5)]); axis square;
        ax = gca;
        set(ax, 'YTick', [log(0.5) log(1) log(2)]);
        ax.YTickLabel=({'0.5','1','2'});
        set(gca,'Tickdir','OUT','Box','OFF')
%         ylim([ 0.25 , 2]) 
        
        subplot(2,10, 4 + subplot_count);
        plot(lhi2, lpi2,'ko','markersize',marker_size); lsline;
        title(sprintf('  R = %.4f \nP = %.4f',rr(4),pp(4)));
        xlabel('LHI'), ylabel('LPI');   xlim([0 1]); ylim([-0.05 1.05]);axis square;
        set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,10, 5 + subplot_count);
        plot(lhi2, cv2,'ko','markersize',marker_size); lsline;
        title(sprintf('  R = %.4f \nP = %.4f',rr(5),pp(5)));
        %title(['R = ', num2str(r), ' p = ', num2str(p)]);
        xlabel('LHI'), ylabel('CV'); ylim([0 1]);   xlim([0 1]);
        axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        % annotation('textbox',[0.03 0.87 0.08 0.08],'String',eyepref,'EdgeColor','none','fontsize',32)
    end
end 

function [rr,pp] = generate_panel_c(data_contra,data_ipsi,onoff_map,electrode_rows_c,range2_col,subplot_count,show_fig)
    
    CV2_contra = data_contra.cv_map_interpolated;
    LPI2_contra = data_contra.LPI_intepolated;
    SF50_2_contra = data_contra.sf50_map_interpolated;
    LHI2_contra = data_contra.LHI_interpolated;
    
    CV2_ipsi = data_ipsi.cv_map_interpolated;
    LPI2_ipsi = data_ipsi.LPI_intepolated;
    SF50_2_ipsi = data_ipsi.sf50_map_interpolated;
    LHI2_ipsi = data_ipsi.LHI_interpolated;
        
    n_points = length(electrode_rows_c) * length(range2_col); 
    cv2 = zeros(1,n_points); 
    lpi2 = zeros(1,n_points); 
    sf50_2 = zeros(1,n_points); 
    lhi2 = zeros(1,n_points); 
    cc = 0;
    ncontra = 0;
    nipsi = 0; 
    for ii = electrode_rows_c
        for jj = range2_col
            cc = cc + 1; 
            if onoff_map(ii,jj) == 1 
                cv2(cc) = CV2_contra(ii,jj);
                lpi2(cc) = LPI2_contra(ii,jj);
                sf50_2(cc) = SF50_2_contra(ii,jj);
                lhi2(cc) = LHI2_contra(ii,jj);
                ncontra = ncontra + 1 ;
            elseif onoff_map(ii,jj) == -1 
                cv2(cc) = CV2_ipsi(ii,jj);
                lpi2(cc) = LPI2_ipsi(ii,jj);
                sf50_2(cc) = SF50_2_ipsi(ii,jj);
                lhi2(cc) = LHI2_ipsi(ii,jj);
                nipsi = nipsi + 1 ;
            end 
            
        end 
    end 
    

%     [~,rect2,~,size_rect2] = remove_boundary(OriSFdata,boundary);
%     cvcrop2 = imcrop(CV2,rect2);
%     lpicrop2 = imcrop(LPI2,rect2);
%     sf50crop2 = imcrop(SF50_2,rect2);
%     lhicrop2 = imcrop(LHI2,rect2);

    % NumRandPoint = round(size_rect2 * boundary); 
    % selecting points randomly
    % rand_points_interpol = randperm(size_rect2, NumRandPoint);

%     cv2 = cvcrop2(rand_points_interpol);
%     lpi2 = lpicrop2(rand_points_interpol);
%     sf50_2 = sf50crop2(rand_points_interpol);
%     lhi2 = lhicrop2(rand_points_interpol);
    
%%
    rr = zeros(5,1);
    pp = zeros(5,1);
    [rr(1),pp(1)]=corr(cv2', sf50_2');
    [rr(2),pp(2)]=corr(cv2', lpi2');
    [rr(3),pp(3)]=corr(lhi2',sf50_2');
    [rr(4),pp(4)]=corr(lhi2', lpi2');
    [rr(5),pp(5)]=corr(lhi2', cv2');
    
    if show_fig == 1
        subplot_count = (subplot_count - 1) * 5;
        marker_size = 4;
        %   ****** Correlation for Interpolated version ****************
        subplot(2,10, 1 + subplot_count);
        plot(cv2, log(sf50_2),'ko','markersize',marker_size); lsline;
        title(sprintf('  R = %.4f \nP = %.4f',rr(1),pp(1)));
        xlabel('CV'), ylabel('SF50');   xlim([0 1]); ylim([log(.25) log(2.5)]); axis square;
        ax = gca;
        set(ax, 'YTick', [log(0.5) log(1) log(2)]);
        ax.YTickLabel=({'0.5','1','2'});
        set(gca,'Tickdir','OUT','Box','OFF')
%         ylim([ 0.25 , 2]) 
        
        subplot(2,10, 2 + subplot_count);
        plot(cv2, lpi2,'ko','markersize',marker_size); lsline;
        title(sprintf('  R = %.4f \nP = %.4f',rr(2),pp(2)));
        xlabel('CV'), ylabel('LPI');   xlim([0 1]); ylim([-0.05 1.05]); axis square;
        set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,10, 3 + subplot_count);
        plot(lhi2, log(sf50_2),'ko','markersize',marker_size); lsline;
        title(sprintf('  R = %.4f \nP = %.4f',rr(3),pp(3)));
        xlabel('LHI'), ylabel('SF50');   xlim([0 1]); ylim([log(.25) log(2.5)]); axis square;
        ax = gca;
        set(ax, 'YTick', [log(0.5) log(1) log(2)]);
        ax.YTickLabel=({'0.5','1','2'});
        set(gca,'Tickdir','OUT','Box','OFF')
%         ylim([ 0.25 , 2]) 
        
        subplot(2,10, 4 + subplot_count);
        plot(lhi2, lpi2,'ko','markersize',marker_size); lsline;
        title(sprintf('  R = %.4f \nP = %.4f',rr(4),pp(4)));
        xlabel('LHI'), ylabel('LPI');   xlim([0 1]); ylim([-0.05 1.05]); axis square;
        set(gca,'Tickdir','OUT','Box','OFF')
        
        subplot(2,10, 5 + subplot_count);
        plot(lhi2, cv2,'ko','markersize',marker_size); lsline;
        title(sprintf('  R = %.4f \nP = %.4f',rr(5),pp(5)));
        %title(['R = ', num2str(r), ' p = ', num2str(p)]);
        xlabel('LHI'), ylabel('CV'); ylim([0 1]);   xlim([0 1]);
        axis square; set(gca,'Tickdir','OUT','Box','OFF')
        
        % annotation('textbox',[0.03 0.87 0.08 0.08],'String',eyepref,'EdgeColor','none','fontsize',32)
    end
end 


function [rr_average,pp_average,mm_average] = generate_panel_d(data_contra,data_ipsi,onoff_map,electrode_rows,range2_columns)
    %%  Correlation CV/SF LHI/SF
    
    num_repeat = length(electrode_rows);
    rr_all = zeros(4,num_repeat);
    pp_all = zeros(4,num_repeat);
    mm_all = zeros(4,num_repeat);
    for ii = 1:length(electrode_rows)
        [rr,pp,mm] = generate_panel_d_fit(data_contra,data_ipsi,onoff_map,electrode_rows(ii),range2_columns);
        rr_all(:,ii) = rr; 
        pp_all(:,ii) = pp; 
        mm_all(:,ii) = mm; 
    end 
    
    rr_average = mean(rr_all,2); 
    pp_average = mean(pp_all,2); 
    mm_average = mean(pp_all,2); 
    
    m_cv_sf50 = mm_all(1,:); 
    m_cv_lpi = mm_all(2,:); 
    m_lhi_sf50 = mm_all(3,:); 
    m_lhi_lpi = mm_all(4,:); 
    
    %% figures 
    %{
    figure; clf;
    set(gcf,'position',[315         418        1314         420])
    subplot(1,2,1)

%     hold on 
%     subplot(2,12,[19:21])
    % subplot(4,11,[43,44])
    bin = .1;%0.05;
    xscale = -2.25:bin:1.5;
    [nx] = histc(m_cv_sf50,xscale);
    [ny] = histc(m_cv_lpi,xscale);
    plot(xscale+bin/2,nx,'color','k')
    hold on
    plot(xscale+bin/2,ny,'color',[0.5,0.5,0.5])
    ylabel('Numble of slopes'); 
    xlabel('Slope')
    box off
    % legend
    % axis square
%     text(0.2,16, 'CV / LPI','color',[0.5,0.5,0.5])
%     text(-1.5,16, 'CV / SF50r','color','k')
    text(-2,max(nx)-10, 'CV / LPI','color',[0.5,0.5,0.5])
    text(-2,max(nx)-15, 'CV / SF50r','color','k')
    % ylim([0 15])
    set(gca,'TickDir','out')
    hold off

    subplot(1,2,2)
    %subplot(2,12,[22:24])
    [nx] = histc(m_lhi_lpi,xscale);
    [ny] = histc(m_lhi_sf50 ,xscale);
    plot(xscale+bin/2,nx,'color','k')
    hold on
    plot(xscale+bin/2,ny,'color',[0.5,0.5,0.5])
    ylabel('Numble of slopes'); 
    xlabel('Slope')
    box off
    % legend
    %axis square
%     text(-1.1,16, 'LHI / LPI','color','k')
%     text(0.3,16, 'LHI / SF50r','color',[0.5,0.5,0.5])
    text(-2,max(ny)-15, 'LHI / LPI','color','k')
    text(-2,max(ny)-10, 'LHI / SF50r','color',[0.5,0.5,0.5])
    %ylim([0 15])
    set(gca,'TickDir','out')
    %}
    
    hold on
    pos1 = [0.55 0.22 0.15 0.12];
    subplot('Position',pos1)
    bin = .25;%0.05;
    xscale = -2.25:bin:1.5;
    [nx] = histc(m_cv_sf50,xscale);
    [ny] = histc(m_cv_lpi,xscale);
    plot( xscale + bin/2 , nx ,'color','k')
    hold on
    plot(xscale+bin/2,ny,'color',[0.5,0.5,0.5])
    ylabel('Numble of slopes'); 
    xlabel('Slope')
    box off
    % legend
    %axis square
%     text(0.2,16, 'CV / LPI','color',[0.5,0.5,0.5])
%     text(-1.5,16, 'CV / SF50r','color','k')
    text(-1.8, 12 , 'CV / LPI','color',[0.5,0.5,0.5])
    text(-1.8, 14 , 'CV / SF50','color','k')
    %ylim([0 15])
    set(gca,'TickDir','out')
    hold off
    xlim([ -2, 2]) 
    ylim([ 0 , 15]) 
    
    pos2 = [0.75 0.22 0.15 0.12];
    subplot('Position',pos2)
    [nx] = histc(m_lhi_lpi,xscale);
    [ny] = histc(m_lhi_sf50 ,xscale);
    plot(xscale+bin/2,nx,'color','k')
    hold on
    plot(xscale+bin/2,ny,'color',[0.5,0.5,0.5])
    ylabel('Numble of slopes'); 
    xlabel('Slope')
    box off
    % legend
    %axis square
%     text(-1.1,16, 'LHI / LPI','color','k')
%     text(0.3,16, 'LHI / SF50r','color',[0.5,0.5,0.5])
    text(-1.8,12, 'LHI / LPI','color','k')
    text(-1.8,14, 'LHI / SF50','color',[0.5,0.5,0.5])
    %ylim([0 15])
    set(gca,'TickDir','out')
    xlim([ -2, 2]) 
    ylim([ 0 , 15]) 
end 

function [rr,pp,mm] = generate_panel_d_fit(data_contra,data_ipsi,onoff_map,electrode_row,range2)
    
    CV2_contra = data_contra.cv_map_interpolated;
    LPI2_contra = data_contra.LPI_intepolated;
    SF50_2_contra = data_contra.sf50_map_interpolated;
    LHI2_contra = data_contra.LHI_interpolated;
    
    CV2_ipsi = data_ipsi.cv_map_interpolated;
    LPI2_ipsi = data_ipsi.LPI_intepolated;
    SF50_2_ipsi = data_ipsi.sf50_map_interpolated;
    LHI2_ipsi = data_ipsi.LHI_interpolated;
    
    range_contra = range2(onoff_map(electrode_row,range2) == 1);
    cv2_contra = CV2_contra(electrode_row,range_contra);
    lpi2_contra = LPI2_contra(electrode_row,range_contra);
    sf50_2_contra = SF50_2_contra(electrode_row,range_contra);
    lhi2_contra = LHI2_contra(electrode_row,range_contra);
     
    range_ipsi = range2(onoff_map(electrode_row,range2) == -1) ;
    cv2_ipsi = CV2_ipsi(electrode_row,range_ipsi);
    lpi2_ipsi = LPI2_ipsi(electrode_row,range_ipsi);
    sf50_2_ipsi = SF50_2_ipsi(electrode_row,range_ipsi);
    lhi2_ipsi = LHI2_ipsi(electrode_row,range_ipsi);
    
    cv2 = [cv2_contra cv2_ipsi]; 
    lpi2 = [lpi2_contra lpi2_ipsi]; 
    sf50_2 = [sf50_2_contra sf50_2_ipsi];
    lhi2 = [lhi2_contra lhi2_ipsi]; 
    
%     cv2 = CV2(electrode_row,range2);
%     lpi2 = LPI2(electrode_row,range2);
%     sf50_2 = SF50_2(electrode_row,range2);
%     lhi2 = LHI2(electrode_row,range2);

    rr = zeros(4,1);
    pp = zeros(4,1);
    mm = zeros(4,1);
    [rr(1),pp(1)]=corr(cv2', sf50_2');
    [rr(2),pp(2)]=corr(cv2', lpi2');
    [rr(3),pp(3)]=corr(lhi2',sf50_2');
    [rr(4),pp(4)]=corr(lhi2', lpi2');
    
    pp1 =polyfit(cv2', sf50_2', 1);
    pp2 =polyfit(cv2', lpi2', 1);
    pp3 =polyfit(lhi2',sf50_2', 1);
    pp4 =polyfit(lhi2', lpi2', 1);
    
    mm(1) = pp1(1); 
    mm(2) = pp2(1); 
    mm(3) = pp3(1); 
    mm(4) = pp4(1); 
end 

function [rect1,rect2,size_rect1,size_rect2] = remove_boundary(OriSFdata,boundary)
    % remove boundary
    LHI = OriSFdata.LHI_map; 
    LHI_interpolated = OriSFdata.LHI_interpolated; 
    bound1 = round(boundary * size(LHI,1));
    bound2 = round(boundary * size(LHI_interpolated,1));
    rect1 = [bound1 bound1 size(LHI,1)-bound1 size(LHI,2)-bound1];
    rect2 = [bound2 bound2 size(LHI_interpolated,1)-bound2 size(LHI_interpolated,2)-bound2];
    
    size_rect1 = (size(LHI,1)- 2*bound1)^2 ;
    size_rect2 = (size(LHI_interpolated,1)- 2*bound2)^2 ;
end
