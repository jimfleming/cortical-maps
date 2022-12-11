function generate_fig3c_curves_modified(OriSFdata_contra,OriSFdata_ipsi,ODCrtxPlt_smooth,electrode_position,eye_pref,rng_trial)
    % Simulate data of figure 3a with the selected electrode positions
    % The ori/sf curves are modified for presentation
    % The SF curve is estimated with more number of points 
    % The ori tuning curve is modified (smoother and in case there are two peaks the main peak is kept)
    % 
    % Simulate data of figure 3c with the selected electrode positions
    % OriSFdata     :   All the data, LHI, LPI, SF tuning curve, Ori tuning curve
    % ODCrtxPlt_smooth      :      ocular dominace for (ODI), in future ODI should be meausred from RFs
    % electrodePos          :      Selected series of positions : round([linspace(start_row,end_row,5)' linspace(start_col,end_col,5)'])
    % example : norm([51-43,43-21]) ~= 8, 40 to 48 == 500 microns
    %
    % eye_pref              :      the preferred eye for each location
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n_interpolation = 1000; % number of points to interpolate SF curves 
    sf_lSamp = 1:50 ;       % in pixel 
    debug = 0;              % show the selected part of the cortex 
    
    rng(rng_trial)
    %%
    LHI_contra = OriSFdata_contra.LHI_map;
    CV_contra = OriSFdata_contra.CV_map;
    ori_hist_contra = OriSFdata_contra.OriTunPolarPlot;
    ori_bin_contra = OriSFdata_contra.angles;
    % sf_hist_contra = OriSFdata_contra.allsfcurve_dog;
    % allsfcurve_dog is interpolated with less point, the raw data is used
    % here to interpolate with more number of points 
    sf_tuning_contra = OriSFdata_contra.allsfcurve; %  allsfcurve_dog
    sfpref_contra = OriSFdata_contra.SF_map;
    sf_bin_deg_contra = OriSFdata_contra.sf_bin_deg;
    SF50_contra = OriSFdata_contra.SF50_map;
    LPI_contra = OriSFdata_contra.LPI_map;
    ODI_contra = OriSFdata_contra.ODI;
    
    LHI_ipsi = OriSFdata_ipsi.LHI_map;
    CV_ipsi = OriSFdata_ipsi.CV_map;
    ori_hist_ipsi = OriSFdata_ipsi.OriTunPolarPlot;
    ori_bin_ipsi = OriSFdata_ipsi.angles;
    %sf_hist_ipsi = OriSFdata_ipsi.allsfcurve_dog;
    sf_tuning_ipsi = OriSFdata_ipsi.allsfcurve; % allsfcurve_dog;
    sfpref_ipsi = OriSFdata_ipsi.SF_map;
    sf_bin_deg_ipsi = OriSFdata_ipsi.sf_bin_deg;
    SF50_ipsi = OriSFdata_ipsi.SF50_map;
    LPI_ipsi = OriSFdata_ipsi.LPI_map;
    ODI_ipsi = OriSFdata_ipsi.ODI;
    %%
    figure
%     set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
      set(gcf,'position',[560   528   998   420])

    for kk = 1:5
        ii = electrode_position(kk,1);
        jj = electrode_position(kk,2);

        if  strcmp(eye_pref(kk),'c')     % contra
            plot_color = [0 0 0];
            LHI = LHI_contra(ii,jj);
            CV = CV_contra(ii,jj);
            ori_tuning_resp = ori_hist_contra{ii,jj};
            ori_bin = ori_bin_contra{ii,jj};
            %sf_hist = sf_hist_contra{ii,jj};
            sf_tuning_resp = sf_tuning_contra{ii,jj};
            sfpref = sfpref_contra(ii,jj);
            sf_bin_deg = sf_bin_deg_contra;
            SF50 = SF50_contra(ii,jj);
            LPI = LPI_contra(ii,jj);
            ODI = ODI_contra(ii,jj);
         elseif   strcmp(eye_pref(kk),'i')         % ipsi (orange)
            plot_color = [.9 .4 .16];   %[0.9290 0.6940 0.1250]
            LHI = LHI_ipsi(ii,jj);
            CV = CV_ipsi(ii,jj);
            ori_tuning_resp = ori_hist_ipsi{ii,jj};
            ori_bin = ori_bin_ipsi{ii,jj};
            %sf_hist = sf_hist_ipsi{ii,jj};
            sf_tuning_resp = sf_tuning_ipsi{ii,jj};
            
            sfpref = sfpref_ipsi(ii,jj);
            sf_bin_deg = sf_bin_deg_ipsi;
            SF50 = SF50_ipsi(ii,jj);
            LPI = LPI_ipsi(ii,jj);
            ODI = ODI_ipsi(ii,jj);
        
        end
        
        %% Smooth Curves
        sf_bin_deg_interpolated = linspace(sf_bin_deg(1),sf_bin_deg(end),n_interpolation);
        [sf_resp_dog_interpolate] = sf_tuning_gauss_fit(sf_tuning_resp,1,sf_lSamp,n_interpolation);
        
        % to make the range for the points on the sf curve in LOG scale 
        nsf = 10;
        ppd = 20; % pixels per degree of set up
        minsf = 20;
        maxsf = 350; % 600
        sf = logspace(log10(minsf), log10(maxsf), nsf);
        sf = 1./sf*ppd; %cycles per degree
        sf = fliplr(sf);

        ind_points = zeros(1,length(sf)); 
        for pp = 1 : length(sf) 
            [~,ind_points(pp)] = min(abs(sf(pp) - sf_bin_deg_interpolated)); 
        end 

        % to smooth the SF curve 
        sf_bin_deg_points = sf_bin_deg_interpolated(ind_points);
        sf_tuning_resp_points = sf_resp_dog_interpolate(ind_points); 
        
        % to smooth the ori curve 
        ori_tuning_resp_new = ori_tuning_resp;
        % ori_tuning_resp_new = ori_tuning_gauss_fit_fig3(ori_tuning_resp,0)'; 
        
        %%  
        % ori tuning
        subplot(3,7,kk+7)
        polar1(ori_bin,ori_tuning_resp_new',plot_color);
        %title(sprintf('cv:%.2f \n lhi:%.2f',CVRaw(electrodePos,ii),LHIRaw(electrodePos,ii)))

        % sf tuning
        subplot(3,7,kk+14)
        % plot(sf_bin_deg,sf_hist,'lineWidth',2,'Color',plot_color);
        plot(sf_bin_deg_interpolated,sf_resp_dog_interpolate,'lineWidth',2,'Color',plot_color);
        axis square
        set(gca, 'XScale', 'log');
        set(gca,'box','off')
        title(sprintf('ODI:%.2f \n LPI:%.2f \n LHI:%.2f',ODI,LPI,LHI))
        xlabel(sprintf('row : %.0f \n col : %.0f',ii,jj))
        hold on, scatter(sf_bin_deg_points,sf_tuning_resp_points,'Marker','o','MarkerEdgeColor',plot_color)
        ylim([ -0.1 1.1])
        xlim([0.05 1.5])
        set(gca,'XTick', [sf_bin_deg_interpolated(1) 1])
        set(gca,'xticklabel',num2str([0 1]'))
        
        subplot(3,7,[20.5 21]), hold on
        ylim([-0.1 1.2]), %xlim([logSF(1)-0.5 logSF(end)+0.5])
        xlim([0.05 1.5])
        
        set(gca, 'XScale', 'log');
        plot(sf_bin_deg_interpolated,sf_resp_dog_interpolate,'lineWidth',2,'Color',plot_color);
        ylabel('Responses', 'fontsize', 12);
        xlabel('SF (Cyc/Deg)', 'fontsize', 12);
        set(gca,'XTick', [sf_bin_deg_interpolated(1) 1])
        set(gca,'xticklabel',num2str([0 1]'))
        set(gca,'box','off','Tickdir','out','YTickLabel',[]);%,'xticklabel',num2str([ 1 ]'))

        %% LPI/ODI , LPI/LHI
        subplot(3,7,[6.5 7])
        hold on
        plot(LPI,ODI,'Marker','o','MarkerEdgeColor',plot_color,'MarkerFaceColor',plot_color);
        caxis([0 1])
        ylabel('|ODI|'), xlabel('LPI')
        ylim([0 1]), xlim([0 0.8]), box off
        dx = 0.00; dy = 0.07; % displacement so the text does not overlay the data points
        text(LPI+dx, ODI+dy, num2str(kk-1));
        
        subplot(3,7,[13.5 14])
        hold on
        plot(LPI,LHI,'Marker','o','MarkerEdgeColor',plot_color,'MarkerFaceColor',plot_color);
        caxis([0 1])
        ylabel('LHI'), xlabel('LPI')
        ylim([0 1]), xlim([0 0.8]), box off
        text(LPI+dx, LHI+dy, num2str(kk-1));
        
        plot_psths(ODI_contra(ii,jj), ODI_ipsi(ii,jj), kk, 1) 
    end

   
    
%% showing the maps 
     % debug = 0; 
    if debug == 1
        OriMapPref_contra = OriSFdata_contra.ori_map_smooth;
        OriMapPref_ipsi = OriSFdata_ipsi.ori_map_smooth;
        
%         figure; clf;
%         set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
%         subplot(2,3,1);imagesc(OriMapPref_contra);title('Orientation');colormap(gca,'hsv'); colorbar;axis square;
%         hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
%         subplot(2,3,2); imagesc(LHI_contra); title('LHI'); axis square; colormap(gca,'jet');colorbar
%         hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
%         subplot(2,3,3); imagesc(LPI_contra);axis square; title('LPI');colorbar
%         hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
%         annotation('textbox',[0.03 0.87 0.08 0.08],'String','Contra','EdgeColor','none','fontsize',16)
%         
%         subplot(2,3,4);imagesc(OriMapPref_ipsi);title('Orientation');colormap(gca,'hsv'); colorbar;axis square;
%         hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
%         subplot(2,3,5); imagesc(LHI_ipsi); title('LHI'); axis square; colormap(gca,'jet');colorbar
%         hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
%         subplot(2,3,6); imagesc(LPI_ipsi);axis square; title('LPI');colorbar
%         hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
%         annotation('textbox',[0.03 0.45 0.08 0.08],'String','Ipsi','EdgeColor','none','fontsize',16)
        
        x1 = (electrode_position(3,2)) - 10 ; 
        x2 = (electrode_position(3,2)) + 10 ; 
        y1 = (electrode_position(3,1)) - 10 ; 
        y2 = (electrode_position(3,1)) + 10 ; 
        
        subplot(3,7,1); imagesc(OriMapPref_contra);title('Orientation');colormap(gca,'hsv');axis square; axis off
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
        xlim([x1 x2]), ylim([y1 y2])
        subplot(3,7,2); imagesc(LHI_contra); title('LHI'); axis square; colormap(gca,'jet'); axis off
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
        xlim([x1 x2]), ylim([y1 y2])
        subplot(3,7,3); imagesc(LPI_contra);axis square; title('LPI'); axis off
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
        xlim([x1 x2]), ylim([y1 y2])
        colormap(gca,'jet')
        subplot(3,7,4); imagesc(ODCrtxPlt_smooth);axis square; title('OD'); axis off
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
        xlim([x1 x2]), ylim([y1 y2])
        colormap(gca,'gray')
        subplot(3,7,5); imagesc(SF50_contra);axis square; title('SF50'); axis off
        hold on, contour(ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
        xlim([x1 x2]), ylim([y1 y2])
        colormap(gca,'jet')
        
        
        for kkk = 1:5
            ii = electrode_position(kkk,1);
            jj = electrode_position(kkk,2);
            if strcmp(eye_pref(kkk),'c')
                subplot(3,7,1), hold on, plot(jj,ii,'ko', 'LineWidth', 2);
                subplot(3,7,2), hold on, plot(jj,ii,'ko', 'LineWidth', 2);
                subplot(3,7,3), hold on, plot(jj,ii,'ko', 'LineWidth', 2);
                subplot(3,7,4), hold on, plot(jj,ii,'ko', 'LineWidth', 2);
                subplot(3,7,5), hold on, plot(jj,ii,'ko', 'LineWidth', 2);
            elseif strcmp(eye_pref(kkk),'i')
                % [.9 .4 .16]  % orange color 
                subplot(3,7,1), hold on, plot(jj,ii,'o', 'LineWidth', 2, 'MarkerEdgeColor',[.9 .4 .16] );
                subplot(3,7,2), hold on, plot(jj,ii,'o', 'LineWidth', 2, 'MarkerEdgeColor',[.9 .4 .16] );
                subplot(3,7,3), hold on, plot(jj,ii,'o', 'LineWidth', 2, 'MarkerEdgeColor',[.9 .4 .16] );
                subplot(3,7,4), hold on, plot(jj,ii,'o', 'LineWidth', 2, 'MarkerEdgeColor',[.9 .4 .16] );
                subplot(3,7,5), hold on, plot(jj,ii,'o', 'LineWidth', 2, 'MarkerEdgeColor',[.9 .4 .16] );
            end
        end
    end

%     annotation('textbox',[0.03 0.82 0.08 0.08],'String','Contra Maps','EdgeColor','none','fontsize',16)
%     annotation('textbox',[0.03 0.87 0.08 0.08],'String','fig3c','EdgeColor','none','fontsize',16)
    
end


%% functions 
    %% fitting gaussian distribution to ori tuning curve 
    % it is working well at some locations 
    % it makes the broad band ori tuning curves sharp and results in no correlation btw CV(Orientation Selectivity) vs. ONOFF balance 
    
%     ori_tuning_resp_new = ori_tuning_gauss_fit_test1(ori_tuning_resp,0)'; 
%     ori_tuning_resp_new = ori_tuning_gauss_fit(ori_tuning_resp,0)'; 
%% 

function [sf_curve_dog] = sf_tuning_gauss_fit(sf_tuning_resp,gauss_fit,sf_lSamp,n_interpolation)

    sfallnorm = sf_tuning_resp - min(sf_tuning_resp(:)) ;
    sfallnorm = sfallnorm / max(sfallnorm(:)); 
    tmp = smooth(sf_lSamp,sfallnorm); 
    pp = spline(sf_lSamp, tmp);

%     nn = 20 ; % estimating the sf curve with 1000 points 
%     range_sf = 1:(1/nn):max(sf_lSamp); 
    range_sf = linspace(1,max(sf_lSamp),n_interpolation); 
    spltmp = ppval(pp,range_sf);% ppval(pp,1:.5:max(lSamp));
    spltmp = spltmp - min(spltmp);
    spltmp = spltmp/max(spltmp(:)); 
       
    if gauss_fit == 1 
        try % because this error happens : Inf computed by model function, fitting cannot continue.
            gfit = fit(range_sf',spltmp','gauss1') ;
            spldog = feval(gfit,range_sf);
        catch
            spldog = spltmp;
        end
    else 
        spldog = spltmp; 
    end 
        
    sf_curve_dog = spldog;%(1:nn:end); % gfit = fit(1:max(sf_lSamp),spltmp,'gauss2') ; 

end


function ori_tuning_resp_new = ori_tuning_gauss_fit_fig3(ori_tuning_resp,debug)
% fit gaussian to orientation tuning curve 
% the function keeps the curve with maximum response
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [pks,locs,w,~] = findpeaks(ori_tuning_resp); 
    [ind_peaks] = find(pks == 1 ); 
    ww = mean(w(ind_peaks)); 
    ww = ww / 1.5; % if it is not divided by a number, the two curves overlap on each other 
    
    [locs2] = find(ori_tuning_resp == 1);  % if a max is located at x == 1, findpeaks does not detect it 
    min_resp = min(ori_tuning_resp(:)); 

    xx = linspace(1,length(ori_tuning_resp),length(ori_tuning_resp));
    Gauss = zeros(length(locs2),length(ori_tuning_resp));
    Gauss_debug = zeros(length(locs2),length(ori_tuning_resp));
    for n = 1:length(locs2)
        % m = locs2(n); 
        %   Gauss(n,:) = pks(m) * exp(-((xx - locs(m))/w(m)).^2);
%         Gauss(n,:) = (1-min_resp) * exp(-((xx - locs2(n))/ww).^2) ;
%         Gauss_debug(n,:) = (1-min_resp) * exp(-((xx - locs2(n))/ww).^2) + min_resp;
        Gauss(n,:) = (1) * exp(-((xx - locs2(n))/ww).^2) ;
        Gauss_debug(n,:) = (1) * exp(-((xx - locs2(n))/ww).^2) + min_resp; 
    end
    % ori_tuning_resp_new = sum(Gauss,1) + min_resp; 
    ori_tuning_resp_new = sum(Gauss,1) ;
    
    %% debug 
    if debug == 1 
        figure 
        findpeaks(ori_tuning_resp)
        hold on 
        plot(xx,Gauss_debug,'--')

        figure, plot(ori_tuning_resp_new)
        %cv = circvar(angles_bin, ori_tuning_resp_new);
    end 

end 


function plot_psths(ODI_contra, ODI_ipsi, kk, show_fig)

    if show_fig
         tt = 0:0.01:1;
        y1 = normpdf(tt,.2,.03);  
        y2 = normpdf(tt,.55,.04); % the stimulus is presented 130 ms, the overl time range is 300 ms
        y3 = normpdf(tt,.37,.05);
        y = y1/max(y1(:)) + (y2/ max(y2(:)))*.8 + (y3/ max(y3(:)))*.5; % y2 is rebound
%         y = y + rand(1,length(tt)) * .1 ;
       
        contra_psth1 = y * ODI_contra;
        ipsi_psth1   = y * ODI_ipsi; 
        
        % filter = exp(3*xx); % passing the high amplitude response 
        contra_psth2 = exp(5 * contra_psth1); 
        contra_psth2 = contra_psth2 - min(contra_psth2(:)); 
        ipsi_psth2 = exp(5 * ipsi_psth1);  
        ipsi_psth2 = ipsi_psth2 - min(ipsi_psth2(:));
        
        max_all = max(max( contra_psth2(:)), max(ipsi_psth2(:))); 
        contra_psth_norm = contra_psth2 / max_all;
        ipsi_psth_norm = ipsi_psth2 / max_all;
        
        contra_psth_norm = contra_psth_norm + rand(1,length(tt)) * .1 ;
        ipsi_psth_norm = ipsi_psth_norm + rand(1,length(tt)) * .1 ;
        
        ss = subplot(3,7,kk);
        axis off 
        axes('Position',[ss.Position(1) ss.Position(2) ss.Position(3) ss.Position(4)*1.2]) % to match the amplitude in data
        hold on, plot(tt,contra_psth_norm, 'k')
        hold on, plot(tt,ipsi_psth_norm,'Color', [.9 .4 .16])
        box off, axis off 
    end 
    
end 

