function debug_RF(rf_input,angles_bin,ori_tuning_resp,ori_pref,sf_tuning_resp,sf_lSamp,sf_pref_deg,pix2deg)
    %debug
    figure;
    set(gcf,'Units','Normalized','Outerposition',[0 0 .8 .8])
    subplot(131),imagesc(rf_input),colormap jet;colorbar,axis square
    
    subplot(132),polar1(angles_bin,ori_tuning_resp','k');hold on;axis square
    title(sprintf('Ori : %.2f ',ori_pref))
    
    sfallnorm = sf_tuning_resp - min(sf_tuning_resp(:)) ;
    sfallnorm = sfallnorm / max(sfallnorm(:));
    YCyclePix = size(rf_input,1) ./ (1:max(sf_lSamp)); % Cycle in pixel
    YCyclePerDeg = YCyclePix ./ pix2deg;
    YCyclePerDeg = 1 ./YCyclePerDeg;
    subplot(133),plot(YCyclePerDeg,sfallnorm,'lineWidth',2);axis square
    title(sprintf('SF : %.2f CPD',sf_pref_deg))
    xlabel('Cyc/Deg')
    set(gca,'box','off','Tickdir','out');
end
