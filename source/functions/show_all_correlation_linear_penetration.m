function [rr,pp] = show_all_correlation_linear_penetration(OriSFdata,electrode_row,electrode_row_inperp,eyepref,show_fig)


CV = OriSFdata.CV_map;
SF_map = OriSFdata.SF_map;
LPI_map = OriSFdata.LPI_map;
SF50_map = OriSFdata.SF50_map;
LHI = OriSFdata.LHI_map;

CV2 = OriSFdata.cv_map_interpolated;
SF2 = OriSFdata.sf_map_interpolated;
LPI2 = OriSFdata.LPI_intepolated;
SF50_2 = OriSFdata.sf50_map_interpolated;
LHI2 = OriSFdata.LHI_interpolated;

%%
alpha = .1; % do not consider 10% of the maps from each side 
bound1 = round(size(CV,1) * alpha); 
range1 = bound1 : (size(CV,1)-bound1); 

bound2 = round(size(CV2,1) * alpha); 
range2 = bound2 : (size(CV2,1)-bound2); 

cv_sample = CV(electrode_row,range1);
sf_sample = SF_map(electrode_row,range1);
lpi_sample = LPI_map(electrode_row,range1);
sf50_sample = SF50_map(electrode_row,range1);
lhi_sample = LHI(electrode_row,range1);

cv2_sample = CV2(electrode_row_inperp,range2);
sf2_sample = SF2(electrode_row_inperp,range2);
lpi2_sample = LPI2(electrode_row_inperp,range2);
sf50_2_sample = SF50_2(electrode_row_inperp,range2);
lhi2_sample = LHI2(electrode_row_inperp,range2);

rr = zeros(12,1);
pp = zeros(12,1);

[rr(1),pp(1)]=corr(cv_sample', sf_sample');
[rr(2),pp(2)]=corr(cv_sample', lpi_sample');
[rr(3),pp(3)]=corr(cv_sample', sf50_sample');
[rr(4),pp(4)]=corr(lhi_sample', sf_sample');
[rr(5),pp(5)]=corr(lhi_sample', lpi_sample');
[rr(6),pp(6)]=corr(lhi_sample', sf50_sample');
[rr(7),pp(7)]=corr(cv2_sample', sf2_sample');
[rr(8),pp(8)]=corr(cv2_sample', lpi2_sample');
[rr(9),pp(9)]=corr(cv2_sample', sf50_2_sample');
[rr(10),pp(10)]=corr(lhi2_sample', sf2_sample');
[rr(11),pp(11)]=corr(lhi2_sample', lpi2_sample');
[rr(12),pp(12)]=corr(lhi2_sample',sf50_2_sample');

if show_fig == 1
    figure;
    set(gcf , 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
    
    subplot(2,6,1);
    plot(cv_sample, sf_sample,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f', rr(1), pp(1)));
    xlabel('CV'), ylabel('SF Pref');
    xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,2);
    plot(cv_sample, lpi_sample,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(2),pp(2)));
    xlim([0 1]);
    xlabel('CV'), ylabel('LPI');
    xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,3);
    plot(cv_sample, sf50_sample,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(3),pp(3)));
    xlim([0 1]); axis square;
    xlabel('CV'), ylabel('SF50');
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,4);
    plot(lhi_sample, sf_sample,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(4),pp(4)));
    xlim([0 1]); axis square;
    xlabel('LHI'), ylabel('SF Pref');
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,5);
    plot(lhi_sample, lpi_sample,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(5),pp(5)));
    xlabel('LHI'), ylabel('LPI');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,6);
    plot(lhi_sample, sf50_sample,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(6),pp(6)));
    xlabel('LHI'), ylabel('SF50');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    %   ****** Correlation for Interpolated version ****************
    subplot(2,6,7);
    plot(cv2_sample, sf2_sample,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(7),pp(7)));
    xlabel('CV'), ylabel('SF Pref');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,8);
    plot(cv2_sample, lpi2_sample,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(8),pp(8)));
    xlabel('CV'), ylabel('LPI');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,9);
    plot(cv2_sample, sf50_2_sample,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(9),pp(9)));
    xlabel('CV'), ylabel('SF50');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,10);
    plot(lhi2_sample, sf2_sample,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f', rr(10),pp(10)));
    xlabel('LHI'), ylabel('SF Pref');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,11);
    plot(lhi2_sample, lpi2_sample,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(11),pp(11)));
    xlabel('LHI'), ylabel('LPI');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,12);
    plot(lhi2_sample, sf50_2_sample,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(12),pp(12)));
    xlabel('LHI'), ylabel('SF50');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    annotation('textbox',[0.03 0.87 0.08 0.08],'String',['Linear(' eyepref ')'],'EdgeColor','none','fontsize',32)
end
