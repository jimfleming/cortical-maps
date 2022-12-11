function [rr,pp] = show_all_correlation(OriSFdata,rect1,rect2,rand_points,rand_points_interpol,eyepref,show_fig)


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

cvcrop = imcrop(CV,rect1);
cv = cvcrop(rand_points);
sfcrop = imcrop(SF_map,rect1);
sf = sfcrop(rand_points);
lpicrop = imcrop(LPI_map,rect1);
lpi = lpicrop(rand_points);
sf50crop = imcrop(SF50_map,rect1);
sf50 = sf50crop(rand_points);
lhicrop = imcrop(LHI,rect1);
lhi = lhicrop(rand_points);

cvcrop2 = imcrop(CV2,rect2);
cv2 = cvcrop2(rand_points_interpol);
sfcrop2 = imcrop(SF2,rect2);
sf2 = sfcrop2(rand_points_interpol);
lpicrop2 = imcrop(LPI2,rect2);
lpi2 = lpicrop2(rand_points_interpol);
sf50crop2 = imcrop(SF50_2,rect2);
sf50_2 = sf50crop2(rand_points_interpol);
lhicrop2 = imcrop(LHI2,rect2);
lhi2 = lhicrop2(rand_points_interpol);

rr = zeros(12,1);
pp = zeros(12,1);

[rr(1),pp(1)]=corr(cv', sf');
[rr(2),pp(2)]=corr(cv', lpi');
[rr(3),pp(3)]=corr(cv', sf50');
[rr(4),pp(4)]=corr(lhi', sf');
[rr(5),pp(5)]=corr(lhi', lpi');
[rr(6),pp(6)]=corr(lhi', sf50');
[rr(7),pp(7)]=corr(cv2', sf2');
[rr(8),pp(8)]=corr(cv2', lpi2');
[rr(9),pp(9)]=corr(cv2', sf50_2');
[rr(10),pp(10)]=corr(lhi2', sf2');
[rr(11),pp(11)]=corr(lhi2', lpi2');
[rr(12),pp(12)]=corr(lhi2',sf50_2');

if show_fig == 1
    figure;
    set(gcf , 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
    
    subplot(2,6,1);
    plot(cv, sf,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f', rr(1), pp(1)));
    xlabel('CV'), ylabel('SF Pref');
    xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,2);
    plot(cv, lpi,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(2),pp(2)));
    xlim([0 1]);
    xlabel('CV'), ylabel('LPI');
    xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,3);
    plot(cv, sf50,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(3),pp(3)));
    xlim([0 1]); axis square;
    xlabel('CV'), ylabel('SF50');
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,4);
    plot(lhi, sf,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(4),pp(4)));
    xlim([0 1]); axis square;
    xlabel('LHI'), ylabel('SF Pref');
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,5);
    plot(lhi, lpi,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(5),pp(5)));
    xlabel('LHI'), ylabel('LPI');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,6);
    plot(lhi, sf50,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(6),pp(6)));
    xlabel('LHI'), ylabel('SF50');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    %   ****** Correlation for Interpolated version ****************
    subplot(2,6,7);
    plot(cv2, sf2,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(7),pp(7)));
    xlabel('CV'), ylabel('SF Pref');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,8);
    plot(cv2, lpi2,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(8),pp(8)));
    xlabel('CV'), ylabel('LPI');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,9);
    plot(cv2, sf50_2,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(9),pp(9)));
    xlabel('CV'), ylabel('SF50');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,10);
    plot(lhi2, sf2,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f', rr(10),pp(10)));
    xlabel('LHI'), ylabel('SF Pref');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,11);
    plot(lhi2, lpi2,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(11),pp(11)));
    xlabel('LHI'), ylabel('LPI');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    subplot(2,6,12);
    plot(lhi2, sf50_2,'ko'); lsline;
    title(sprintf('  R = %.4f \nP = %.4f',rr(12),pp(12)));
    xlabel('LHI'), ylabel('SF50');   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
    
    annotation('textbox',[0.03 0.87 0.08 0.08],'String',eyepref,'EdgeColor','none','fontsize',32)
end
