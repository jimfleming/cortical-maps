function plot_receptive_field(OriSFdata,cortex_RF,electrodePos,reference_ori_map)
% LHIRaw = OriSFdata.LHI_map;
% CVRaw = OriSFdata.CV_map;
    
allCXrfONReference = cortex_RF.allCXrfON ;
allCXrfOFFReference = cortex_RF.allCXrfOFF ;

OriHistRFReference = OriSFdata.OriTunPolarPlot;
OriBinRFReference = OriSFdata.angles;

SFHistRFReference = OriSFdata.allsfcurve_dog;
sfpref = OriSFdata.SF_map; 
sf_bin_deg = OriSFdata.sf_bin_deg; 
SF50Raw = OriSFdata.SF50_map;

ori_pref_rf = OriSFdata.ori_pref_rf;

figure;clf
set(gcf,'position',[10         677        1673         300])
%set(fg11,'position',[10         677        1673         225])
%set(fg11,'position',[10         300        1673         450])
xs = linspace(0.02,0.99-(1/32),32);
width = 1/(32*1.3);
height = 1/(3*1.35);
caxisVal = 1;

for ii = 2:2:60
    RFrefON = allCXrfONReference{electrodePos,ii}; 
    RFrefOFF = allCXrfOFFReference{electrodePos,ii}; 
    RFreONOFF = RFrefON + RFrefOFF; 
    MaxValrf = max(max(abs(RFrefON(:))),max(abs(RFrefOFF(:)))); 
    
    axes('Position',[xs(ii/2+1),0.7,width,height])
    RFreONOFFnorm = RFreONOFF / MaxValrf; 
    imagesc(RFreONOFFnorm)
    %xlim([100 300]),ylim([100 300])
    title(num2str(ii))
    caxis([-caxisVal caxisVal])
    axis off
    axis square
    
    axes('Position',[xs(ii/2+1),0.55,width,height])
    RFrefONnorm = RFrefON / MaxValrf; 
    imagesc(RFrefONnorm)
    caxis([-caxisVal caxisVal])
    axis off
    axis square
    
    axes('Position',[xs(ii/2+1),0.4,width,height])
    RFrefOFFnorm = RFrefOFF / MaxValrf; 
    imagesc(RFrefOFFnorm)
        caxis([-caxisVal caxisVal])
    axis off
    axis square
    
    axes('Position',[xs(ii/2+1),0.2,width,height])
    ori_tuning = OriHistRFReference{electrodePos,ii};
    angles = OriBinRFReference{electrodePos,ii};%{electrodePos,ii};
    ph = polar1(angles,ori_tuning','k');
    %title(sprintf('cv:%.2f \n lhi:%.2f',CVRaw(electrodePos,ii),LHIRaw(electrodePos,ii)))
    title(sprintf('%.2f \n %.2f',ori_pref_rf(electrodePos,ii),reference_ori_map(electrodePos,ii)))
    
    axes('Position',[xs(ii/2+1),0.00,width,height])
    SFPlot = SFHistRFReference{electrodePos,ii};
    %SFPlot = (SFHistRFReference(electrodePos,ii,:));
    plot(SFPlot,'lineWidth',2);
    axis square
    set(gca, 'XScale', 'log');
    set(gca,'box','off','Tickdir','out','YTickLabel',[],'XTickLabel',[])
    title(sprintf('%.2f',SF50Raw(electrodePos,ii)))
    
    annotation('textbox',[0.01 0.88 0.1 0.1],'String',['row = ' num2str(electrodePos)],'EdgeColor','none','fontsize',14)
end

colormap('jet')
end 