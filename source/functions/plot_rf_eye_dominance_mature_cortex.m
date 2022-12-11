function plot_rf_eye_dominance_mature_cortex(OriSFdata,cortex_RF,ODCrtxPlt_smooth,electrodePos,eyepref)

    % ODI = OriSFdata.ODI; 
    LHIRaw = OriSFdata.LHI_map;
    CVRaw = OriSFdata.CV_map;
    LPI_map = OriSFdata.LPI_map;

    allCXrfONReference = cortex_RF.allCXrfON ;
    allCXrfOFFReference = cortex_RF.allCXrfOFF ;

    OriHistRFReference = OriSFdata.OriTunPolarPlot;
    OriBinRFReference = OriSFdata.angles;

    SFHistRFReference = OriSFdata.allsfcurve_dog;
    sfpref = OriSFdata.SF_map; 
    sf_bin_deg = OriSFdata.sf_bin_deg; 
    SF50Raw = OriSFdata.SF50_map;

    CxOnOffBal = OriSFdata.CxOnOffBal;

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

        axes('Position',[xs(ii/2+1),0.6,width,height])
        RFreONOFFnorm = RFreONOFF / MaxValrf; 
        imagesc(RFreONOFFnorm)
        %xlim([100 300]),ylim([100 300])

        caxis([-caxisVal caxisVal])
        axis off
        axis square
        if ODCrtxPlt_smooth(electrodePos,ii) < .5 % ipsi
            title(num2str(ii),'color',[0.9290 0.6940 0.1250]) % orange 
        else 
            title(num2str(ii),'color','k')
        end 

        axes('Position',[xs(ii/2+1),0.35,width,height])
        ori_tuning = OriHistRFReference{electrodePos,ii};
        angles = OriBinRFReference{electrodePos,ii};%{electrodePos,ii};
        ph = polar1(angles,ori_tuning','k');
    %     title(sprintf('cv:%.2f \n lhi:%.2f \n od:%.2f ',...
    %         CVRaw(electrodePos,ii),LHIRaw(electrodePos,ii),ODI(electrodePos,ii)))
        title(sprintf('cv:%.2f \n lhi:%.2f \n bal:%.2f ',...
            CVRaw(electrodePos,ii),LHIRaw(electrodePos,ii),CxOnOffBal(electrodePos,ii)))

        axes('Position',[xs(ii/2+1),0.1,width,height])
        SFPlot = SFHistRFReference{electrodePos,ii};
        %SFPlot = (SFHistRFReference(electrodePos,ii,:));
        plot(SFPlot,'lineWidth',2);
        axis square
        set(gca, 'XScale', 'log');
        set(gca,'box','off','Tickdir','out','YTickLabel',[],'XTickLabel',[])
        title(sprintf('sf:%.2f \n 50:%.2f\n lpi:%.2f',sfpref(electrodePos,ii), SF50Raw(electrodePos,ii),LPI_map(electrodePos,ii)))

        annotation('textbox',[0.01 0.88 0.1 0.1],'String',[eyepref ',row = ' num2str(electrodePos)],'EdgeColor','none','fontsize',14)
    end

    colormap('jet')
end 