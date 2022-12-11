function plot_rf_linear_lhi(OriSFdata,cortex_RF,ODCrtxPlt_smooth,electrodePos,eyepref)

LHI_map_2d = OriSFdata.LHI_map;
CVRaw = OriSFdata.CV_map;
ori_pref_rf = OriSFdata.ori_pref_rf;

allCXrfONReference = cortex_RF.allCXrfON ;
allCXrfOFFReference = cortex_RF.allCXrfOFF ;

OriHistRFReference = OriSFdata.OriTunPolarPlot;
OriBinRFReference = OriSFdata.angles;

%% LHI linear measurement 
channel_range = 2:2:60; 
oripreference = ori_pref_rf(electrodePos,channel_range); 
[LHI_linear] = calculate_linear_LHI(length(channel_range), oripreference,180); 
lhi_vec_linear = nan(length(channel_range),1);
lhi_vec_linear(channel_range) = LHI_linear;


show_lhi_cv_corr(CVRaw,LHI_map_2d,lhi_vec_linear,electrodePos,10:2:50)

%%
figure;clf
set(gcf,'position',[10         677        1673         300])
%set(fg11,'position',[10         677        1673         225])
%set(fg11,'position',[10         300        1673         450])
xs = linspace(0.02,0.99-(1/32),32);
width = 1/(32*1.3);
height = 1/(3*1.35);
caxisVal = 1;

for ii = channel_range
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

    axes('Position',[xs(ii/2+1),0.2,width,height])
    ori_tuning = OriHistRFReference{electrodePos,ii};
    angles = OriBinRFReference{electrodePos,ii};%{electrodePos,ii};
    ph = polar1(angles,ori_tuning','k');
    %LHI2  = LHI_linear(round(ii/2)); 
    LHI2 = lhi_vec_linear(ii);
    title(sprintf('cv:%.2f \n lhi:%.2f \n %.2f \n ori:%.0f',CVRaw(electrodePos,ii),LHI_map_2d(electrodePos,ii),LHI2,ori_pref_rf(electrodePos,ii)))
    
    annotation('textbox',[0.01 0.88 0.1 0.1],'String',[eyepref ',row = ' num2str(electrodePos)],'EdgeColor','none','fontsize',14)
end

colormap('jet')

end 
 

function [lhi] = calculate_linear_LHI(num_channels, oripreference,sigma)
    %sigma=180;
    
    lhi = nan(1,num_channels);
    channels = 1:num_channels; 
    c1=1; 
    for ich=4:num_channels-3 
        distance = abs(channels-channels(ich))*100;    
        distanceweighting = exp(-(distance.^2)./(2*sigma^2))*c1;
        multiplier = 1/sum(distanceweighting);
        distanceweighting = distanceweighting * multiplier;
        Magnitude = abs( sum( distanceweighting .* exp(2*1i*(oripreference)) ) );
        lhi(ich) = Magnitude;
    end
end 


function show_lhi_cv_corr(CVRaw,LHI_map_2d,lhi_vec_linear,electrodePos,corr_range)

    cv = CVRaw(electrodePos,corr_range);
    lhi = LHI_map_2d(electrodePos,corr_range);
    lhi2 = lhi_vec_linear(corr_range);

    figure; clf;
    subplot(1,2,1);
    plot(lhi, cv,'ko'); lsline;
    [r1,p1]=corr(lhi', cv');
    title(['LHI map  R = ', num2str(r1), ' p = ', num2str(p1)]);
    xlabel('LHI'), ylabel('CV');   ylim([0 1]);   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')

    subplot(1,2,2);
    plot(lhi2, cv,'ko'); lsline;
    [r2,p2]=corr(lhi2, cv');
    title(['LHI linear  R = ', num2str(r2), ' p = ', num2str(p2)]);
    xlabel('LHI'), ylabel('CV');   ylim([0 1]);   xlim([0 1]); axis square;
    set(gca,'Tickdir','OUT','Box','OFF')
end 

