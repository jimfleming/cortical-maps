function generate_fig3e_temporal_filter(rf_mature,data_mature,electrode_position,beta_dist_val)
% fig3e
% figure 3 bottom panel (applying temporal filter to RFs)
% electrode_position = [25,32; 33,44;41,20; 36,38 ];
% beta_dist_val = [2,2,1,1.5];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fig3e 

OriHistRFReference = data_mature.OriTunPolarPlot;
OriBinRFReference = data_mature.angles;

figure(91), clf
% set(gcf,'position',[10         677        1673         300])
set(gcf,'position',[10         677        1200         300])

xs1 = linspace(0.01,.85,15);
xs2 = linspace(0.01,0.8,4);
caxisVal = 1;
width1 = .05;
height1 =.18 ;
width2 = .17 ;

height2 = .45;
sss = .06 ; 

for kk = 1 : 4
    ii = electrode_position(kk,1);
    jj = electrode_position(kk,2);
    
    rfON = rf_mature.allCXrfON{ii,jj};
    rfOFF = rf_mature.allCXrfOFF{ii,jj};
    RFSpaceSim = rfON + rfOFF;
    SingleRFNorm = RFSpaceSim / max(abs(RFSpaceSim(:)));
    
    %% adding temporal filter
    tt = 0:0.01:1;
    y1 = betapdf(tt,beta_dist_val(kk),4);
    y1 = y1 / max(y1(:));
%     figure(90),
%     plot(tt*80,y1)
%     set(gca,'box','off','tickdir','out')
%     xlabel('milisecond','fontsize',20)
%     ylabel('resp','fontsize',20)
    %%
    tt1 = linspace(0.05,.6,4);
    
    figure(92)
    for qq = 1 : 4
        [~,ind] = (min(abs(tt1(qq) - tt)));
        ww = y1(ind);
        %rf_temporal{qq} = SingleRFNorm * ww ;
        subplot(4,4,(qq-1) * 4 + kk)
        imagesc(SingleRFNorm * ww ),caxis([-caxisVal caxisVal]),axis square,
        axis off
    end
    colormap('jet')
    
    %%
    figure(91)

    col_move = [-1 0 1]; 
    for kkk = 1 : 3
        cc = col_move(kkk); 
        axes('Position',[xs2(kk) + sss * (kkk-1) ,0.7,width1,height1])
        ori_tuning = OriHistRFReference{ii,jj + cc};
        angles = OriBinRFReference{ii,jj};%{electrodePos,ii};
        ph = polar1(angles,ori_tuning','k');
    end
     
    axes('Position',[xs2(kk),0.2,width2,height2])
    imagesc(SingleRFNorm)
    caxis([-caxisVal caxisVal])
    axis off
    axis square
    colormap(gca,'jet')
%     set(gca,'box','off','Tickdir','out','YTickLabel',[],'xticklabel',[])
%     xlabel(sprintf('row : %.0f \n col : %.0f',ii,jj))
    
end

figure(91)
text1 = ['fig3e, Contra, points : [25,32; 33,44;41,20; 36,38 ]'];
annotation('textbox',[0.03 0.9 0.98 0.08],'String',text1,'EdgeColor','none','fontsize',10)

figure(92)
text2 = ['fig3e, Contra' ', points : [25,32; 33,44;41,20; 36,38 ],  betapdf(tt,beta dist val(kk),4), beta dist val = [2,2,1,1.5]' ];
annotation('textbox',[0.03 0.93 0.98 0.08],'String',text2,'EdgeColor','none','fontsize',10)


%% to reduce the size of the rf file 
% elect = [25,32; 33,44; 41,20; 36,38]; 
% 
% 
% rf_on_fig3e{60,60} = []; 
% rf_off_fig3e{60,60} = []; 
% 
% 
% for gg = 1 : 4 
%      ii = elect(gg,1); 
%      jj = elect(gg,2);
%      
%      rf_on_fig3e{ii,jj} = rf_contra_mature.allCXrfON{ii,jj};
%      rf_off_fig3e{ii,jj} = rf_contra_mature.allCXrfOFF{ii,jj};
%      
% end 
