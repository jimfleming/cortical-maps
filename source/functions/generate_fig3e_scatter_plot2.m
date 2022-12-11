function generate_fig3e_scatter_plot(data_contra_mature,row_range,col_range,n_ori_smooth)
% fig3 e scatter plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% row_cv_onoffbal = [100 200]; % rng 4
% range_cv_onoffbal = 20:220;  % rng 4
row_cv_onoffbal = row_range; 
range_cv_onoffbal = col_range;  
cvmap_contra_map = imresize(data_contra_mature.CV_map,n_ori_smooth);
onoffbal_contra_map = imresize(data_contra_mature.CxOnOffBal,n_ori_smooth);

cvmap_contra = []; 
onoffbal_contra = []; 

for qq = 1 : length(row_cv_onoffbal)
    cvmap_contra = [cvmap_contra cvmap_contra_map(row_cv_onoffbal(qq),range_cv_onoffbal)]; 
    onoffbal_contra = [onoffbal_contra onoffbal_contra_map(row_cv_onoffbal(qq),range_cv_onoffbal)]; 
end 

onoffbal_contra = 1 - onoffbal_contra ; 
% figure, plot(onoffbal_contra, cvmap_contra,'ko'); lsline;
% [r,p]=corr(onoffbal_contra', cvmap_contra');
figure, plot(cvmap_contra, onoffbal_contra,'ko','MarkerSize',20); lsline;
[r,p]=corr(cvmap_contra', onoffbal_contra');
title(['R = ', num2str(r), ' p = ', num2str(p)]);
axis square, ylabel('ONOFF balance'), xlabel('CV(OS)'),
set(gca,'box' ,'off','TickDir','OUT')
xlim([-0.05 1.1])
ylim([-0.05 1.1])

text = ['Contra (CV vs. ONOFF balance)' ', Num points : ' num2str(length(cvmap_contra)) ', row a : ' num2str(row_cv_onoffbal(1)) ...
    ', range : ' num2str(range_cv_onoffbal(1)) ':' num2str(range_cv_onoffbal(end ))];
annotation('textbox',[0.03 0.93 0.98 0.08],'String',text,'EdgeColor','none','fontsize',10)

