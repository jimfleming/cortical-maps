function generate_fig3e_scatter_plot3(data_mature_stats,row_range,col_range,n_ori_smooth)
% fig3e scatter plot
% the onoff balance is measured based on the final max values of RFon and RFoff
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% row_cv_onoffbal = [100 200]; % rng 4
% range_cv_onoffbal = 20:220;  % rng 4
row_cv_onoffbal = row_range; 
range_cv_onoffbal = col_range;  
cvmap_contra_map = imresize(data_mature_stats.CV_map,n_ori_smooth);
% onoffbal_contra_map = imresize(data_contra_mature.CxOnOffBal,n_ori_smooth);

max_on = data_mature_stats.max_on_response_norm; 
max_off = data_mature_stats.max_off_response_norm; 
onoff_bal = abs(max_on + (-1) * max_off) ./ (max_on + max_off); 
onoff_bal = 1 - onoff_bal; 
if max(onoff_bal(:) > 1) 
    % onoff balance should be less than 1 always 
    error('There is a value more than 1')
end 
onoff_bal_interp = imresize( onoff_bal , n_ori_smooth );
onoff_bal_interp( onoff_bal_interp > 1 ) = 1;  % there are some cases, that interpolation causes the onoff balance to go beyond 1, although the uninterpolated values are under 1 

% in case there are more than one row is selected 
% cvmap_contra = []; 
% onoffbal_contra = []; 
% for qq = 1 : length(row_cv_onoffbal)
%     cvmap_contra = [cvmap_contra cvmap_contra_map(row_cv_onoffbal(qq),range_cv_onoffbal)]; 
%     onoffbal_contra = [onoffbal_contra onoff_bal_interp(row_cv_onoffbal(qq),range_cv_onoffbal)]; 
% end 

cvmap_datapoints = cvmap_contra_map(row_cv_onoffbal,range_cv_onoffbal);
onoffbal_datapoints = onoff_bal_interp(row_cv_onoffbal,range_cv_onoffbal);

figure, plot(cvmap_datapoints, onoffbal_datapoints,'ko','MarkerSize',20); lsline;
[r,p]=corr(cvmap_datapoints', onoffbal_datapoints');
title(['R = ', num2str(r), ' p = ', num2str(p)]);
axis square, ylabel('ONOFF balance'), xlabel('CV(OS)'),
set(gca,'box' ,'off','TickDir','OUT')
xlim([-0.05 1.1])
ylim([-0.05 1.1])

text = ['Contra (CV vs. ONOFF balance)' ', Num points : ' num2str(length(cvmap_datapoints)) ', row a : ' num2str(row_cv_onoffbal(1)) ...
    ', range : ' num2str(range_cv_onoffbal(1)) ':' num2str(range_cv_onoffbal(end ))];
annotation('textbox',[0.03 0.93 0.98 0.08],'String',text,'EdgeColor','none','fontsize',10)

