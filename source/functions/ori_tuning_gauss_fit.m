function ori_tuning_resp_new = ori_tuning_gauss_fit(ori_tuning_resp,debug)
% fit gaussian to orientation tuning curve 
% the function keeps the curve with maximum response
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [pks,locs,w,~] = findpeaks(ori_tuning_resp); 
    [ind_peaks] = find(pks == 1 ); 
    ww = mean(w(ind_peaks)); 
    ww = ww / 2;
    
    [locs2] = find(ori_tuning_resp == 1);  % if a max is located at x == 1, findpeaks does not detect it 
    min_resp = min(ori_tuning_resp(:)); 

    xx = linspace(1,length(ori_tuning_resp),length(ori_tuning_resp));
    Gauss = zeros(length(locs2),length(ori_tuning_resp));
    Gauss_debug = zeros(length(locs2),length(ori_tuning_resp));
    for n = 1:length(locs2)
        % m = locs2(n); 
        %   Gauss(n,:) = pks(m) * exp(-((xx - locs(m))/w(m)).^2);
        Gauss(n,:) = (1-min_resp) * exp(-((xx - locs2(n))/ww).^2) ;
        Gauss_debug(n,:) = (1-min_resp) * exp(-((xx - locs2(n))/ww).^2) + min_resp; 
    end
    ori_tuning_resp_new = sum(Gauss,1) ;%+ min_resp; 

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

