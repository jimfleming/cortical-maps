function [ori_tuning_resp, sf_tuning_resp, sf_tuning_resp_dog, cv, ori_pref_deg, sf_pref_deg, sf50_deg, LPI, angles_bin, num_cycles_bin_deg] = ... 
                    fft_ori_tuning_uf(rf_input,sf_lSamp)

    %the size of RF should be ODD for fft calculation
    if mod(size(rf_input,1),2) == 0
        rf_input = rf_input(1:end-1,1:end-1);
    end

    %   lSamp = [1 2 3 4 5 6 7 8 10 12 14 16 18 20 24 28 32 36 42 50];   % spatial frequency
    img_fft = fftshift(fft2(rf_input));
    img_fft_abs = abs(img_fft);
    img_fft_abs = interp2(img_fft_abs,1);
    img_fft_abs_norm_sf = img_fft_abs/max(img_fft_abs(:));
    %img_fft_abs_norm_sf = img_fft_abs_norm_sf .^ .5 ; 
    img_fft_absLog = log(img_fft_abs + 1);
    
    img_fft_abs2 = img_fft_abs .^ 3; % to make the orientation tuning curve sharper
    img_fft_abs_norm_ori = img_fft_abs2/max(img_fft_abs2(:));
    %img_fft_abs_norm_sf2 = log(img_fft_abs2 + 1);
    
    center_fft = round(size(img_fft_abs2)/2);

    %% 
    [ori_tuning_resp, cv, ori_pref_deg, ori_pref_rad, angles_bin] = measure_ori_tuning_curve(img_fft_abs_norm_ori, center_fft, sf_lSamp); 

    sf_tuning_resp = 0; 
    sf_tuning_resp_dog = 0; 
    sf_pref_deg = 0; 
    sf50_deg = 0; 
    LPI = 0; 
    num_cycles_bin_deg = 0; 
end

function [ori_tuning_resp_new, cv, ori_pref_deg, ori_pref_rad, angles_bin] = measure_ori_tuning_curve(img_fft_abs_norm_ori,center_fft,sf_lSamp)

    nc = 2*16;
    rtest= 0:max(sf_lSamp(:)-1); % rtest= 1:max(sf_lSamp(:));
    lr = length(rtest);
    [xc,yc,~] = cylinder(rtest,nc);
    angles_bin = atan2(yc(end,:),xc(end,:));
    xr = xc+center_fft(1);
    yr = yc+center_fft(2);
    xr = round(xr);
    yr = round(yr);
    ori_tuning_resp = zeros([nc+1,1]);
    
    for i=1:length(ori_tuning_resp)
        for j = 1:lr
            if yr(j,i) && xr(j,i)
                ori_tuning_resp(i) = ori_tuning_resp(i)+(img_fft_abs_norm_ori(yr(j,i),xr(j,i)));
            end
        end
        ori_tuning_resp(i) = ori_tuning_resp(i)/lr;
    end
    % ori_tuning_resp(8:24) = .2; 
    
    % to make the tuning curve sharper 
    % ori_tuning_resp = smooth(1:length(ori_tuning_resp),ori_tuning_resp); 
    % ori_tuning_resp = ori_tuning_resp .^ 2 ; 
    
    angles_bin = -1*angles_bin; %
    angles_bin = mod(angles_bin, 2*pi);
    
    ori_tuning_resp = ori_tuning_resp / max(ori_tuning_resp);
    ori_tuning_resp_new = ori_tuning_resp; 
    
    %% fitting gaussian distribution to ori tuning curve 
    % it is working well at some locations 
    % it makes the broad band ori tuning curves sharp and results in no correlation btw CV(Orientation Selectivity) vs. ONOFF balance 
    
%     ori_tuning_resp_new = ori_tuning_gauss_fit_test1(ori_tuning_resp,0)'; 
%     ori_tuning_resp_new = ori_tuning_gauss_fit(ori_tuning_resp,0)'; 



%% 
%    cv = circvar(angles_bin, ori_tuning_resp_new');
     cv = 0; 
%     if cv > .8 
%         rand_ori = (rand(length(ori_tuning_resp),1) - .5) * .1; 
%         ori_tuning_resp = ori_tuning_resp + rand_ori; 
%         ori_tuning_resp = ori_tuning_resp / max(ori_tuning_resp);
%     end 
    [~, ind ] = max(ori_tuning_resp_new);
    ori_pref_rad = angles_bin(ind);
    ori_pref_deg = ori_pref_rad * (180/pi);

    if ori_pref_deg > 180
        ori_pref_deg = ori_pref_deg - 180;
    elseif ori_pref_deg< 0
        ori_pref_deg = ori_pref_deg + 180;
    end
end 

function sf_tuning_resp = measure_sf_tuning_curve(img_fft_abs_norm_sf,ori_pref_rad,center_fft,sf_lSamp)

    xSamp = sf_lSamp.*cos(ori_pref_rad)+ center_fft(1);
    ySamp = sf_lSamp.*sin(ori_pref_rad)+ center_fft(2);
    % y is increasing downward
    % ySamp = size(img_fft_abs_norm_sf,1) - ySamp;
    ySamp = 2*center_fft(2) - ySamp;

    xSamp = round(xSamp);
    ySamp = round(ySamp);
    sf_tuning_resp = zeros(1,length(sf_lSamp)); %% sf response container
    for is = 1:length(sf_lSamp)
        if ySamp(is) &&  xSamp(is)
            %sf_tuning_resp(is) = img_fft_absLog(ySamp(is), xSamp(is));% figure,plot(sfall)
            sf_tuning_resp(is) = img_fft_abs_norm_sf(ySamp(is), xSamp(is));% figure,plot(sf_tuning_resp)
        end
    end
end 

function [sf_pref_deg,sf50_deg,LPI,sf_curve_dog,num_cycles_vec_deg] = ...
    sf_measurement(sf_tuning_resp,rf_input_size,gauss_fit,sf_lSamp,pix2deg)
    % the average distance between cells is 100 Micron == 5 pixels 
    % 1 deg = 250 Microns in Cat retina, therefore 1 deg = 12.5 pixels , and the values in lSamp should be devided by 12.64 (158/12.5) to transfer the scale to cycle per degree 
    % cycle_deg_ratio = 12.64 ; 
    % sf_tuning_resp      :  fft power response along the major aixs of RF
    % rf_input_size       :  width of RF for sf measurement in deg 
    % gauss_fit           :  if equals 1, gauss2 is fit to the data (it runs slower)
    % sf_lSamp            :  sampling points in pixel from the center of fft 
    % pix2deg             :  12 pix2deg == 1 degree = 12 pixels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sfallnorm = sf_tuning_resp - min(sf_tuning_resp(:)) ;
    sfallnorm = sfallnorm / max(sfallnorm(:)); 
    tmp = smooth(sf_lSamp,sfallnorm); 
    pp = spline(sf_lSamp, tmp);

    nn = 2 ; 
    range_sf = 1:(1/nn):max(sf_lSamp); 
    spltmp = ppval(pp,range_sf);% ppval(pp,1:.5:max(lSamp));
    spltmp = spltmp - min(spltmp);
    spltmp = spltmp/max(spltmp(:)); 
       
    if gauss_fit == 1 
        try % because this error happens : Inf computed by model function, fitting cannot continue.
            gfit = fit(range_sf',spltmp','gauss1') ;
            spldog = feval(gfit,range_sf);
            
        catch
            spldog = spltmp;
        end
    else 
        spldog = spltmp; 
    end 
    
    [~, ix ]= max(spldog);
    freq_sf_max = range_sf(ix); 
    
    [~,in] = min(abs(spldog(ix:end) - .5));
    freq_sf50 = range_sf(in);
%     if ix == max(sf_lSamp(:))
%         in = 0;
%     end
    if ~isempty(in)
        %sf50 = (ix+in)/cycle_deg_ratio;
        CyclePix_sf50 = size(rf_input_size,1) / (freq_sf_max+freq_sf50); % Cycle in pixel
        sf50_deg = pix2deg / CyclePix_sf50 ;
    end
    LPI = spldog(1)/max(spldog);
    %   sf_pref = ix/cycle_deg_ratio;
    sf_curve_dog = spldog(1:nn:end); % gfit = fit(1:max(sf_lSamp),spltmp,'gauss2') ; 

    CyclePix = size(rf_input_size,1) / freq_sf_max; % Cycle in pixel
    sf_pref_deg = pix2deg / CyclePix;  % cycle per degree 
    
    num_cycles_vec_pix = size(rf_input_size,1) ./ (1:max(sf_lSamp)); % Cycle in pixel
    num_cycles_vec_deg = pix2deg ./ num_cycles_vec_pix ;
end

function sf_pref_deg2 = measure_sf_from_radius(img_fft_abs_norm_ori,center_fft,pix2deg)
    % find the peaks and select the one which is most distant from the center
    [y,x] = find(img_fft_abs_norm_ori == max(img_fft_abs_norm_ori(:)));
    point = [x(1), y(1)];
    RMaxFreq = norm(point-center_fft);
    CyclePix = size(img_fft_abs_norm_ori,1) / RMaxFreq; % Cycle in pixel
    sf_pref_deg2 = pix2deg / CyclePix;
end 
