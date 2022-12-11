function [ori_result_deg,coverage_mod,coverage_old] = modify_ori_coverage2(primordial_ori_map,ONOFFODLabelSorted,NumPinwheel,CenterPinwheelSorted,acf_parameters_coverage)

% Convolve ori within each island to bring all the ori within each island 
% choose the one that has the pinwheel closest to center of an island (afferent cluster)

% Inputs : 
% NumPinwheel
% ONOFFODLabelSorted
% primordial_ori_map
% CenterPinwheelSorted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
niter = 10; 

acf_center_x = acf_parameters_coverage.acf_center_x;
acf_center_y = acf_parameters_coverage.acf_center_y;
acf_surround_x = acf_parameters_coverage.acf_surround_x; 
acf_surround_y = acf_parameters_coverage.acf_surround_y; 

ori_map_input =  primordial_ori_map;
ori_map_complex = degree_to_complexNum(ori_map_input); 
ori_map_real = real(ori_map_complex);
ori_map_imag = imag(ori_map_complex);

% wSize = size(primordial_ori_map,1) ; 
% kern = gen_diff_gauss(wSize,wSize, 0.717, 0.433, .5, 12.86, 25.72);   % from Swindale (1992)
% kern_phase = gen_diff_gauss_phase(wSize,wSize, 0.717, 0.717, .5, 12.86, 12.86,2.5);   % from Swindale (1992)
% kern = functionDoGfilter(15,15,30,30,1);  % (10,10,14,14,.9)
kern = functionDoGfilter(acf_center_x,acf_center_y,acf_surround_x,acf_surround_y,1);  % (10,10,14,14,.9)
% kern = functionDoGfilter(acf_center,acf_center,acf_surround,acf_surround,1);  % (10,10,14,14,.9)

kern_angle = linspace(0,135,4); 

%%
coverage_old = zeros(NumPinwheel,1);
coverage_mod = zeros(NumPinwheel,1);

ori_map_real_input = ori_map_real; 
ori_map_imag_input = ori_map_imag; 

n_kern_rotate = 3 ; 
for ii = 1 : NumPinwheel
    
    SelectedIsland = ONOFFODLabelSorted == ii; 
    coverage_old(ii) = calculate_coverage(ori_map_input,SelectedIsland);
    % SelectedOriIsland = ori_map_input .* SelectedIsland;
    [pw_locations] = get_pinwheel_location_and_charge_coverage(ori_map_input,0,0); 
    pw_dist_original = min(sqrt((pw_locations(:,1) - CenterPinwheelSorted{ii}(2)).^2 + (pw_locations(:,2) - CenterPinwheelSorted{ii}(1)).^2 ));
    
    clear ori_real_all ori_imag_all ori_deg_all
    ori_real_all{length(kern_angle)*n_kern_rotate} = []; 
    ori_imag_all{length(kern_angle)*n_kern_rotate} = [];
    ori_deg_all{length(kern_angle)*n_kern_rotate} = [];
    count = 0; 
    flag = 0; 
    for jj = 1:length(kern_angle)
        kern_angle_rotate = kern_angle(jj); 
        kern_real = imrotate(kern,kern_angle_rotate,'crop');
        for kk = 1 : n_kern_rotate
            count = count + 1 ; 
            kern_imag = imrotate(kern,kern_angle_rotate+45*kk,'crop');
            
            [ori_real_conv_new,~] = calculate_ori_conv(ori_map_real_input,kern_real,SelectedIsland,niter);
            [ori_imag_conv_new,~] = calculate_ori_conv(ori_map_imag_input,kern_imag,SelectedIsland,niter);
            ori_conv_new = ori_real_conv_new + ori_imag_conv_new * 1i ;
            ori_result_deg = complexNum_to_degree(ori_conv_new);
            % SelectedOriIsland2 = ori_result_deg .* SelectedIsland;
            coverage_new(count) = calculate_coverage(ori_result_deg,SelectedIsland);
            [pw_locations] = get_pinwheel_location_and_charge_coverage(ori_result_deg,0,0); 
            
            if ~isempty(pw_locations)
                pw_dist = sqrt((pw_locations(:,1) - CenterPinwheelSorted{ii}(2)).^2 + (pw_locations(:,2) - CenterPinwheelSorted{ii}(1)).^2 );
            else 
                continue 
            end 
            
            dist_pw_vec = min(pw_dist(:));
            if dist_pw_vec < pw_dist_original
                ind_max = count ; 
                pw_dist_original = dist_pw_vec; 
                
                flag = 1; 
            end
            
            if 0 == 1 % debug
                debug_conv(ori_map_real_input,ori_map_imag_input,ori_real_conv_new,ori_imag_conv_new,SelectedIsland,ori_map_input,ori_result_deg)
            end
            
            ori_real_all{count} = ori_real_conv_new;
            ori_imag_all{count} = ori_imag_conv_new;
            ori_deg_all{count} = ori_result_deg;
        end
    end

    if flag
        if 0 == 1 % debug
            debug_conv(ori_map_real_input,ori_map_imag_input,ori_real_all{ind_max},ori_imag_all{ind_max},SelectedIsland,ori_map_input,ori_deg_all{ind_max})
        end
        ori_map_real_input = ori_real_all{ind_max};
        ori_map_imag_input = ori_imag_all{ind_max};
        ori_map_input = ori_deg_all{ind_max}; 
        coverage_mod(ii) = coverage_new(ind_max); 
    else 
        coverage_mod(ii) =  coverage_old(ii); 
    end 
end 

ori_result = ori_map_real_input + ori_map_imag_input * 1i ;
ori_result_deg = complexNum_to_degree(ori_result);

end 


%% functions

function [ori_conv_island,mse_ori] = calculate_ori_conv(ori_map,kern,island,niter)
    ori_conv = ori_map;
    for ii = 1 : niter
        conv1 = conv2(ori_conv, kern, 'same');
        conv1 = conv1 / max(abs(conv1(:)));
        ori_conv = conv1 .* island + ori_map .* ~island;
    end
    ori_conv_island = ori_conv ;
    mse_ori = norm(ori_map(island(:)) - ori_conv_island(island(:))) / sum(island(:)); 
end 

function coverage = calculate_coverage(ori_map_deg,SelectedIsland)
    ori_deg_island = ori_map_deg(SelectedIsland(:));
    [uni_ori,~] = unique(floor(ori_deg_island/10));
    coverage = length(uni_ori) / 18 ; % 180 / 10 = 18 conditions 
end

function debug_conv(ori_map_real_input,ori_map_imag_input,ori_real_conv_new,ori_imag_conv_new,SelectedIsland,ori_map_input,ori_conv_new_deg)
    ori_imag_input_island = ori_map_imag_input .* SelectedIsland;
    ori_real_input_island = ori_map_real_input .* SelectedIsland;
    ori_imag_output_island = ori_imag_conv_new .* SelectedIsland;
    ori_real_output_island = ori_real_conv_new .* SelectedIsland;

    figure
    subplot(232),imagesc(ori_imag_input_island), colormap gray, axis('square'), colorbar, title('Real input'),
    subplot(233),imagesc(ori_real_input_island), colormap gray, axis('square'), colorbar, title('Imag input'),
    subplot(235),imagesc(ori_imag_output_island), colormap gray, axis('square'), colorbar, title('Real output'),
    subplot(236),imagesc(ori_real_output_island), colormap gray, axis('square'), colorbar, title('Imag output'),

    subplot(231),imagesc(ori_map_input), colormap(gca,'hsv'), axis('square'),colorbar, title('Ori input')
    hold on, hold on, contour(SelectedIsland, 1, 'k', 'LineWidth', 5);
    subplot(234),imagesc(ori_conv_new_deg), colormap(gca,'hsv'), axis('square'),colorbar, title('Ori output')
    hold on, hold on, contour(SelectedIsland, 1, 'k', 'LineWidth', 5);
end 
