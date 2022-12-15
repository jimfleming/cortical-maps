def orientation_map_development2():
    pass


# function [developed_ori_map,OriReal,OriImag,Selectivity] = orientation_map_development2(primordial_ori_map, acf_parameters_smooth, num_iter, show_fig)
#     % Modified from Swindale (1992)
#     % 1 pixel = 50 microns
#     %
#     % orientation_map: the orientation map in degree
#     % num_iter: the number of iterations for developing the orientation_map
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#     acf_center_x = acf_parameters_smooth.acf_center_x;
#     acf_center_y = acf_parameters_smooth.acf_center_y;
#     acf_surround_x = acf_parameters_smooth.acf_surround_x;
#     acf_surround_y = acf_parameters_smooth.acf_surround_y;

#     z = degree_to_complexNum(primordial_ori_map);
#     z_noise = z .* normrnd(0, 0.03, size(z));
#     % to reduce the orientation selectivity (initially is 1, which is the maximum) to 20% and add noise without changing the orientations
#     z = z / 5 + z_noise;
#     % kern_z = gen_diff_gauss(size(orientation_map, 2), size(orientation_map, 1), 0.717, 0.433, 1.0, 12.86, 25.72);   % from Swindale (1992)
#     % based on Swindlae paper, Ad1 < Bd2 , if it is transforred to ASF termination, Wcenter < Wsurround
#     % d1 = sigma_center^2, d2 = sigma_surround^2
#     % kern_z = functionDoGfilter(15,15,30,30,1);  % (10,10,14,14,.9)
#     kern_z = functionDoGfilter(acf_center_x,acf_center_y,acf_surround_x,acf_surround_y,1);  % (10,10,14,14,.9)

#     for ii = 1:num_iter
#         ori_selectivity = (real(z).^2 + imag(z).^2) .^ 0.5;
#         z_conv = conv2(z, kern_z, 'same');
#         z_conv_abs = (real(z_conv).^2 + imag(z_conv).^2) .^ 0.5;
#         z_conv = z_conv / max(max(z_conv_abs));
#         change = z_conv .* (1 - ori_selectivity);
#         z = z + change;
#     end
#     developed_ori_map = complexNum_to_degree(z);
#     OriReal = real(z);
#     OriImag = imag(z);
#     Selectivity = (OriReal.^2 + OriImag.^2) .^ 0.5;

#     if show_fig == 1
#         figure; subplot(1,2,1); imagesc(primordial_ori_map); title('Primordial Orientation map','FontSize',20); axis square; colormap hsv; colorbar;
#         subplot(1,2,2); imagesc(developed_ori_map); title('Mature Orientation map','FontSize',20); axis square; colormap hsv; colorbar;
#     end
# end
