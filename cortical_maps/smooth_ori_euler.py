# function [ori_map_smooth,ori_map_interpolated] = smooth_ori_euler(developed_reference_ori_map,ODCrtxPlt_smooth,n_ori_smooth,debug,show_fig)
#     %   This function plots
#     %   INPUT   :   OriSpaceAll
#     %   OUTPUT  :

#     %%
#     OriSpaceAll = developed_reference_ori_map;
#     OriMapTemp = zeros(size(OriSpaceAll));
#     Nangles = linspace(22.5, 180, 8);
#     % To measure fft , the norm of Ori map are calculated based on being 1 for
#     % the desired point and gradually becomes zero for a location selective to
#     % stilmuli 90 degrees rotated. This interval can be reduced to make the periodicity
#     % (fft ring shape) more pronounced

#     deg2rad = (pi/180);
#     rad2deg = 1/deg2rad;
#     angleRange = 1:180;
#     clear OriMapNorm
#     OriMapNorm{length(Nangles)} = [];
#     for ii = 1:length(Nangles)
#         phaseDeg = Nangles(ii);
#         phaseDeg = 180 - phaseDeg;
#         probAngleCos = cos(2 * (angleRange+ phaseDeg ) .* deg2rad ) ;  %   figure,plot(angleRange,probAngleCos)
#         probAngleCosNorm = (probAngleCos + 1 ) /2;%   figure,plot(angleRange,probAngleCosNorm)

#         for jj = 1:size(OriSpaceAll,1)
#             for kk = 1:size(OriSpaceAll,2)
#                 OriPixVal = round(OriSpaceAll(jj,kk)); % figure,imagesc(OriSpaceAll), colormap('hsv')
#                 if OriPixVal <= 0
#                     OriPixVal = OriPixVal + 180;
#                 elseif OriPixVal > 180
#                     OriPixVal = OriPixVal - 180;
#                 end

#                 OriMapTemp(jj,kk) = probAngleCosNorm(OriPixVal) ;
#             end
#         end
#         OriMapNorm{ii} = OriMapTemp;
#     end

#     %%
#     test22_202 =  OriMapNorm{1};
#     test45_225 = OriMapNorm{2};
#     test67_247 = OriMapNorm{3};
#     test90_270 = OriMapNorm{4};
#     test112_292 =  OriMapNorm{5};
#     test135_315 = OriMapNorm{6};
#     test157_337 = OriMapNorm{7};
#     test0_180 = OriMapNorm{8};

#     fsig = 1.2;
#     myfilter = fspecial('gaussian', [round(fsig)*6+1 round(fsig)*6+1], fsig);

#     test0_180 = filter2(myfilter,test0_180,  'same');
#     test45_225 = filter2(myfilter,test45_225,  'same');
#     test90_270 = filter2(myfilter,test90_270,  'same');
#     test135_315 = filter2(myfilter,test135_315,  'same');
#     test22_202 = filter2(myfilter,test22_202,  'same');
#     test67_247 = filter2(myfilter,test67_247,  'same');
#     test112_292 = filter2(myfilter,test112_292, 'same');
#     test157_337 = filter2(myfilter,test157_337,  'same');

#     X0 = exp(2*1i*360*(pi/180)).*test0_180; %   figure,imagesc()
#     X22 = exp(2*1i*22.5*(pi/180)).*test22_202;
#     X45 = exp(2*1i*45*(pi/180)).*test45_225;
#     X67 = exp(2*1i*67.5*(pi/180)).*test67_247;
#     X90 = exp(2*1i*90*(pi/180)).*test90_270;
#     X112 = exp(2*1i*112.5*(pi/180)).*test112_292;
#     X135 = exp(2*1i*135*(pi/180)).*test135_315;
#     X157 = exp(2*1i*157.5*(pi/180)).*test157_337;
#     %sumX = (X0 + X45 + X90 + X135)/4;
#     sumX = (X0 +X22+ X45 + X67+ X90 + X112 + X135 + X157)/8;
#     ori_map_smooth = angle(sumX)/2*(rad2deg) ; %
#     ori_map_smooth(ori_map_smooth<=0) = ori_map_smooth(ori_map_smooth<=0) + 180;

#     %% Smoothed Version
#     ResizeNum   = n_ori_smooth ;
#     test0_180S  = imresize(test0_180,ResizeNum);
#     test22_202S = imresize(test22_202,ResizeNum);
#     test45_225S = imresize(test45_225,ResizeNum);
#     test67_247S = imresize(test67_247,ResizeNum);
#     test90_270S = imresize(test90_270,ResizeNum);
#     test112_292S = imresize(test112_292,ResizeNum);
#     test135_315S = imresize(test135_315,ResizeNum);
#     test157_337S = imresize(test157_337,ResizeNum);

#     X0S  = exp(2*1i*360*(pi/180)).*test0_180S;
#     X22S = exp(2*1i*22.5*(pi/180)).*test22_202S;
#     X45S = exp(2*1i*45*(pi/180)).*test45_225S;
#     X67S = exp(2*1i*67.5*(pi/180)).*test67_247S;
#     X90S = exp(2*1i*90*(pi/180)).*test90_270S;
#     X112S = exp(2*1i*112.5*(pi/180)).*test112_292S;
#     X135S = exp(2*1i*135*(pi/180)).*test135_315S;
#     X157S = exp(2*1i*157.5*(pi/180)).*test157_337S;
#     sumXS = (X0S +X22S+ X45S + X67S+ X90S + X112S + X135S + X157S)/8;
#     ori_map_interpolated = angle(sumXS)/2*(rad2deg);
#     ori_map_interpolated(ori_map_interpolated<=0) = ori_map_interpolated(ori_map_interpolated<=0) + 180;
# end

# function result = fft_measure(input)
#     % Measure fft of input image
#     % input     :   input orientation map component
#     % result    :   modified fft for presentation
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     fft_input = fft2(input);
#     fft_input(1) = 0; % assign dc componenet to zero
#     fft_log = log(1+abs(fftshift(fft_input)));
#     fft_log_sharp = fft_log .^ 3;
#     fft_log_intrep = interp2(fft_log_sharp, 2);
#     result = fft_log_intrep;
# end
