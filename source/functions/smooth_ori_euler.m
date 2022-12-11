function [ori_map_smooth,ori_map_interpolated] = smooth_ori_euler(developed_reference_ori_map,ODCrtxPlt_smooth,n_ori_smooth,debug,show_fig)
%   This function plots 
%
%   INPUT   :   OriSpaceAll
%
%   OUTPUT  :   

%% 
OriSpaceAll = developed_reference_ori_map; 
OriMapTemp = zeros(size(OriSpaceAll)); 
Nangles = linspace(22.5, 180, 8);
% To measure fft , the norm of Ori map are calculated based on being 1 for
% the desired point and gradually becomes zero for a location selective to
% stilmuli 90 degrees rotated. This interval can be reduced to make the periodicity
% (fft ring shape) more pronounced 

deg2rad = (pi/180);
rad2deg = 1/deg2rad; 
angleRange = 1:180; 
clear OriMapNorm
OriMapNorm{length(Nangles)} = [];
for ii = 1:length(Nangles)
    phaseDeg = Nangles(ii);
    phaseDeg = 180 - phaseDeg;
    probAngleCos = cos(2 * (angleRange+ phaseDeg ) .* deg2rad ) ;  %   figure,plot(angleRange,probAngleCos)
    probAngleCosNorm = (probAngleCos + 1 ) /2;%   figure,plot(angleRange,probAngleCosNorm)
    
    for jj = 1:size(OriSpaceAll,1)
        for kk = 1:size(OriSpaceAll,2)
            OriPixVal = round(OriSpaceAll(jj,kk)); % figure,imagesc(OriSpaceAll), colormap('hsv')
            if OriPixVal <= 0 
                OriPixVal = OriPixVal + 180; 
            elseif OriPixVal > 180
                OriPixVal = OriPixVal - 180; 
            end 
            
            OriMapTemp(jj,kk) = probAngleCosNorm(OriPixVal) ; 
        end 
    end 
    OriMapNorm{ii} = OriMapTemp; 
end

%%
test22_202 =  OriMapNorm{1};
test45_225 = OriMapNorm{2};
test67_247 = OriMapNorm{3};
test90_270 = OriMapNorm{4};
test112_292 =  OriMapNorm{5};
test135_315 = OriMapNorm{6};
test157_337 = OriMapNorm{7};
test0_180 = OriMapNorm{8};

fsig = 1.2;
myfilter = fspecial('gaussian', [round(fsig)*6+1 round(fsig)*6+1], fsig);

test0_180 = filter2(myfilter,test0_180,  'same');
test45_225 = filter2(myfilter,test45_225,  'same');
test90_270 = filter2(myfilter,test90_270,  'same');
test135_315 = filter2(myfilter,test135_315,  'same');
test22_202 = filter2(myfilter,test22_202,  'same');
test67_247 = filter2(myfilter,test67_247,  'same');
test112_292 = filter2(myfilter,test112_292, 'same');
test157_337 = filter2(myfilter,test157_337,  'same');

X0 = exp(2*1i*360*(pi/180)).*test0_180; %   figure,imagesc()
X22 = exp(2*1i*22.5*(pi/180)).*test22_202;
X45 = exp(2*1i*45*(pi/180)).*test45_225;
X67 = exp(2*1i*67.5*(pi/180)).*test67_247;
X90 = exp(2*1i*90*(pi/180)).*test90_270;
X112 = exp(2*1i*112.5*(pi/180)).*test112_292;
X135 = exp(2*1i*135*(pi/180)).*test135_315;
X157 = exp(2*1i*157.5*(pi/180)).*test157_337;
%sumX = (X0 + X45 + X90 + X135)/4;
sumX = (X0 +X22+ X45 + X67+ X90 + X112 + X135 + X157)/8;
ori_map_smooth = angle(sumX)/2*(rad2deg) ; %  
ori_map_smooth(ori_map_smooth<=0) = ori_map_smooth(ori_map_smooth<=0) + 180; 

%% Smoothed Version 
ResizeNum   = n_ori_smooth ; 
test0_180S  = imresize(test0_180,ResizeNum);
test22_202S = imresize(test22_202,ResizeNum);
test45_225S = imresize(test45_225,ResizeNum);
test67_247S = imresize(test67_247,ResizeNum);
test90_270S = imresize(test90_270,ResizeNum);
test112_292S = imresize(test112_292,ResizeNum);
test135_315S = imresize(test135_315,ResizeNum);
test157_337S = imresize(test157_337,ResizeNum);

X0S  = exp(2*1i*360*(pi/180)).*test0_180S;
X22S = exp(2*1i*22.5*(pi/180)).*test22_202S;
X45S = exp(2*1i*45*(pi/180)).*test45_225S;
X67S = exp(2*1i*67.5*(pi/180)).*test67_247S;
X90S = exp(2*1i*90*(pi/180)).*test90_270S;
X112S = exp(2*1i*112.5*(pi/180)).*test112_292S;
X135S = exp(2*1i*135*(pi/180)).*test135_315S;
X157S = exp(2*1i*157.5*(pi/180)).*test157_337S;
sumXS = (X0S +X22S+ X45S + X67S+ X90S + X112S + X135S + X157S)/8;
ori_map_interpolated = angle(sumXS)/2*(rad2deg);
ori_map_interpolated(ori_map_interpolated<=0) = ori_map_interpolated(ori_map_interpolated<=0) + 180; 

%%
if show_fig == 1
    figure
    subplot(4,4,1); imagesc(test0_180, [0 1]); title('0 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(4,4,2); imagesc(test45_225, [0 1]); title('45 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(4,4,5); imagesc(test90_270, [0 1]); title('90 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(4,4,6); imagesc(test135_315, [0 1]); title('135 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(2,2,2);
    imagesc(ori_map_smooth); colormap(gca,'hsv'); colorbar ;axis square;
    axis off
    
    subplot(4,4,9); imagesc(test0_180S, [0 1]); title('0 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(4,4,10); imagesc(test45_225S, [0 1]); title('45 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(4,4,13); imagesc(test90_270S, [0 1]); title('90 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(4,4,14); imagesc(test135_315S, [0 1]); title('135 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(2,2,4);
    imagesc(ori_map_interpolated); colormap(gca,'hsv'); colorbar;axis square;
    axis off
end




%%  fft ori measurement 
if show_fig == 1
    figure
    subplot(2,4,1); imagesc(test0_180, [0 1]); title('0 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(2,4,2); imagesc(test45_225, [0 1]); title('45 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(2,4,5); imagesc(test90_270, [0 1]); title('90 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(2,4,6); imagesc(test135_315, [0 1]); title('135 Deg'); colormap(gca,'gray'); axis square;axis off
    
    fft0_180 = log(1+abs(fftshift(fft2(test0_180))));
    fft45_225 = log(1+abs(fftshift(fft2(test45_225)))); 
    fft90_270 = log(1+abs(fftshift(fft2(test90_270))));
    fft135_315 = log(1+abs(fftshift(fft2(test135_315)))); 

    subplot(2,4,3); imagesc(fft0_180); title('0 Deg'); colormap(gca,'jet'); axis square;axis off
    subplot(2,4,4); imagesc(fft45_225); title('45 Deg'); colormap(gca,'jet'); axis square;axis off
    subplot(2,4,7); imagesc(fft90_270); title('90 Deg'); colormap(gca,'jet'); axis square;axis off
    subplot(2,4,8); imagesc(fft135_315); title('135 Deg'); colormap(gca,'jet'); axis square;axis off
end

%%      Debug : smoothing Ori map 
if show_fig == 1
    figure
    set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
    subplot(121),imagesc(ori_map_interpolated),colormap('hsv'),colorbar, axis square,  title('Smooth Orientation map','FontSize',20)
    od_contour_levels1 = 1;
    hold on, contour(imresize(ODCrtxPlt_smooth,n_ori_smooth), od_contour_levels1, 'k', 'LineWidth', 5);
    ori_contour_levels = 15;
    subplot(122),imagesc(ori_map_interpolated); colormap hsv; colorbar ;axis square; axis off, title('Orientation map (contours)','FontSize',20)
    hold on ,plot_merged_contours(ori_map_interpolated, ori_contour_levels, 0.2, 'k', 2);
    od_contour_levels2 = 1;
    hold on, contour(imresize(ODCrtxPlt_smooth,n_ori_smooth), od_contour_levels2, 'k', 'LineWidth', 5);
end

%% fft resized for presentation 
if show_fig == 1

    figure
    % raw maps 
    subplot(4,4,1); imagesc(test0_180, [0 1]); title('0 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(4,4,2); imagesc(test45_225, [0 1]); title('45 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(4,4,5); imagesc(test90_270, [0 1]); title('90 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(4,4,6); imagesc(test135_315, [0 1]); title('135 Deg'); colormap(gca,'gray'); axis square;axis off
    
    range_plot = 70; % pixel 
    range_plot_half =  round(range_plot/2) ; 
    fft_plot1 = fft_measure(test0_180); 
    plt_range1 = round(size(fft_plot1,1)/2 - range_plot_half) ; 
    plt_range2 = round(size(fft_plot1,1)/2 + range_plot_half) ; 
    subplot(4,4,3); imagesc(fft_measure(test0_180)); title('0 Deg'); colormap(gca,'jet'); axis square;axis off, xlim([plt_range1, plt_range2]), ylim([plt_range1, plt_range2])
    subplot(4,4,4); imagesc(fft_measure(test45_225)); title('45 Deg'); colormap(gca,'jet'); axis square;axis off, xlim([plt_range1, plt_range2]), ylim([plt_range1, plt_range2])
    subplot(4,4,7); imagesc(fft_measure(test90_270)); title('90 Deg'); colormap(gca,'jet'); axis square;axis off, xlim([plt_range1, plt_range2]), ylim([plt_range1, plt_range2])
    subplot(4,4,8); imagesc(fft_measure(test135_315)); title('135 Deg'); colormap(gca,'jet'); axis square;axis off, xlim([plt_range1, plt_range2]), ylim([plt_range1, plt_range2])
     annotation('textbox',[0.03 0.87 0.98 0.08],'String','Raw maps','EdgeColor','none','fontsize',12)
     
    % interpolated maps 
    subplot(4,4,9); imagesc(test0_180S, [0 1]); title('0 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(4,4,10); imagesc(test45_225S, [0 1]); title('45 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(4,4,13); imagesc(test90_270S, [0 1]); title('90 Deg'); colormap(gca,'gray'); axis square;axis off
    subplot(4,4,14); imagesc(test135_315S, [0 1]); title('135 Deg'); colormap(gca,'gray'); axis square;axis off
    
    fft_plot2 = fft_measure(test0_180S); 
    plt_range3 = round(size(fft_plot2,1)/2 - range_plot_half) ; 
    plt_range4 = round(size(fft_plot2,1)/2 + range_plot_half) ; 
    subplot(4,4,11); imagesc(fft_measure(test0_180S)); title('0 Deg'); colormap(gca,'jet'); axis square; axis off, xlim([plt_range3, plt_range4]), ylim([plt_range3, plt_range4])
    subplot(4,4,12); imagesc(fft_measure(test45_225S)); title('45 Deg'); colormap(gca,'jet'); axis square; axis off, xlim([plt_range3, plt_range4]), ylim([plt_range3, plt_range4])
    subplot(4,4,15); imagesc(fft_measure(test90_270S)); title('90 Deg'); colormap(gca,'jet'); axis square; axis off, xlim([plt_range3, plt_range4]), ylim([plt_range3, plt_range4])
    subplot(4,4,16); imagesc(fft_measure(test135_315S)); title('135 Deg'); colormap(gca,'jet'); axis square; axis off, xlim([plt_range3, plt_range4]), ylim([plt_range3, plt_range4])
    annotation('textbox',[0.03 0.47 0.98 0.08],'String','Interpolated maps','EdgeColor','none','fontsize',12)
    annotation('textbox', [.7 .9 .9 .08], 'string', 'cropped fft', 'edgecolor', 'none', 'fontsize', 12)
end

end 

function result = fft_measure(input)
% Measure fft of input image 
% input     :   input orientation map component 
% result    :   modified fft for presentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fft_input = fft2(input);
    fft_input(1) = 0; % assign dc componenet to zero 
    fft_log = log(1+abs(fftshift(fft_input)));
    fft_log_sharp = fft_log .^ 3; 
    fft_log_intrep = interp2(fft_log_sharp, 2); 
    result = fft_log_intrep; 
end 


