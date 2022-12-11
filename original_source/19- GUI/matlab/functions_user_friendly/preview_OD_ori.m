function [OD, ori] = preview_OD_ori(OD, ori, kern_std, n_interp, smooth_OD_border, ori_contour_levels, label_contours)
%     close all;
    OD = denoise_OD(OD, kern_std, n_interp, smooth_OD_border);
%     figure(1); imagesc(round(OD)); colormap(jet); hold on;
%     c_od = contour(OD, [0.5 0.5], 'w');
%     if label_contours
%         clabel(c_od);
%     end
%     if ~isreal(ori)
%         ori = angle(ori) ./ 2 .* (180./pi) + 90;
%     end
%     figure(2); imagesc(ori); colormap(hsv); hold on;
%     c_ori = contour(ori, ori_contour_levels, 'k', 'LineWidth', 2);
%     if label_contours
%         clabel(c_ori);
%     end
end

function OD_interp = denoise_OD(OD, kern_std, n_interp, use_raw_values)
    if kern_std
        kern = fspecial('gaussian', kern_std*6, kern_std);
        OD_smoothed = conv2(OD, kern, 'same');
        if use_raw_values
            OD_smoothed = (OD_smoothed - min(min(OD_smoothed))) ./ (max(max(OD_smoothed)) - min(min(OD_smoothed)));
            OD_interp = interp2(OD_smoothed, n_interp);
        else
            OD_smoothed = round((OD_smoothed - min(min(OD_smoothed))) ./ (max(max(OD_smoothed)) - min(min(OD_smoothed))));
            OD_interp = round(interp2(OD_smoothed, n_interp));
        end
    else
        OD_interp = OD;
    end
end