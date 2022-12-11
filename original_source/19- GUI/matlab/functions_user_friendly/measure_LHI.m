function LHI = measure_LHI(newOriMap1,sigmamn)

mnppx = 50;
%sigmamn = 45; % microns ---standard deviation for LHI calculation

distx = 1:(sigmamn/mnppx)*6;
distx = distx - ceil(length(distx)/2);
[X,Y] = meshgrid(distx,distx);
dist = sqrt(X.^2+Y.^2);
distW = exp(-(dist.^2)./(2*(sigmamn/mnppx)^2));%    figure,imagesc(distW>0.01),colorbar
mult = 1/sum(distW(:));

temp1 = conv2(exp(newOriMap1*(pi/180)*2*1i), distW, 'same'); %,'valid'); %  figure,imagesc(LHI),colorbar
LHI = abs(temp1)*mult;