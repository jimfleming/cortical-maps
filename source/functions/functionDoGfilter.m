function mask1 = functionDoGfilter(c1,c2,s1,s2,CSRatio)
% using multivariate 2D mask inhibitory and excitatory effect
% mask of different angles affecting same input

% input parameter
%   Mask1
center_width1 = c1;
center_width2 = c2;
%   Mask2
surround_width1 = s1;
surround_width2 = s2;

orientation = 0; % orientation in degree, the masks parameters must be asymmetric
%orientation =  90 - orientation1;

mu = [0 0]; %center of mask
rx1 = 20; %15
rx2 = 20; %15
x1 = -rx1:rx1;
x2 = -rx2:rx2;
[X1,X2] = meshgrid(x1,-x2);

variance1 = [center_width1 0;0 center_width2];
theta1 = orientation * (pi/180);
rotation_mat1 = [-cos(theta1) sin(theta1);sin(theta1) cos(theta1)];
covariance1 = rotation_mat1 * variance1 * rotation_mat1' ;

F1 = mvnpdf([X1(:) X2(:)], mu, covariance1);
%F1 = center_strength * ( sqrt(det(Sigma1))*2*(pi) ) * reshape(F1,length(x2),length(x1));
F1 =  reshape(F1,length(x2),length(x1));

variance2 = [surround_width1 0; 0 surround_width2];
theta2 = orientation * (pi/180);
rotation_mat2 = [-cos(theta2) sin(theta2);sin(theta2) cos(theta2)];
covariance2 = rotation_mat2 * variance2 * rotation_mat2';

F2 = mvnpdf([X1(:) X2(:)], mu, covariance2);
%F2 = surround_strength * ( sqrt(det(Sigma2))*2*(pi) ) * reshape(F2,length(x2),length(x1));
F2 =  reshape(F2,length(x2),length(x1));
mask1 = (CSRatio*F1-F2) ;
% 
% AreaC = mask1 > 0;  %   figure,imagesc(AreaC)
% AreaS = mask1 < 0;  %   figure,imagesc(AreaS)
% 
% w = CSRatio ;
% mask1 = (mask1 .* AreaC) * w + mask1 .* AreaS; 
