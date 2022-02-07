function orientim = orient_img_calc_new(normim ,gradientsigma)
normim=normim/max(max(abs(normim)));
[rows,cols] = size(normim);
% Calculating image gradients
sze = fix(6*gradientsigma);
if ~mod(sze,2);
    sze = sze+1;
end
f = fspecial('gaussian', sze, gradientsigma); % Generate Gaussian filter.
[fx,fy] = gradient(f);                        % Gradient of Gausian.

Gx = filter2(fx,normim); % Gradient of the image in x
Gy = filter2(fy,normim); % and y
% Estimate the local orientation of each block
%
%Gy=fy;

%f=1;
%f=ones(60,60);
%Gx=fx;
f = fspecial('gaussian', 6*22, 22); % Generate Gaussian filter.
orientim= atan2(filter2(f,2*Gy.*Gx),filter2(f,Gy.^2-Gx.^2))/2;%-pi/2;



%orientim=mod(orientim,pi)+pi/2;

