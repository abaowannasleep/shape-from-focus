clear, clc, close all
addpath(genpath('Functions'));
addpath(genpath('Filters'));
%%
focus = (0.02: 0.02: 0.88);
simu_stone = zeros(2000, 2000, size(focus, 2));
for i = 1: size(focus, 2)
    if i<10
        imagePath = strcat('./simu_data/im_0', num2str(i),'.bmp');
    else
        imagePath = strcat('./simu_data/im_', num2str(i),'.bmp'); 
    end
    im = imread(imagePath); 
    simu_stone(:, :, i) = double(im);
end
fms = GLVM(simu_stone, 9, 0);
[I, z, s, A] = gauss3P(focus, fms);
%% err
ssdErr = ssd(fms, focus, A, z, s);
niccErr = nicc(fms, simu_stone);
ksErr = kstest(fms, focus, I, A, z, s);
%% Carve depthmap by removing error pixels:
mask = imread('./simu_data/mask.png');
mask(mask==255) = 1;
mask = double(mask);
z0 = z .* mask;

% r1
mask1 = ones(size(z));
mask1(ksErr>35) = 0;
mask1 = mask1.*mask;
z1 = z .* mask1;
% r2
mask2 = ones(size(z));
mask2(isnan(niccErr)) = 0;
mask2(niccErr<0) = 0;
mask2 = mask2.*mask;
z2 = z .* mask2;
% r_ours
mask3 = depthFilter(z);
mask3(ssdErr>1000) = 0;
mask3 = mask3 .* mask;
z3 = z .* mask3;

%% accuracy
t = 0.0172;
m = abs(z0 - dmap);
m1 = abs(z1 - dmap);
m2 = abs(z2 - dmap);
m3 = abs(z3 - dmap);

vp = zeros(size(m));
vp1 = zeros(size(m));
vp2 = zeros(size(m));
vp3 = zeros(size(m));
vp(m>t) = 1;
vp1(m1>t) = 1;
vp2(m2>t) = 1;
vp3(m3>t) = 1;
% vp(m<=t) = 1;
% vp1(m1<=t) = 1;
% vp2(m2<=t) = 1;
% vp3(m3<=t) = 1;

sum(sum(vp.*mask))/sum(sum(mask))
sum(sum(vp1.*mask1))/sum(sum(mask))
sum(sum(vp2.*mask2))/sum(sum(mask))
sum(sum(vp3.*mask3))/sum(sum(mask))
%% rmse
e = rmse(z0, dmap, mask)
e1 = rmse(z1, dmap, mask1)
e2 = rmse(z2, dmap, mask2)
e3 = rmse(z3, dmap, mask3)
% ee = rmse(z_result, dmap, mask)

%% me
e = me(z0, dmap, mask)
e1 = me(z1, dmap, mask1)
e2 = me(z2, dmap, mask2)
e3 = me(z3, dmap, mask3)
% ee = me(z_result, dmap, mask)

%% Display the result:
displayR1(z0, z1, z2, z3, mask, mask1, mask2, mask3);
%%

alpha = 0.01;
d_max = max(max(z3));
d_min = min(min(z3));
d = d_min + (d_max-d_min) * rand(size(z3), 'double');
d = d .* mask; 

zGD = z3;
mask_gd = mask3;
for i = 1:3
    [zGD, err] = gradient_descent(zGD, d, mask_gd, alpha);
    mask_gd = depthFilter(zGD);
    zGD = zGD .* mask_gd;
end

%% other three filters
% median filter
zMF = medfilt2(z3, [9 9]);
zMF = zMF.*mask;

% WGIF 【code from https://data.mendeley.com/datasets/69n24h6ydh/1】 
g = mask3;
radius=2;
zWGIF=WeightedGuidedImageFilter(z3, g, radius, 0.10, 1);

%% rmse
e = rmse(z0, dmap, mask)
e = rmse(zGD, dmap, mask)
e = rmse(zMF, dmap, mask)
e = rmse(zWGIF, dmap, mask)

%% me
e = me(z0, dmap, mask)
e = me(zGD, dmap, mask)
e = me(zMF, dmap, mask)
e = me(zWGIF, dmap, mask)

%% Display the result2:
displayR2(z0, zMF, zWGIF, zGD, mask);
