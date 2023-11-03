function q = gguidedfilter(I, p, r, eps)
%   Code for
%   [1] "Gradient Domain Guided Image Filtering", Fei Kou, Weihai Chen, Changyun Wen, Zhengguo Li,
%   IEEE Transactions on Image Processing, Vol.24,  No.11, pp. 4528 - 4539, 2015.

%   Code Author: Fei KOU
%   Email: koufei@hotmail.com
%   Date: 9/21/2015
%   Modified from the code provied by Dr. Kaiming He et. al. for their paper 
%   [2] Guided Image Filtering, by Kaiming He, Jian Sun, and Xiaoou Tang, in ECCV 2010 
%  
%   The code and the algorithm are for non-comercial use only.

%   - guidance image: I (should be a gray-scale/single channel image)
%   - filtering input image: p (should be a gray-scale/single channel image)
%   - local window radius: r
%   - regularization parameter: eps

[hei, wid] = size(I);
N = boxfilter(ones(hei, wid), r); % the size of each local patch; N=(2r+1)^2 except for boundary pixels.

mean_I = boxfilter(I, r) ./ N;
mean_p = boxfilter(p, r) ./ N;
mean_Ip = boxfilter(I.*p, r) ./ N;
cov_Ip = mean_Ip - mean_I .* mean_p; % this is the covariance of (I, p) in each local patch.
mean_II = boxfilter(I.*I, r) ./ N;
var_I = mean_II - mean_I .* mean_I;

r2=1;    
N2 = boxfilter(ones(hei, wid), r2); 
mean_I2 = boxfilter(I, r2) ./ N2;
mean_II2 = boxfilter(I.*I, r2) ./ N2;
var_I2 = mean_II2 - mean_I2 .* mean_I2; %the 3*3 variance of image I 

var = (var_I2.*var_I).^0.5;
eps0 = (0.001)^2;
varfinal = (var+eps0)*sum(sum(1./(var+eps0)))/(hei*wid);    %equ.9

minV = min(var(:));
meanV = mean(var(:));
alpha = meanV;  
kk = -4/(minV-alpha);
w = 1-1./(1+exp(kk*(var-alpha)));    %equ.11

a = (cov_Ip +eps*w./varfinal) ./ (var_I + eps./varfinal);   %equ.12
b = mean_p - a .* mean_I;   %equ.13

mean_a = boxfilter(a, r) ./ N;
mean_b = boxfilter(b, r) ./ N;

q = mean_a .* I + mean_b; 
end