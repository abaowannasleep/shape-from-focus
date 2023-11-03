close all;

I = double(imread('tulips.bmp')) / 255;

p = I;
r = 16;
eps = 0.1^2;

q = zeros(size(I));
q(:, :, 1) = guidedfilter(I(:, :, 1), p(:, :, 1), r, eps);
q(:, :, 2) = guidedfilter(I(:, :, 2), p(:, :, 2), r, eps);
q(:, :, 3) = guidedfilter(I(:, :, 3), p(:, :, 3), r, eps);
I_enhanced = (I - q) * 5 + q;

figure('Name','Guided Filtering');
imshow([I, q, I_enhanced], [0, 1]);

q2 = zeros(size(I));
q2(:, :, 1) = gguidedfilter(I(:, :, 1), p(:, :, 1), r, eps);
q2(:, :, 2) = gguidedfilter(I(:, :, 2), p(:, :, 2), r, eps);
q2(:, :, 3) = gguidedfilter(I(:, :, 3), p(:, :, 3), r, eps);
I_enhanced2 = (I - q2) * 5 + q2;

figure('Name','Proposed Filtering');
imshow([I, q2, I_enhanced], [0, 1]);
figure('Name','Enhanced images and detail images');
imshow([I_enhanced, I_enhanced2, 0.5+I-q, 0.5+I-q2], [0, 1]);