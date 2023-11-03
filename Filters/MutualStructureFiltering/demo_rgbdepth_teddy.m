close all
clear

rgb = im2double(imread('testdata/rgb_depth_teddy/rgb_input.png'));
depth = im2double(imread('testdata/rgb_depth_teddy/depth_input.png'));

eps_I = 1e-4;
eps_G = 5e-4;
lambda_I = 1;
lambda_G = 30;
maxiter = 20;
r = 2;
[I, G] = mutual_structure_joint_filtering(rgb, repmat(depth,[1 1 3]), r, eps_I, eps_G, lambda_I, lambda_G, maxiter);

figure, imshow(I)
figure, imshow(medfilt2(rgb2gray(G),[3 3]))