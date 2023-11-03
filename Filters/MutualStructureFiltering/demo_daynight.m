close all
clear

day = im2double(imread('testdata/day_night/day.png'));
night = im2double(imread('testdata/day_night/night.png'));

day = imresize(day,0.6);
night = imresize(night,0.6);

eps_I = 1e-3;
eps_G = 2e-4;
lambda_I = 3;
lambda_G = 1;
maxiter = 20;
r = 4;

[I, G] = mutual_structure_joint_filtering(day, night , r, eps_I, eps_G, lambda_I, lambda_G, maxiter);

figure, imshow(I)
figure, imshow(G)