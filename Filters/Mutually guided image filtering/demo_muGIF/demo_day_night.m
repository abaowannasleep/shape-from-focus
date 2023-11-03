clear all;close all;clc
addpath(genpath('./'))

%% Read images
I = imresize(im2double(imread('day.png')),1);
I2 = imresize(im2double(imread('night.png')),1);
figure;imshow(I);
figure;imshow(I2)
%% Mutual structure extraction
alpha_t = .06;
alpha_r = .06;
N_iter =10;
mode = 0; %mode = 0 dynamic/dynamic 
          % mode = 1; % static/dynamic
          % mode = 2; % dynamic only
tic
[T,R] = muGIF(I,I2,alpha_t,alpha_r,N_iter,mode)  ;
toc
figure;imshow(T);
figure;imshow(R);

