clear all;close all;clc
addpath(genpath('./'))

%% Read images
I = imresize(im2double(imread('depth_input.png')),1);
I2 = imresize(im2double(imread('rgb_input.png')),1);
figure;imshow(I,[]);colormap('hot');
figure;imshow(I2)

%% Do the job -- dynamic only on depth
alpha_t = .001;% can be tuned
alpha_r = 0;
N_iter =10;
mode = 2; %mode = 0 dynamic/dynamic 
          % mode = 1; % static/dynamic
          % mode = 2; % dynamic only
tic
[T,~] = muGIF(I,I,alpha_t,alpha_r,N_iter,mode)  ;
toc
figure;imshow(T,[]);colormap('hot');


%% Do the job -- static/dynamic 
alpha_t = .05; % can be tuned
alpha_r = 0;
N_iter =10;
mode = 1; 

tic
[T,~] = muGIF(I2,I,alpha_t,alpha_r,N_iter,mode)  ;
toc
figure;imshow(T,[]);colormap('hot');

%% Do the job -- dynamic/dynamic
alpha_t = .005; % can be tuned
alpha_r = 0.02; % can be tuned
N_iter =10;
mode = 0; 

tic
[T,R] = muGIF(I,I2,alpha_t,alpha_r,N_iter,mode)  ;
toc
figure;imshow(T,[]);colormap('hot');
figure;imshow(R)
