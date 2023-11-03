clear all;close all;clc
addpath(genpath('./'))
%% Read the image
I = im2double(imread('colorband.jpg'));
figure;imshow(I)
%% Do the job
alpha_t = .12;
N_iter =10;
mode = 2; %mode = 0 dynamic/dynamic 
          % mode = 1; % static/dynamic
          % mode = 2; % dynamic only
tic
[T,R] = muGIF(I,I,alpha_t,0,N_iter,mode)  ;
toc
figure;imshow(T,[]);colormap('hot');

