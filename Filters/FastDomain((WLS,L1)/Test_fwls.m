clc 
clear
addpath('graph')

img = (imread('flower.bmp'));
sigma_r=0.03;
lambda=20^2;

X1 = fwls(double(img), double(img), sigma_r,lambda, 5, 1.2);

tic
[m,n,c]=size(img); N=m*n; test=reshape(double((img)),N,c);
[~, edges]=lattice(m,n,0);
weight=(sum((test(edges(:,1),:)-test(edges(:,2),:)).^2,2));
weight=exp(-sqrt(weight)/(sigma_r*255)); weight=weight(:,[1 1 1]);
W=adjacency(edges,weight(:,1)); D=spdiags(sum(W,2),0,N,N);
L=D-W; I=spdiags(ones(N,1),0,N,N);
gt=(I+lambda*L)\test; gt=reshape(gt,[m,n,c]);
toc
figure, imshow([uint8(img) uint8(gt) uint8(X1)]);