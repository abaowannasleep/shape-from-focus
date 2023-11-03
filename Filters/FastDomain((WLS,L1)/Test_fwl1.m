clc 
clear
img = (imread('lamp.jpg')); 

sigma_r=0.1;
lambda=20^2;

X1 = fwl1(double(img) ,double(img), lambda, 5,5, sigma_r); 
X1 = reshape(X1,size(img,1),size(img,2),size(img,3));
figure, imshow(uint8(X1))

