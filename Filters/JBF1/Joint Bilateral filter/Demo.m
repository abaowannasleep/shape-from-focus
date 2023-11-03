%Demo
%Author: Mahmoud Afifi, York University

close all
clear all
clc
%Bilateral filter
I1=imread('image1.jpg');
I2=I1; %use the same image for traditional bilateral filter
result=bfilter2(I1,I2,15,2,0.15);
figure;
imshow(result);
title('Bilateral filter'); 
imwrite(result,'result1.jpg');

%Joint-bilateral filter
I2=imread('image2.jpg');
result=bfilter2(I1,I2,15,2,0.15);
figure;
imshow(result);
title('Joint-bilateral filter');
imwrite(result,'result2.jpg');

