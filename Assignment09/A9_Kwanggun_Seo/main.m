%% Automatic Computation of F

close all
clear
clc
%% setup 
% input two images
run('vlfeat-0.9.20/toolbox/vl_setup')
vl_version verbose

I1 = imread('Images(A9)/1.jpg');
I2 = imread('Images(A9)/2.jpg');

I3 = imread('Images(A9)/3.jpg');
I4 = imread('Images(A9)/4.jpg');
% figure
% subplot(2,2,1); imshow(I1);
% subplot(2,2,2); imshow(I2);
% subplot(2,2,3); imshow(I3);
% subplot(2,2,4); imshow(I4);

iteration = 10000;
AutoFundamental(I1,I2,1,iteration)
AutoFundamental(I3,I4,1,iteration)

AutoFundamentalHarrisCorner(I1,I2,1,iteration)
AutoFundamentalHarrisCorner(I3,I4,1,iteration)
