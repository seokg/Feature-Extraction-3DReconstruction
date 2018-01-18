%% Assignment 10
% Two View 3D Reconstruction
clc 
close all
clear

showImg = false;
%% load image / data
run('vlfeat-0.9.20/toolbox/vl_setup')
vl_version verbose

I1 = imread('A10/1.jpg');
I2 = imread('A10/2.jpg');

if showImg
imshow([I1 I2]);
end

%-- Focal length:
fc = [ 2892.843725329502400 ; 2882.249450476587300 ];

%-- Principal point:
cc = [ 824.425157504919530 ; 605.187152104484080 ];

%-- Skew coefficient:
alpha_c = 0.000000000000000;

K = [diag(fc) cc; zeros(1,3)];
K(end) = 1;
K(1,2) = alpha_c;
%% compute F 
disp('Computing F')
iteration = 500;
[F,pts] = AutoFundamental(I1,I2,5,iteration);
%% compute E
E = K.' * F * K;


 
%% find R T
disp('Finding R T ...')
camera =  zeros(3,4,2);
camera(:,:,1) = K * [diag(ones(3,1)), zeros(3,1)];
% motionFromF(F,K,pts(1:2,:)',pts(3:4,:)');
[camera(:,:,2),R,T,P] = findRT(E, camera(:,:,1), pts,K);
%% Triangulate 01
disp('Triangulating ...')
X = triangulatePts(pts(1:2,:), pts(3:4,:) ,camera(:,:,1),camera(:,:,2));
% X = triangulatePts(normalizePts(pts(1:2,:)), normalizePts(pts(3:4,:)) ,[diag(ones(3,1)), zeros(3,1)],[R,T]);

%% Make Point Cloud
disp('Making PLY ...')
rgb = zeros(3,length(pts));
for i = 1: length(pts)
rgb(:,i) = reshape(I1(round(pts(2,i)),round(pts(1,i)),:),[3,1]);
end
plyWrite(X(1:3,:), rgb, 'output_true.ply');

%% Triangulate 01
disp('Triangulating ...')
X = triangulatePts(pts(1:2,:), pts(3:4,:) ,camera(:,:,1),P(:,:,1));
% X = triangulatePts(normalizePts(pts(1:2,:)), normalizePts(pts(3:4,:)) ,[diag(ones(3,1)), zeros(3,1)],[R,T]);

%% Make Point Cloud
disp('Making PLY ...')
rgb = zeros(3,length(pts));
for i = 1: length(pts)
rgb(:,i) = reshape(I1(round(pts(2,i)),round(pts(1,i)),:),[3,1]);
end
plyWrite(X(1:3,:), rgb, 'output1.ply');
%% Triangulate 02
disp('Triangulating ...')
X = triangulatePts(pts(1:2,:), pts(3:4,:) ,camera(:,:,1),P(:,:,2));
% X = triangulatePts(normalizePts(pts(1:2,:)), normalizePts(pts(3:4,:)) ,[diag(ones(3,1)), zeros(3,1)],[R,T]);

%% Make Point Cloud
disp('Making PLY ...')
rgb = zeros(3,length(pts));
for i = 1: length(pts)
rgb(:,i) = reshape(I1(round(pts(2,i)),round(pts(1,i)),:),[3,1]);
end
plyWrite(X(1:3,:), rgb, 'output2.ply');

disp('two view reconstruction finished!!!')

%% Triangulate 03
disp('Triangulating ...')
X = triangulatePts(pts(1:2,:), pts(3:4,:) ,camera(:,:,1),P(:,:,3));
% X = triangulatePts(normalizePts(pts(1:2,:)), normalizePts(pts(3:4,:)) ,[diag(ones(3,1)), zeros(3,1)],[R,T]);

%% Make Point Cloud
disp('Making PLY ...')
rgb = zeros(3,length(pts));
for i = 1: length(pts)
rgb(:,i) = reshape(I1(round(pts(2,i)),round(pts(1,i)),:),[3,1]);
end
plyWrite(X(1:3,:), rgb, 'output3.ply');

%% Triangulate 04
disp('Triangulating ...')
X = triangulatePts(pts(1:2,:), pts(3:4,:) ,camera(:,:,1),P(:,:,3));
% X = triangulatePts(normalizePts(pts(1:2,:)), normalizePts(pts(3:4,:)) ,[diag(ones(3,1)), zeros(3,1)],[R,T]);

%% Make Point Cloud
disp('Making PLY ...')
rgb = zeros(3,length(pts));
for i = 1: length(pts)
rgb(:,i) = reshape(I1(round(pts(2,i)),round(pts(1,i)),:),[3,1]);
end
plyWrite(X(1:3,:), rgb, 'output4.ply');


