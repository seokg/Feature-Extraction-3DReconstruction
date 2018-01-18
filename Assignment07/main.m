%% Assignment 08 Robust estimation: Automatic computation of H
% EE838 Feature Extration and 3D reconstruction 
% 20164332 Kwanggun Seo

% A. Objective: Compute homography between two images 
%     i. H1-ex1.png ? H1-ex2.png,  H2-ex1.png ? H2-ex2.png 
% B. Interest points: Compute interest points in each image 
% C. Putative correspondences: Compute a set of interest point matches based on some similarity measure 
% D. RANSAC robust estimation: Repeat for N samples 
%     i. Select 4 correspondences and compute H 
%     ii. Calculate the distance for each putative match 
%     iii. Compute the number of inliers consistent with H 
% E. Optimal estimation: Re-estimate H from all inliers by minimizing ML cost function with Levenberg-Marquardt 
% F. Guided matching: Determine more matches using prediction by computed H

%% SETUP
close all, clc, clear

% load image
boardImg1 = imread('01.jpeg');
boardImg2 = imread('03.jpeg');
boardImg3 = imread('02.jpeg');
boardImg4 = imread('04.jpeg');
boardImg1 = imresize(boardImg1,0.2);
boardImg2 = imresize(boardImg2,0.2);
boardImg3 = imresize(boardImg3,0.2);
boardImg4 = imresize(boardImg4,0.2);
boardImg1 = permute(boardImg1,[2 1 3]);
boardImg2 = permute(boardImg2,[2 1 3]);
boardImg3 = permute(boardImg3,[2 1 3]);
boardImg4 = permute(boardImg4,[2 1 3]);
% boardImg1 = imread('Homography/H1_ex1.png');
% boardImg2 = imread('Homography/H1_ex2.png');
% boardImg3 = imread('Homography/H2_ex1.png');
% boardImg4 = imread('Homography/H2_ex2.png');
[y,x,~] = size(boardImg1);
figure
subplot(2,1,1)
imshow([boardImg1 boardImg2]);
title('board image 1 and board image 2 montage')
subplot(2,1,2)
imshow([boardImg3 boardImg4]);
title('board image 3 and board image 4 montage')


% find interest point
% load vl_feat
run('vlfeat-0.9.20/toolbox/vl_setup')
vl_version verbose

%% FEATURE DETECTION: HARRIS CORNER

fa = harrisCorner(boardImg1,10000,9);
fb = harrisCorner(boardImg2,10000,9);
fa2 = harrisCorner(boardImg3,10000,9);
fb2 = harrisCorner(boardImg4,10000,9);

figure
subplot(2,1,1)
imshow([boardImg1 boardImg2]);
hold on 
plot(fa(1,:),fa(2,:),'ro');
plot(fb(1,:)+x,fb(2,:),'ro');
title('features in board image 1 and board image 2 montage')
hold off
subplot(2,1,2)
imshow([boardImg3 boardImg4]);
hold on 
plot(fa2(1,:),fa2(2,:),'ro');
plot(fb2(1,:)+x,fb2(2,:),'ro');
title('features in board image 3 and board image 4 montage')

hold off

%% FEATURE MATCHING: NNC
matches = nncMatching(single(rgb2gray(boardImg1)),single(rgb2gray(boardImg2)),fa,fb,15);
matches2 = nncMatching(single(rgb2gray(boardImg3)),single(rgb2gray(boardImg4)),fa2,fb2,15);

step = length(matches)/50;
step2 = length(matches2)/50;
if step < 1
    step = 1;
end
if step2 < 1
    step2 = 1;
end
figure; 
subplot(2,1,1)
showMatchedFeatures(boardImg1,boardImg2,fa(1:2,matches(1,1:step:end))',fb(1:2,matches(2,1:step:end))','montage');
title('matched features in board image 1 and board image 2 montage')
subplot(2,1,2)
showMatchedFeatures(boardImg3,boardImg4,fa2(1:2,matches2(1,1:step2:end))',fb2(1:2,matches2(2,1:step2:end))','montage');
title('matched features in board image 3 and board image 4 montage')

%% SIFT MATCHING
% % f has a column for each frame. A frame is a disk of center f(1:2), scale f(3) and orientation f(4)
% [fa, da] = vl_sift(single(rgb2gray(boardImg1))) ;
% [fb, db] = vl_sift(single(rgb2gray(boardImg2))) ;
% % matches [image 1 feature point idx; image 2 feature point idx]
% [matches, scores] = vl_ubcmatch(da, db) ;
% 
% % sort matches in descending order by score
% [sortScores, idx]= sort(scores,'descend');
% sortMatches = [matches(1,idx);matches(2,idx)];
% clear idx;
%% RANSAC: compute H
[finalH,finalIdx] = RANSAC_H(fa(1:2,matches(1,:)),fb(1:2,matches(2,:)),5);
% [finalH2,finalIdx2] = RANSAC_H(fa2(1:2,matches2(1,:)),fb2(1:2,matches2(2,:)),5);

p1 = fa(1:2,matches(1,:));
p2 = fb(1:2,matches(2,:));
p3 = fa2(1:2,matches2(1,:));
p4 = fb2(1:2,matches2(2,:));



step = floor(length(finalIdx)/50);
step2 = floor(length(finalIdx)/50);
if step < 1
    step = 1;
end
if step2 < 1
    step2 = 1;
end
%% figure; 
% subplot(2,1,1);
figure
imshow([boardImg1, boardImg2])
hold on 
% outlier 
mask = zeros(1,length(p1));
mask(finalIdx) = 1;
outlierIdx = find(mask==0);
paddedPtsOut = p2(1:2,outlierIdx(1:end))' + [ones(length(outlierIdx),1)*x zeros(length(outlierIdx),1)];
showMatch(p1(:,outlierIdx)',paddedPtsOut,'r')
% inlier
paddedPtsIn = p2(1:2,finalIdx)' + [ones(length(finalIdx),1)*x zeros(length(finalIdx),1)];
showMatch(p1(1:2,finalIdx)',paddedPtsIn,'g')
% showMatchedFeatures(boardImg1,boardImg2,p1(1:2,finalIdx(1:step:end))',p2(1:2,finalIdx(1:step:end))','montage');
title('inlier outlier matches')
hold off 
% showMatchedFeatures(boardImg3,boardImg4,p3(1:2,finalIdx2(1:step2:end))',p4(1:2,finalIdx2(1:step2:end))','montage');

%% infuse image

newCoordinate = finalH * [p1(1:2,finalIdx(1:end)); ones(1,length(finalIdx))];
newCoordinate = newCoordinate./repmat(newCoordinate(3,:),[3,1]);
newCoordinatePrime = p2(1:2,finalIdx(1:end));
output = deformImg(finalH,boardImg1);
output2 = deformImg(inv(finalH),boardImg2);

output = deformImg(finalH,boardImg2);
output2 = deformImg(inv(finalH),boardImg1);
figure
subplot(2,1,1)
imshow([output boardImg2])
title('wrapped image using H')
subplot(2,1,2)
imshow([output2 boardImg1])
title('wrapped image using inverse H')

