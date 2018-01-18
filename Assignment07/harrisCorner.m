function [pt] = harrisCorner(img,thresh,window)
%% feature extraction from the image
k = 0.04;
pt = [];
% sobel operator
grayImg = single(rgb2gray(img));
gradX = funcSobel(grayImg,'x');
gradY = funcSobel(grayImg,'y');

% Gaussien Filter
g = fspecial('gaussian'); 

Sx2 = conv2(gradX.^2, g, 'same');  
Sy2 = conv2(gradY.^2, g, 'same');
Sxy = conv2(gradX.*gradY, g,'same');

[row, col, ~] = size(img);
outputImg = zeros(row,col);
R = zeros(row,col);
for i = 1:row
    for j = 1:col
       
        H = [Sx2(i, j) Sxy(i, j); Sxy(i, j) Sy2(i, j)];
        R(i,j) = det(H) - k*trace(H).^2;
        if R(i,j) > thresh
           outputImg(i,j) = R(i,j);
           pt = [pt [j;i]];
        end
    end 
end
% non-max suppression
mask= ones(window,window);
mask(ceil(window/2),ceil(window/2)) = 0;
B = ordfilt2(outputImg,80,mask);
outputImg = outputImg > B;

% imshow(outputImg);


end
