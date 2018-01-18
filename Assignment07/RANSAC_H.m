function [finalH,finalIdx] = RANSAC_H(pts1,pts2,threshold)

p = 0.99; % at least one point is inlier
% w = 0.5; % prob inlier ratio
minNumPts = 4;
% N = log(1-p)/log(1-w^minNumPts);
N = 10000;

% initialize
iter = 1;
maxNumInlier = -1;
minStd = 100000;
finalH = [];
finalIdx = [];
rng('shuffle','twister');

% homogeneous coordinate
ptsH1 = [pts1; ones(1,length(pts1))];
ptsH2 = [pts2; ones(1,length(pts2))];
% RANSAC
while iter  < ceil(N)
    % select radom point
    randIdx = ceil(rand([1,4])*length(pts1));
    pts1Sel = pts1(:,randIdx);
    pts2Sel = pts2(:,randIdx);
    
    
%     checkColinearity1 = zeros(4,1);
%     checkColinearity2 = zeros(4,1);
%     for idx = 0:3
%        checkColinearity1(idx+1) = det([...
%             pts1Sel(1,mod(idx,4)+1) pts1Sel(2,mod(idx,4)+1) 1;...
%             pts1Sel(1,mod(idx+1,4)+1) pts1Sel(2,mod(idx+1,4)+1) 1;...
%             pts1Sel(1,mod(idx+2,4)+1) pts1Sel(2,mod(idx+2,4)+1) 1]);
%        checkColinearity2(idx+1) = det([...
%             pts2Sel(1,mod(idx,4)+1) pts2Sel(2,mod(idx,4)+1) 1;...
%             pts2Sel(1,mod(idx+1,4)+1) pts2Sel(2,mod(idx+1,4)+1) 1;...
%             pts2Sel(1,mod(idx+2,4)+1) pts2Sel(2,mod(idx+2,4)+1) 1]);  
%     end
%     if nnz(abs(checkColinearity1) < 1e-2)>0 || nnz(abs(checkColinearity2) < 1e-2)>0
%         disp('points are colinear!!!')
%         continue
%     end    
%     
    % normalize points
    [pts1Sel T1] = normalizePts(pts1Sel);
    [pts2Sel T2] = normalizePts(pts2Sel);
    
    
    % compute model using DLT
    A = computeA(pts1Sel,pts2Sel);
    [U,S,V] = svd(A);
    vectorH = double(vpa(V(:,end)));
%     vectorH = V(:,end);
    H = reshape(vectorH,[3,3])';
    
    H = inv(T2) * H * T1;
%     H = H / H(end);
    % compute error
    % change into homogeneous equation

    
    tempPts2hat = H * ptsH1;
    pts2hat = tempPts2hat ./ repmat(tempPts2hat(3,:),[3,1]);
    tempPts1hat = H\ptsH2;
    pts1hat = tempPts1hat ./ repmat(tempPts1hat(3,:),[3,1]);
    
    [std,distance1,distance2] = standardDeviation(ptsH1,ptsH2,pts1hat,pts2hat);
    [numInlier, inlierIdx]= computeInlier(distance1,distance2,1);
    
    % checking for best possible model
    if numInlier > maxNumInlier ||  (numInlier == maxNumInlier && std < minStd)
        disp('updating H matrix')
       finalH = H;
%        finalH = inv(T2)*H*T1;
       finalIdx = inlierIdx; 
       ptsH1(:,finalIdx);
       ptsH2(:,finalIdx);
       pts1hat(:,finalIdx);
       pts2hat(:,finalIdx);
       maxNumInlier = numInlier;
       minStd = std;
       
       w = numInlier / length(pts1);
       if w ~=0
           N = log(1-p)/log(1-w^minNumPts);
           fprintf('inlier: %i\n',maxNumInlier);
           fprintf('iteration: %i\n',iter);
           fprintf('updated number of iteration: %i\n',N);
       end
    end
    iter = iter + 1;
    
    %% optimize H from the inlier
    if iter  > (N)
%        disp('unoptimized H')
%        disp(finalH)
        
        [normPts1, T1] = normalizePts(pts1(:,finalIdx));
        [normPts2, T2] = normalizePts(pts2(:,finalIdx));
        
        
        A = computeA(normPts1,normPts2);
        [U,S,V] = svd(A);
        vectorH = double(vpa(V(:,end)));
        finalH = reshape(vectorH,[3,3])';
        finalH = inv(T2) * finalH * T1;  
        
        tempPts2hat = finalH * ptsH1;
        pts2hat = tempPts2hat ./ repmat(tempPts2hat(3,:),[3,1]);
        tempPts1hat = finalH\ptsH2;
        pts1hat = tempPts1hat ./ repmat(tempPts1hat(3,:),[3,1]);
        [std,distance1,distance2] = standardDeviation(ptsH1,ptsH2,pts1hat,pts2hat);
        [numInlier, finalIdx]= computeInlier(distance1,distance2,1);
        disp('final inlier')
        disp(numInlier)
        disp('optimized H')
        disp(finalH)

    end
        
end



end


function A = computeA(pts1,pts2)
%% compute A for DLT algorithm
% points are formated as
% [x1 x2 x3 x4
%  y1 y2 y3 y4]
A = [];
for i = 1:length(pts1)
A = [A;
     pts1(1,i) pts1(2,i) 1 0 0 0 -(pts1(1,i) * pts2(1,i)) -(pts2(1,i) * pts1(2,i)) -pts2(1,i);
     0 0 0 pts1(1,i) pts1(2,i) 1 -(pts1(1,i) * pts2(2,i)) -(pts2(2,i) * pts1(2,i)) -pts2(2,i)];
end

end

function [normPts, T] = normalizePts(pts)
%% normalize points         
    meanX = sum(pts(1,:)) / length(pts);
    meanY = sum(pts(2,:)) / length(pts);
    value = sum(sqrt((pts(1,:) - meanX).^2 + (pts(2,:) - meanY).^2)) / length(pts);
    s = sqrt(2) / value;
    tx = - s * meanX;
    ty = - s * meanY;
    
    T = [s 0 tx ;...
         0 s ty ;...
         0 0  1 ];
     
     normPts = T * [pts; ones(1,length(pts))];
end

function [std,distance1,distance2] = standardDeviation(ptsH1,ptsH2,pts1hat,pts2hat)
   
%     totalError = sum( sum( (pts2hat - ptsH2).^2 ,1) ) + sum( sum( (pts1hat - ptsH1).^2 ,1) );
    distance2 = ( sum( (pts2hat - ptsH2).^2 ,1) ).^0.5;
    distance1 = ( sum( (pts1hat - ptsH1).^2 ,1) ).^0.5;

    % comput the standard deviation of the inlier distance
    avg = sum(distance1 + distance2) / length(pts1hat);
    std = (sum((avg - (distance1 + distance2)).^2) / length(pts1hat)).^0.5;
%     fprintf('standard deviation: %i\n', std);
    

end

function [numInlier, inlierIdx]= computeInlier(distance1,distance2,distThreshold)
%% compute the inlier points given threshold     
% compute inlier number    
    inlierTrueMatrix = (distance1 < distThreshold) & (distance2 < distThreshold);
    inlierIdx = find(inlierTrueMatrix==1);
    numInlier = length(inlierIdx);
%     fprintf('number of inlier: %i\n', numInlier);

end