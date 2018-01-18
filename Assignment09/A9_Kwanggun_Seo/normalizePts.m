function [normPts, T] = normalizePts(pts)
%% normalize points         
    meanX = mean(pts(1,:));
    meanY = mean(pts(2,:));
    value = mean(sqrt((pts(1,:) - meanX).^2 + (pts(2,:) - meanY).^2));
    s = sqrt(2) / value;
    tx = - s * meanX;
    ty = - s * meanY;
    
    T = [s 0 tx ;...
         0 s ty ;...
         0 0  1 ];
     
    normPts = T * [pts; ones(1,length(pts))];
end

