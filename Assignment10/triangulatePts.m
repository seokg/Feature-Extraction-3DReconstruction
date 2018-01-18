function X = triangulatePts(pts1, pts2, camera1, camera2)
%% Triangulation


%% 
X = zeros(4,length(pts1));
for j = 1: length(pts1)
    
        A = computeAtriangulation(pts1(:,j),pts2(:,j),camera1,camera2);
        [~,~,V] = svd(A);
        X(:,j) = V(:,end);
        
end
X(1:4,:) = X(1:4,:) ./ X([4 4 4 4],:);

end

function A = computeAtriangulation(pts1,pts2,P1,P2)
    
A = [pts1(1).*P1(3,:)-P1(1,:);...
    pts1(2).*P1(3,:)-P1(2,:);...
    pts2(1).*P2(3,:)-P2(1,:);...
    pts2(2).*P2(3,:)-P2(2,:);...
    ];
end