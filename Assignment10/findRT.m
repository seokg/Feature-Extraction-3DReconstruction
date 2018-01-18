function [finalP,R,T,P]= findRT(E, camera01,pts,K)
%% compute R T using essential matrix
% out of four solution find one

%% Reset Essential Matrix
[U,~,V] = svd(E);
if det(U*V.') < 0.98
   E = -E;
end
[U,S,V] = svd(E);

d = (S(1,1) + S(2,2))/2;
S = diag([d d 0]);
E = U * S * V';

% check if E is valid
if abs(det(E)) > power(10,-5)
    disp('ERROR: determine of Essential matrix is not 0')
   return 
end

[U,~,V] = svd(E);

u3 = U(:,3);

W=[0 -1 0; ...
    1 0 0; ...
    0 0 1];

Z = [0 1 0; ...
    -1 0 0; ...
     0 0 0];
 
P = zeros(3,4,4); 

% Set Rotation Matrix
R1 = U*W*V.';
R2 = U*W.'*V.';
if det(R1) < 0
    R1 = -R1;
end
if det(R2) < 0
    R2 = -R2;
end
% Tx = U*Z*U.';
% t = -[Tx(3, 2), Tx(1, 3), Tx(2, 1)];
u3 = U* [0,0,1].';

% lamda = 1
P(:,:,1) = K*[R1 u3];
P(:,:,2) = K*[R2 -u3];
% lamda = -1
P(:,:,3) = K*[R2 u3];
P(:,:,4) = K*[R2 -u3];


%% data initializatio 
pts1 = pts(1:2,:);
pts2 = pts(3:4,:);
X = zeros(4,length(pts1),4);
count = zeros(1,4);
%% Triangulate to determine the best P
for i = 1: 4
%     p = P(:,:,i);
    for j = 1: length(pts)
        A = computeAtriangulation(pts1(:,j),pts2(:,j),camera01,P(:,:,i));
        [~,~,V] = svd(A);
        X(:,j,i) = V(:,end);
        
        
    end
    X(1:4,:,i) = X(1:4,:,i)./ X([4 4 4 4],:,i);
    feasiblePtsNum = nnz(X(4,:,i) > 0);
    disp('feasible points')
    disp(feasiblePtsNum);
    Rt = K\P(:,:,i);
    % find out if last element in the image point is positive
    dprd = Rt(3,1:3) * ((X(1:3,:,i) - repmat(Rt(1:3,4),1,size(X,2))));
    % find out if the X is positive or not
    count(i) = sum(X(3,:,i)>0 & dprd > 0);
    disp(count(i))

end

%% 
[~,idx ] = max(count);
finalP = P(:,:,idx);
temp = K\finalP;
R = temp(:,1:3);
T = temp(:,end);


end

function A = computeAtriangulation(pts1,pts2,P1,P2)
    
A = [pts1(1).*P1(3,:)-P1(1,:);...
    pts1(2).*P1(3,:)-P1(2,:);...
    pts2(1).*P2(3,:)-P2(1,:);...
    pts2(2).*P2(3,:)-P2(2,:);...
    ];
end