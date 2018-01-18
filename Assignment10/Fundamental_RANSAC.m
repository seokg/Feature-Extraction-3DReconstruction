function[finalF, finalIdx] = Fundamental_RANSAC(p,q,thresh,iteration)
%% compute the fundmental matrix using ransac
% input 
%   p,q: points from image 1 and 2
%   thres: threshold value
% output
%   F: fundamental matrix
%   inlierIdx: inlier index

%% computing iteration
prob = 0.99; % at least one point is inlier
w = 0.01; % prob inlier ratio
minNumPts = 7;
N = log(1-prob)/log(1-w^minNumPts);

%% initialize
iter = 1;
maxNumInlier = -1;
finalF = [];
finalIdx = [];
rng('shuffle','twister');

% change to homogeneous coordinate
pH = [p; ones(1,length(p))];
qH = [q; ones(1,length(q))];

%% RANSAC
while iter  < iteration
    if iter > N
        [finalF,finalIdx] = optimizeF(pH,qH,finalIdx,thresh);
        break;
    end
    % select radom point
    randIdx = ceil(rand([1,minNumPts])*length(p));
    pSel = p(:,randIdx);
    qSel = q(:,randIdx);
    [pSel, T1] = normalizePts(pSel);
    [qSel, T2] = normalizePts(qSel);     
    
    % compute model using DLT
    A = computeA(pSel,qSel);
    [~,~,V] = svd(A,0);
    FF{1} = reshape(V(:,end-1),[3 3]);
    FF{2} = reshape(V(:,end  ),[3 3]);
    
    a = vgg_singF_from_FF(FF);
    F = [];
    for i = 1:length(a)
      Fi = a(i)*FF{1} + (1-a(i))*FF{2};
      if signs_OK(Fi,pSel,qSel)
        F = cat(3, F, Fi);
      end
      
    end
    
    % continue when there i non present
    if isempty(F)
        continue
    end
    
    
    [~,~,chan]= size(F);
    inliertemp = -1; 
    for i = 1:chan
        % unnormalize F
        Ftemp = T2' * F(:,:,i)' * T1;
        Ftemp = Ftemp / norm(Ftemp);
        if Ftemp(end) < 0
          Ftemp = -Ftemp;
        end
        
        % compute error
        [~, dist] = sampsonError(Ftemp, pH, qH);
        inlier = nnz(dist<thresh);
        if inlier > inliertemp
            inliertemp = inlier;
            idxtemp = find(dist<thresh);
            trueF = Ftemp;
        end
    end

        
    % checking for best possible model
    if inliertemp > maxNumInlier
       disp('updating F matrix')
       disp(trueF)
       finalF = trueF;
       finalIdx = idxtemp;
       maxNumInlier = inliertemp;
       
       w = maxNumInlier / length(p);
       if w ~=0
           N = log(1-prob)/log(1-w^minNumPts);
           fprintf('inlier: %i\n',maxNumInlier);
           fprintf('iteration: %i\n',iter);
           fprintf('updated number of iteration: %i\n',N);
       end
    end
    iter = iter + 1;
    
    %% optimize F from the inlier
    if iter  > (N)
        [finalF,finalIdx] = optimizeF(pH,qH,finalIdx,thresh);
    end
        
end




end

function A = computeA(p,q)
%% compute A for DLT algorithm
% points are formated as
% [x1 x2 x3 x4
%  y1 y2 y3 y4]
A = [];
for i = 1:length(p)
A = [A;
     p(1,i)*q(1,i), p(2,i)*q(1,i), q(1,i), p(1,i)*q(2,i), p(2,i)*q(2,i), q(2,i), p(1,i), p(2,i) 1];
end

end
function [finalF,finalIdx] = optimizeF(p,q,finalIdx,thresh)
        [normPts1, T1] = normalizePts(p(1:2,finalIdx));
        [normPts2, T2] = normalizePts(q(1:2,finalIdx));
        
        A = computeA(normPts1,normPts2);
        [~,~,V] = svd(A);
        vectorF = double(vpa(V(:,end)));
        finalF = reshape(vectorF,[3,3])';
        finalF = T2' * finalF * T1;  
        
        finalF = finalF / norm(finalF);
        if finalF(end) < 0
          finalF = -finalF;
        end
        disp('final F')
        disp(finalF)
        % compute error
        [~, dist] = sampsonError(finalF, p, q);
        finalIdx = find(dist<thresh);

end
function [error,dist] = sampsonError(F,p,q)
dist = zeros(1,length(p));
for i = 1 : length(p)
   upper = (q(:,i).' * F * p(:,i)).^2;
   Fp = F*p(:,i);
   Fq = F*q(:,i);
   lower = Fp(1).^2+Fp(2).^2+Fq(1).^2+Fq(2).^2;
   dist(i) = upper/lower; 
    
end
error = sum(dist);

end
%% used from the external source
% http://www.robots.ox.ac.uk/~vgg/hzbook/code/
function OK = signs_OK(F,x1,x2)
[u,s,v] = svd(F');
e1 = v(:,3);
l1 = vgg_contreps(e1)*x1;
s = sum( (F*x2) .* l1 );
OK = all(s>0) | all(s<0);
return
end
function a = vgg_singF_from_FF(F)

% precompute determinants made from columns of F{1}, F{2}
for i1 = 1:2
  for i2 = 1:2
    for i3 = 1:2
      D(i1,i2,i3) = det([F{i1}(:,1) F{i2}(:,2) F{i3}(:,3)]);
    end
  end
end

% Solve The cubic equation for a
a = roots([-D(2,1,1)+D(1,2,2)+D(1,1,1)+D(2,2,1)+D(2,1,2)-D(1,2,1)-D(1,1,2)-D(2,2,2)
            D(1,1,2)-2*D(1,2,2)-2*D(2,1,2)+D(2,1,1)-2*D(2,2,1)+D(1,2,1)+3*D(2,2,2)
            D(2,2,1)+D(1,2,2)+D(2,1,2)-3*D(2,2,2)
            D(2,2,2)]);
a = a(abs(imag(a))<10*eps);

return
end
function Y = vgg_contreps(X)

% vgg_contreps  Contraction with epsilon tensor.
%
% B = vgg_contreps(A) is tensor obtained by contraction of A with epsilon tensor.
% However, it works only if the argument and result fit to matrices, in particular:
%
% - if A is row or column 3-vector ...  B = [A]_x
% - if A is skew-symmetric 3-by-3 matrix ... B is row 3-vector such that A = [B]_x
% - if A is skew-symmetric 4-by-4 matrix ... then A can be interpreted as a 3D line Pluecker matrix
%                                               skew-symmetric 4-by-4 B as its dual Pluecker matrix.
% - if A is row 2-vector ... B = [0 1; -1 0]*A', i.e., A*B=eye(2)
% - if A is column 2-vector ... B = A'*[0 1; -1 0], i.e., B*A=eye(2)
%
% It is vgg_contreps(vgg_contreps(A)) = A.

% werner@robots.ox.ac.uk, Oct 2001

if prod(size(X)) == 3  % get [X]_\times
  Y = [0 X(3) -X(2)
      -X(3) 0 X(1)
       X(2) -X(1) 0];
elseif all(size(X) == [1 2])
  Y = [0 1; -1 0]*X';
elseif all(size(X) == [2 1])
  Y = X'*[0 1; -1 0];
elseif all(size(X) == [3 3]) % get X from [X]_\times
  Y = [X(2,3) X(3,1) X(1,2)];
elseif all(size(X) == [4 4])  % pluecker matrix dual
  Y = [0      X(3,4) X(4,2) X(2,3) 
       X(4,3) 0      X(1,4) X(3,1)
       X(2,4) X(4,1) 0      X(1,2)
       X(3,2) X(1,3) X(2,1) 0     ];
else
  error('Wrong matrix size.')
end
end

