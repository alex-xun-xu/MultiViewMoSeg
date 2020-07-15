% Fundamental matrix fitting code adapted from Peter Kovesi's
% implementation of 8-point fundamental matrix estimation. Need 8 or more
% keypoint correspondences.

function F = fundamental_fit(X)

npts = size(X,2);
x1 = X(1:3,:);
x2 = X(4:6,:);

% Build the constraint matrix
A = [x2(1,:)'.*x1(1,:)'   x2(1,:)'.*x1(2,:)'  x2(1,:)' ...
    x2(2,:)'.*x1(1,:)'   x2(2,:)'.*x1(2,:)'  x2(2,:)' ...
    x1(1,:)'             x1(2,:)'            ones(npts,1) ];
[U,D,V] = svd(A,0); % Under MATLAB use the economy decomposition     
  
% Extract fundamental matrix from the column of V corresponding to
% smallest singular value.
F = reshape(V(:,9),3,3)';

% Enforce constraint that fundamental matrix has rank 2 by performing
% a svd and then reconstructing with the two largest singular values.
[U,D,V] = svd(F,0);
F = U*diag([D(1,1) D(2,2) 0])*V';
    
end

