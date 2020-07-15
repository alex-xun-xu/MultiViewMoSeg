% Calculate the residuals for Fundamental matrix estimation. The code is
% adapted from Peter Kovesi's implementation of Fundemental matrix
% estimation
% (http://www.csse.uwa.edu.au/~pk/research/matlabfns/Robust/ransacfitfundma
% trix7.m).
%
% Note that this code allows for F being a cell array of fundamental matrices of
% which we have to pick the best one. (A 7 point solution can return up to 3
% solutions)

function [dist, P] = fundamental_res(F,X)

x1 = X(1:3,:);    % Extract x1 and x2 from x
x2 = X(4:6,:);
n = size(x1,2);
t = 0.1;

if isempty(F)
    dist = inf*ones(n,1);
    P = [];
    return;
end

if iscell(F)  % We have several solutions each of which must be tested
    
    nF = length(F);   % Number of solutions to test
    P = F{1}(:);      % Initial allocation of best solution
    ninliers = 0;     % Number of inliers
    dist = zerosz(1,n);
    
    for k = 1:nF
        x2tFx1 = zeros(1,length(x1));
        for n = 1:length(x1)
            x2tFx1(n) = x2(:,n)'*F{k}*x1(:,n);
        end
        
        Fx1 = F{k}*x1;
        Ftx2 = F{k}'*x2;
        
        % Evaluate distances
        d = abs(x2tFx1.^2 ./ ...
            (Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2)); 
        
        inliers = find(abs(d) < t);     % Indices of inlying points
        
        if length(inliers) > ninliers   % Record best solution
            ninliers = length(inliers);
            dist = d';
            P = F{k}(:);
        end
    end
    
else     % We just have one solution
    % Ensure the Fundamental matrix is of size 3x3.
    F = reshape(F,3,3);
    x2tFx1 = zeros(1,length(x1));
    for n = 1:length(x1)
        x2tFx1(n) = x2(:,n)'*F*x1(:,n);
    end
    
    Fx1 = F*x1;
    Ftx2 = F'*x2;
    
    % Evaluate distances
    d =  abs(x2tFx1.^2 ./ ...
        (Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2));
    dist = d';
    
    P = F(:); 
end

end