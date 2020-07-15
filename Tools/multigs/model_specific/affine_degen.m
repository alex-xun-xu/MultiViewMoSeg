function r = affine_degen(X)
        
x1 = X(1:3,:);    % Extract x1 and x2 from x
x2 = X(4:6,:);

r = ...
    iscolinear(x1(:,1),x1(:,2),x1(:,3)) | ...
    iscolinear(x2(:,1),x2(:,2),x2(:,3));

end
