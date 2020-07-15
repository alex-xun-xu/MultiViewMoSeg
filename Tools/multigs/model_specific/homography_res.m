function [dist, H] = homography_res(H, X)
    
    H = reshape(H,3,3);
   
    
    x1 = X(1:3,:);   % Extract x1 and x2 from x
    x2 = X(4:6,:);    
    n = size(x1,2);
    
    % Calculate, in both directions, the transfered points    
    Hx1    = H*x1;
    invHx2 = H\x2;
    
    % Normalise so that the homogeneous scale parameter for all coordinates
    % is 1.
    
    x1     = hnormalise(x1);
    x2     = hnormalise(x2);     
    Hx1    = hnormalise(Hx1);
    invHx2 = hnormalise(invHx2); 
    
    dist = sum((x1-invHx2).^2)  + sum((x2-Hx1).^2);
    
    dist = reshape(dist,n,1);
    
    H = H(:);

end

