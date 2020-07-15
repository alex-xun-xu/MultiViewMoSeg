function P = homography_fit(X)

x1 = X(1:3,:);
x2 = X(4:6,:);

H = vgg_H_from_x_lin(x1,x2);

P = H(:);


end


