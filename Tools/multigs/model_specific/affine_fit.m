function P = affine_fit(X)

x1 = X(1:3,:);
x2 = X(4:6,:);

H = vgg_Haffine_from_x_MLE(x1,x2);

P = H(:);


end


