function [ fitfn resfn degenfn psize numpar ] = getModelParam(model_type)

%---------------------------
% Model specific parameters.
%---------------------------

switch lower(model_type)

    case 'affine'
        fitfn = @affine_fit;
        resfn = @homography_res;
        degenfn = @affine_degen;
        psize = 3;
        numpar = 9;
    
    case 'homography'
        fitfn = @homography_fit;
        resfn = @homography_res;
        degenfn = @homography_degen;
        psize = 4;
        numpar = 9;
    case 'fundamental'
        fitfn = @fundamental_fit;
        resfn = @fundamental_res;
        degenfn = @fundamental_degen;
        psize = 8;
        numpar = 9;
    otherwise
        error('unknown model type!');
end

end