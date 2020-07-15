function [ fitfn resfn degenfn psize numpar ] = getModelParam(model_type)

%---------------------------
% Model specific parameters.
%---------------------------

switch model_type

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