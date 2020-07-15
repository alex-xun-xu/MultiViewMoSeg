%-------------------------------------------------------------------------
% Accelerated Hypothesis Generation for Multi-Structure Robust Fitting
%-------------------------------------------------------------------------
% The demo code in this package implements the guided-sampling method for
% multi-structure robust fitting proposed in:
%
% T.-J. Chin, J. Yu and D. Suter
% Accelerated Hypothesis Generation for Multi-Structure Robust Fitting
% In Proc. European Conf. on Computer Vision, Crete, Greece, 2010.
%
% T.-J. Chin, J. Yu and D. Suter
% Accelerated Hypothesis Generation for Multi-Structure Data via Preference Analysis
% To appear in IEEE Trans. on Pattern Analysis and Machine Intelligence.
%
% Copyright (c) 2010 Tat-Jun Chin and Jin Yu
% School of Computer Science, The University of Adelaide, South Australia
% http://www.cs.adelaide.edu.au/~{tjchin,jinyu}
%
% The program is free for non-commercial academic use. Any commercial use
% is strictly prohibited without the authors' consent. Please acknowledge
% the authors by citing the above paper in any academic publications that
% have made use of this package or part of it.
%
% If you encounter any problems or questions please email to 
% tjchin@cs.adelaide.edu.au.
% 
% This program makes use of Peter Kovesi and Andrew Zisserman's MATLAB
% functions for multi-view geometry
% (http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/
%  http://www.robots.ox.ac.uk/~vgg/hzbook/code/).


function [ par res inx tim ] = randomSampling_AcrossSets(lim,data,set1idx,set2idx,M,model_type)
% input: 
% lim (1x1) = Maximum CPU seconds allowed.
% data (dxn) = Input data of dimensionality d.
% M (1x1) = Maximum number of hypotheses to be generated.
% model_type (string) = Type of model to be estimated.
%
% output:
% par (dxM) = Parameters of the putative models.
% res (nxM) = Residuals as measured to the putative models.
% inx (pxM) = Indices of p-subsets
% tim (1xM) = CPU time for generating each model.



%---------------------------
% Model specific parameters.
%---------------------------
[ fitfn resfn degenfn psize numpar ] = getModelParam(model_type);


%------------------------
% Check other parameters.
%------------------------
degenmax = 10;  % Max number of inner loop to avoid degeneracy.

%-----------------
% Prepare storage.
%-----------------
n = size(data,2);
par = zeros(numpar,M);
res = zeros(n,M);
inx = zeros(psize,M);
tim = zeros(1,M);

%-----------------
% Random sampling.
%-----------------
fprintf('Random sampling for %.2f seconds...',lim);
t0 = cputime;
for m=1:M    

    % Sample.
    degencnt = 0;
    isdegen = 1;
    while (isdegen==1)&&(degencnt<=degenmax)
        
        % Pick a p-subset.
        [ pinx ] = randsample(n,psize);
        
        %%% Check Sit Points from Same Cluster
        if sum(ismember(pinx,set1idx))==length(pinx)||sum(ismember(pinx,set2idx))==length(pinx)
            continue;
        end         
        
        % Increment degeneracy count.
        degencnt = degencnt + 1;

        
        psub = data(:,pinx);
        
        if isempty(strfind(model_type,'fundamental'))
            % Fit the model, and check for degeneracy.
            isdegen = feval(degenfn,psub);
        else
            % Check for degeneracy.
            [ isdegen F ] = feval(degenfn,psub);
        end
    end
    if (isdegen==1)
        error('Cannot find a valid p-subset!');
    end    
    
    if isempty(strfind(model_type,'fundamental'))
        % Fit the model on the p-subset
        st = feval(fitfn,psub);
        % Compute residuals.
        ds = feval(resfn,st,data);
    else
        % Compute residuals.
        [ ds st ] = feval(resfn,F,data);
    end
    
    
    % Store.
    par(:,m) = st;
    res(:,m) = ds;
    inx(:,m) = pinx;
    tim(1,m) = cputime-t0;
    
    if res(1,m)==0
        a=1;
    end
    
    
    if tim(1,m)>=lim
        par = par(:,1:m);
        res = res(:,1:m);
        inx = inx(:,1:m);
        tim = tim(:,1:m);
        break;
    end
    
end
fprintf('done (%fs)\n',tim(end));

end
