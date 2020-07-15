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


clear all;
close all;

%----------
% Set path.
%----------
addpath('./sampling_methods');
addpath('./model_specific');

%----------------------
% Compile required dll.
%----------------------
if (exist('computeIntersection')~=3)
    mex computeIntersection.c
end


%------------------
% Select test data.
%------------------
type = input('Type of model? [1 Homography] [2 Fundamental matrix] > ');


%----------
% Get data.
%----------
switch type
    case 1
        model_type = 'homography';
        load data/wadham.mat;
    case 2
        model_type = 'fundamental';
        load data/dinobooks.mat;
end
    
%-------------------------------
% Normalise raw correspondences.
%-------------------------------
dat_img_1 = normalise2dpts(data(1:3,:));
dat_img_2 = normalise2dpts(data(4:6,:));
normalized_data = [ dat_img_1 ; dat_img_2 ];


%-----------------------------------------
% Sample hypotheses and compute residuals.
%-----------------------------------------

% Maximum number of hypotheses.
M = 25000;

% Maximum CPU seconds allowed
lim = 10;

% Number of structures.
numfun = max(label);

% Storage.
par = cell(2,1);
res = cell(2,1);
inx = cell(2,1);
tim = cell(2,1);
hit = cell(2,4);
met = char('Random','Multi-GS');

% Random sampling.
[ par{1} res{1} inx{1} tim{1} ] = randomSampling(lim,normalized_data,M,model_type);

% Guided-sampling using the Multi-GS method.
[ par{2} res{2} inx{2} tim{2} ] = multigsSampling(lim,normalized_data,M,10,model_type);

% Evaluate.
sets = cell(numfun,1);
for i=1:numfun
    sets{i} = find(label==i);
end


inlier_inx = cell(2,1);
inlier_inx{1} = zeros(size(inx{1},1), M);
inlier_inx{2} = zeros(size(inx{2},1), M);

for meth=1:2
    allhit = false;
    hit{meth,1} = size(inx{meth},2);
    hit{meth,2} = zeros(numfun,1);
    for m=1:size(inx{meth},2)
        for i=1:numfun
            if sum(ismember(inx{meth}(:,m),sets{i}))==size(inx{meth}(:,m),1)
                hit{meth,2}(i) = hit{meth,2}(i) + 1;
                inlier_inx{meth}(:,m) = inx{meth}(:,m);
            end
        end
        if (allhit==false)&&(sum(hit{meth,2}>0)==numfun)
            allhit = true;
            hit{meth,3} = m;
            hit{meth,4} = tim{meth}(m);
        end
    end
    if (allhit==false)
        hit{meth,3} = inf;
        hit{meth,4} = inf;
    end
end        
for meth=1:2
    fprintf('%s - total samples = %d - number of all-inlier samples hit = %d - hit all structures in %d steps (%f s).\n', ...
        met(meth,:), ...
        hit{meth,1}, ...
        sum(hit{meth,2}), ...
        hit{meth,3}, ...
        hit{meth,4});
end


%--------------
% Plot results.
%--------------
fmt = {'y+','bs','ro','g<','cd'};

% Plot.
figure;
axes('Position', [0 0 1 1]) ;
subplot(3,1,1);
imshow([ img1 img2 ]);
hold on;

for i = min(label):max(label)
    inx = (label == i);
    plot(data(1,inx), data(2,inx), fmt{i+1},'MarkerSize',10);
    plot(data(4,inx)+size(img1,2), data(5,inx),fmt{i+1},'MarkerSize',10);
end

title('Keypoints (outliers are marked as +)');


for meth=1:2
    % Indices of all-inlier p-subsets.
    inx = inlier_inx{meth}(:);
    inx = inx(inx>0);

    subplot(3,1,meth+1);
    imshow([ img1 img2 ]);
    hold on;
        
    title(sprintf('All-inlier samples given by %s',met(meth,:)));
        
    if isempty(inx)
        continue;
    end
    
    for i = unique(label(inx))
        each_cluster = inx(label(inx)==i);
        plot(data(1,each_cluster), data(2,each_cluster), fmt{i+1},'MarkerSize',10);
        plot(data(4,each_cluster)+size(img1,2), data(5,each_cluster),fmt{i+1},'MarkerSize',10);
    end
    
end

