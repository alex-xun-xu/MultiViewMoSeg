
clear all;
close all;

addpath(genpath('../../Tools/'));


%% Para
FrameGap = 1; % gap between a pair of frames
max_NumHypoPerFrame = 500;  % Max number of hypotheses sampled from each frame pair

%% Load Seq Information
temp = load('../../Data/SeqList.mat');
SeqList = temp.SeqList;

model_type = lower('homography');   % Model name
[ fitfn resfn degenfn psize numpar ] = getModelParam(model_type);

seq_range = 1:length(SeqList);

for s_i = seq_range
    
    SeqName = SeqList{s_i}; % sequence name
    
    %% Kernel result save path
    save_path = fullfile('../../Results/Kernels/',model_type);
    
    if ~exist(save_path,'dir')
        mkdir(save_path);
    end
    
    kernel_filepath = fullfile(save_path,sprintf('ORK_RandomSamp_Sparse_seq-%s_nhypoframe-%d.mat',...
        SeqName,max_NumHypoPerFrame));

    %%% Load Hypotheses
    save_path = fullfile('../../Results/Hypotheses/',model_type);
    
    hypo_filepath = fullfile(save_path,sprintf('Hypo_RandSamp_Sparse_seq-%s_nHypo-%d.mat',SeqName,max_NumHypoPerFrame));
    temp = load(hypo_filepath,'Hypos');
    
    Model = temp.Hypos;
 
    gt_filepath = fullfile('../../Data/',[SeqName,'_Tracks']);
    temp = load(gt_filepath);
    Data = temp.Data;
    
    %% Compute kernel by accumulating all frame pairs
    K = zeros(Data.nSparsePoints);
    
    for f_i = 1:Data.nFrames-FrameGap
        
        mdl_idx = find(Model.r == f_i)';    % all hypotheses in current frame pair
        Res = [];   % residual w.r.t. hypotheses
        
        for h_i = mdl_idx
            
            r = Model.r(h_i);   % first frame
            v = Model.v(h_i);   % second frame
            
            %% Select points visible on both frames
            visible_pts_ind = Data.visibleSparse(:,f_i) & Data.visibleSparse(:,f_i+1);
            
            y1 = Data.ySparse(:,visible_pts_ind,r);
            y2 = Data.ySparse(:,visible_pts_ind,v);
            
            %% Normalise raw correspondences.
            dat_img_1 = normalise2dpts(y1);
            dat_img_2 = normalise2dpts(y2);
            normalized_data = [ dat_img_1 ; dat_img_2 ];
            
            %% Calculate Residual
            Res = [Res feval(resfn,Model.H(:,h_i),normalized_data)];
            
        end
        
        %% Compute ORK kernel
        [ foo, resinx] = sort(Res,2);
        h = round(0.1*size(Res,2));
        
        K_temp = zeros(Data.nSparsePoints);
        K_ORK = computeIntersection(resinx',resinx',h);
        K_temp(visible_pts_ind,visible_pts_ind)  = K_ORK;
        
        K = K + K_temp;
    end
    
    
    %% Save Results
    save(kernel_filepath,'K');
    
    fprintf('Finish %d-th seq\n',s_i);
    
end


