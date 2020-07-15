%% Sample affine hypotheses

clear all;
close all;


%% Set path.
addpath(genpath('../../Tools/'));


%% Para
FrameGap = 1;   % gap between a pair of frames
max_NumHypoPerFrame = 500;  % Max number of hypotheses sampled from each frame pair

%% Load Seq Information
temp = load('../../Data/SeqList.mat');
SeqList = temp.SeqList;

model_type = lower('fundamental');

seq_range = 1:length(SeqList);

%% Sample hypotheses from all sequences
for s_i = seq_range
    
    SeqName = SeqList{s_i}; % sequence name
    
    %%% Load Ground-Truth Data
    gt_filepath = fullfile('../../Data/',[SeqName,'_Tracks.mat']);
    temp = load(gt_filepath);
    DataOrg = temp.Data;
    
    gt_filepath = fullfile('../../Data/',[SeqName,'_Tracks.mat']);
    temp = load(gt_filepath);
    Data = temp.Data;
    
    num_frames = Data.nFrames;
    
    %%% Save Path for hypotheses
    save_path = fullfile('../../Results/Hypotheses/',model_type);
    
    if ~exist(save_path,'dir')
        mkdir(save_path);
    end
    
    hypo_filepath = fullfile(save_path,sprintf('Hypo_RandSamp_Sparse_seq-%s_nHypo-%d.mat',SeqName,max_NumHypoPerFrame));

    %%% Initialize All Hypotheses
    Hypos.H = [];
    Hypos.r = [];
    Hypos.v = [];
    Hypos.supp = [];
    
    
    for f_i = 1:num_frames-FrameGap
        
        %% Prepare candidate data
        r = f_i;    % first frame
        v = r+FrameGap; % second frame
        
        %%% Select points visible on both frames        
        visible_pts_ind = Data.visibleSparse(:,f_i) & Data.visibleSparse(:,f_i+1);
        
        y1 = Data.ySparse(:,visible_pts_ind,r);
        y2 = Data.ySparse(:,visible_pts_ind,v);

        %% Normalise raw correspondences.
        dat_img_1 = normalise2dpts(y1);
        dat_img_2 = normalise2dpts(y2);
        normalized_data = [ dat_img_1 ; dat_img_2 ];
        
        % Maximum CPU seconds allowed
        lim = 20;
        
        % Storage.
        par = cell(2,1);
        res = cell(2,1);
        inx = cell(2,1);
        tim = cell(2,1);
        hit = cell(2,4);
        met = char('Random','Multi-GS');
        
        % Random sampling.
        [ par{1} res{1} inx{1} tim{1} ] = randomSampling(lim,normalized_data,max_NumHypoPerFrame,model_type);
        
        % Guided-sampling using the Multi-GS method. (alternative sampling strategy)
%         [ par{2} res{2} inx{2} tim{2} ] = multigsSampling(lim,normalized_data,max_NumHypoPerFrame,10,model_type);
        
        
        %% Accumulate Hypotheses
        Hypos.H = [Hypos.H  par{1}];
        Hypos.r = [Hypos.r ; r*ones(size(par{1},2),1)];
        Hypos.v = [Hypos.v ; v*ones(size(par{1},2),1)];
        Hypos.supp = [Hypos.supp  inx{1}];

        
    end
    
    %% Save Hypotheses
    save(hypo_filepath,'Hypos');
    
    fprintf('Finish %d-th seq\n',s_i);
    
end


