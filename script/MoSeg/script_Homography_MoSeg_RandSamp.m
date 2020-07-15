%% script to conduct motion segmentation

clear all;
close all;


addpath(genpath('../../Tools/'));


%% Para
FrameGap = 1; % gap between a pair of frames
max_NumHypoPerFrame = 500;  % Max number of hypotheses sampled from each frame pair

AlphaRange = 5:15;  % range of power scaling parameter for eNN sparsification

model_type = lower('homography');

%% Load Seq Information
temp = load('../../Data/SeqList.mat');
SeqList = temp.SeqList;

seq_range = 1:length(SeqList);

%% Evaluate all power scaling parameters
for Alpha = AlphaRange
    
    %%% motion segmentation result save path
    result_path = fullfile('../../Results/MoSeg/',model_type);
    
    if ~exist(result_path,'dir')
        mkdir(result_path);
    end
    
    result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g.mat',...
        max_NumHypoPerFrame,Alpha));

    error = [];
    ClusterIdx = [];
    
    %% Evaluate all sequences
    for s_i = seq_range
        
        SeqName = SeqList{s_i};

        %% Load GroundTruth Data
        gt_filepath = fullfile('../../Data/',[SeqName,'_Tracks']);
        temp = load(gt_filepath);
        Data = temp.Data;
        
        %% Load Kernel Matrix
        save_path = fullfile('../../Results/Kernels/',model_type);
        
        kernel_filepath = fullfile(save_path,sprintf('ORK_RandomSamp_Sparse_seq-%s_nhypoframe-%d.mat',...
            SeqName,max_NumHypoPerFrame));
        
        temp = load(kernel_filepath);
        
        K = temp.K;
        
        %% Normalize affinity matrix by the Cooccurrence of points       
        PtsOcc = double(Data.visibleSparse);    % point occurrence across all frames
        CoocNormalizer = PtsOcc*PtsOcc'+ 0.1;    % points cooccurrence
        
        K = K./CoocNormalizer;
        
        %% Sparsify affinity matrix
        K = func_Adapt_eNN(K,Alpha)+eps;
        
        %% Do spectral clustering
        nMotion = max(Data.GtLabel);
        ClusterIdx{s_i} = SpectralClustering_svd(K,nMotion,'normalized');
        
        if isrow(ClusterIdx{s_i})
            ClusterIdx{s_i} = ClusterIdx{s_i}';
        end
        
        %% Eval Classification Error Rate
        error(s_i) = Misclassification(ClusterIdx{s_i},Data.GtLabel);
        
        fprintf('seq-%d error=%.2f%%\n',s_i,100*error(s_i));

    end
    
    fprintf('\nalpha = %d mean error = %.2f%%\n\n',Alpha,100*mean(error));
    
    %% Save Results
    
    save(result_filepath,'error','ClusterIdx');
    
end

