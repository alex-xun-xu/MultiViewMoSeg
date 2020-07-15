%% script to do CoRegularizaiton motion segmentation

clear all;
close all;


addpath(genpath('../../Tools/'));

%% Para
max_NumHypoPerFrame = 500;
FrameGap = 1;
model_type = 'CoReg';

lambda_range = 1e-2;    % range of lambda to evaluate

Alpha_Range = 5:15; % range of power scaling parameter to evaluate

%% Load Seq Information
temp = load('../../Data/SeqList.mat');
SeqList = temp.SeqList;

seq_range = 1:length(SeqList);

for Alpha = Alpha_Range
    for lambda = lambda_range
        
        %% motion segmentation result save path
        result_path = fullfile('../../Results/MoSeg/',model_type);
        
        if ~exist(result_path,'dir')
            mkdir(result_path);
        end
        
        result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g_lambda-%g.mat',...
            max_NumHypoPerFrame,Alpha,lambda));
        
        %% motion segmentation on all sequences
        error = [];
        ClusterIdx = [];
        
        for s_i = seq_range
            
            SeqName = SeqList{s_i};
            
            %% Load GroundTruth Data
            gt_filepath = fullfile('../../Data/',[SeqName,'_Tracks']);
            temp = load(gt_filepath);
            Data = temp.Data;
            
            %% Normalize affinity matrix by the Cooccurrence of points            
            PtsOcc = double(Data.visibleSparse);    % point occurrence across all frames
            CoocNormalizer = PtsOcc*PtsOcc'+ 0.1;    % points cooccurrence
            
            %% Load Affine Kernel Matrix
            save_path = fullfile('../../Results/Kernels/','affine');
            
            kernel_filepath = fullfile(save_path,sprintf('ORK_RandomSamp_Sparse_seq-%s_nhypoframe-%d.mat',...
                SeqName,max_NumHypoPerFrame));
            
            temp = load(kernel_filepath);
            
            K_A = temp.K;            
            K_A = K_A./(CoocNormalizer+0.1);
            [K_A,~] = func_Adapt_eNN(K_A,Alpha);
            
            %% Load H Kernel Matrix
            save_path = fullfile('../../Results/Kernels/','homography');
            
            kernel_filepath = fullfile(save_path,sprintf('ORK_RandomSamp_Sparse_seq-%s_nhypoframe-%d.mat',...
                SeqName,max_NumHypoPerFrame));
            
            temp = load(kernel_filepath);
            
            K_H = temp.K;
            K_H = K_H./(CoocNormalizer+0.1);
            [K_H,~] = func_Adapt_eNN(K_H,Alpha);
            
            %% Load F Kernel Matrix
            save_path = fullfile('../../Results/Kernels/','fundamental');
            
            kernel_filepath = fullfile(save_path,sprintf('ORK_RandomSamp_Sparse_seq-%s_nhypoframe-%d.mat',...
                SeqName,max_NumHypoPerFrame));
            
            temp = load(kernel_filepath);
            
            K_F = temp.K;
            K_F = K_F./(CoocNormalizer+0.1);
            [K_F,~] = func_Adapt_eNN(K_F,Alpha);
            
            %% CoRegularization motion segmentation
                                  
            K = [];
            K(:,:,1) = K_A+eps;
            K(:,:,2) = K_H+eps;
            K(:,:,3) = K_F+eps;
            
            nMotion = max(Data.GtLabel);
            epsilon = 1e-6;
            MaxItr = 15;
            
            %% CoRegularization
            [U_CoReg, ~ , ~] = func_CoRegularize_eig(K,nMotion,lambda,epsilon,MaxItr);
            
            %% Normalize
            U_All = [];
            for k_i = 1:size(U_CoReg,3)
                U_All = [U_All func_L2Normalize(U_CoReg(:,:,k_i),2)];
            end
            
            %%% Normalize
            U = func_L2Normalize(U_All,2);
            
            ClusterIdx{s_i} = kmeans(U, nMotion, 'replicates',500, 'start', 'cluster', ...
                'EmptyAction', 'singleton');
            
            if isrow(ClusterIdx{s_i})
                ClusterIdx{s_i} = ClusterIdx{s_i}';
            end
            
            %%% Evaluate Classification Error
            error(s_i) = Misclassification(ClusterIdx{s_i},Data.GtLabel);
                        
            fprintf('Sequence %s Error = %.2f%% \n',SeqName,100*error(s_i));
            
        end
        
        fprintf('Overall Miss Classification Rate = %.2f\n',100*mean(error));
        
        %% Save Results
        
        save(result_filepath,'ClusterIdx','error');
        
    end
    
end

