%% check the co-regularized motion segmentation results

clear all;
close all;


max_NumHypoPerFrame = 500;
FrameGap = 1;
model_type = 'CoReg';

lambda_range = [1e-2];
Alpha_Range = 5:15;


%% Load Seq Information
temp = load('../../Data/SeqList.mat');
SeqList = temp.SeqList;

seq_range = 1:length(SeqList);

error_avg = [];
error_med = [];
valid_lambda = [];

for  alf_i = 1:length(Alpha_Range)
    
    Alpha = Alpha_Range(alf_i);
    
    for ld_i = 1:length(lambda_range)
        
        lambda = lambda_range(ld_i);
        
        Para = [];
        Para.MaxItr = 15;   % max iteration for CoReg
        %     Para.lambda = lambda;
        Para.lambda = lambda;
        
        result_path = fullfile('../../Results/MoSeg/',model_type);
        
        result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g_lambda-%g.mat',...
            max_NumHypoPerFrame,Alpha,lambda));
        
        if ~exist(result_filepath,'file')
            continue;
        end        
        
        temp = load(result_filepath);
        
        error = temp.error;
        
        
        fprintf('Overall Miss Classification Rate = %.2f\n',100*mean(error(end,:)));

        error_avg(alf_i,ld_i) = mean(error);
        error_med(alf_i,ld_i) = median(error);

    end
    valid_lambda = [valid_lambda lambda];
    
end

%% Display all view results
%% Concatenate view
figure;
subplot(1,2,1);
imagesc(error_avg);
colorbar;
title('CoRegularization Mean Error');

xlabel('lambda');
ylabel('alpha');

%%% Median view
subplot(1,2,2);
imagesc(error_med);
colorbar;
title('CoRegularization Median Error');

