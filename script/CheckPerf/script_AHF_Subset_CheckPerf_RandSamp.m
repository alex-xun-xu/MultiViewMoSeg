%% Check the subset constrained motion segmentation results


clear all;
close all;


max_NumHypoPerFrame = 500;
FrameGap = 1;
model_type = 'Subset';

gamma_range = [1e-2];
Alpha_Range = 5:15;


%% Load Seq Information
temp = load('../../Data/SeqList.mat');
SeqList = temp.SeqList;

seq_range = 1:length(SeqList);

error_avg = [];
error_med = [];
valid_gamma = [];

for  alf_i = 1:length(Alpha_Range)
    
    Alpha = Alpha_Range(alf_i);
    
    for ld_i = 1:length(gamma_range)
        
        gamma = gamma_range(ld_i);
        
        Para = [];
        Para.MaxItr = 8;   % max iteration
        %     Para.gamma = gamma;
        Para.gamma = gamma;
        
        result_path = fullfile('../../Results/MoSeg/',model_type);
        
        result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g_gamma-%g.mat',...
            max_NumHypoPerFrame,Alpha,gamma));
        
        if ~exist(result_filepath,'file')
            continue;
        end        
        
        temp = load(result_filepath);
        
        error = temp.error;
        
        
        fprintf('Overall Miss Classification Rate = %.2f\n',100*mean(error(end,:)));

        error_avg(alf_i,ld_i) = mean(error);
        error_med(alf_i,ld_i) = median(error);

    end
    valid_gamma = [valid_gamma gamma];
    
end

%% Display all view results
%% Concatenate view
figure;
subplot(1,2,1);
imagesc(error_avg);
colorbar;
title('Subset Mean Error');

xlabel('gamma');
ylabel('alpha');

%%% Median view
subplot(1,2,2);
imagesc(error_med);
colorbar;
title('Subset Median Error');

