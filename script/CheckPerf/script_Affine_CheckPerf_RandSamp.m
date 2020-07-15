%% script to check motion segmentation performance

clear; close all;

max_NumHypoPerFrame = 500;
FrameGap = 1;
model_type = lower('affine');

AlphaRange = 1:50;

%% Load Seq Information
temp = load('../../Data/SeqList.mat');
SeqList = temp.SeqList;

seq_range = 1:length(SeqList);

error = [];
valid_alpha = [];

for Alpha = AlphaRange
    
    %% Load Segmentation Result
    result_path = fullfile('../../Results/MoSeg/',model_type);

    result_filepath = fullfile(result_path,sprintf('Error_RandSamp_nhpf-%d_alpha-%g.mat',...
        max_NumHypoPerFrame,Alpha));
    
    if ~exist(result_filepath,'file')
        continue;
    else
        temp = load(result_filepath);
    end
        
    error = [error ; temp.error];
    
    valid_alpha = [valid_alpha Alpha];
     
end

avg_error = mean(error,2);
med_error = median(error,2);

figure;
colors = colormap(hsv(3));
plot(valid_alpha,100*avg_error,'color',colors(1,:),'marker','x'); hold on;
plot(valid_alpha,100*med_error,'color',colors(2,:),'marker','x'); hold on;

xlabel('alpha');
ylabel('error (%)');

grid on;

legend({'Mean Error','Median Error'});

title('Affine Motion Segmentation Performance');
