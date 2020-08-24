% Simple code to read ground truth of the UBFC_DATASET
% If you use the dataset, please cite:
% 
% S. Bobbia, R. Macwan, Y. Benezeth, A. Mansouri, J. Dubois, 
% Unsupervised skin tissue segmentation for remote photoplethysmography, 
% Pattern Recognition Letters, Elsevier, 2017.
% 
% yannick.benezeth@u-bourgogne.fr

clear all;
close all;
clc;

% dataset folder
root        =   '~/cbppg/2018_12_UBFC_Dataset/measurements';

% get folder list
dirs     = {
    'subject36',
    'subject15',
    'subject37',
    'subject42',
    'subject41',
    'subject35',
    'subject11',
    'subject16',
    'subject9',
    'subject34',
    'subject13',
    'subject45',
    'subject49',
    'subject17',
    'subject30',
    'subject5',
    'subject48',
    'subject24',
    'subject10',
    'subject8',
    'subject18'
 };

%%
for i=1:size(dirs)
    fprintf('Processing subject %s \n', dirs{i});
    
    vidFolder   =   [root '/' dirs{i}];    
    
    % roiDetector to use - UpperBody or FrontalFaceCART
    roiDetAlg = 'UpperBody';
%     roiDetAlg = 'FrontalFaceCART';
    
    % load ground truth
    ground_truth = dlmread( [vidFolder '/ground_truth.txt' ] );
    gt_trace = ground_truth( 1, : );
    gt_HR = ground_truth( 2, : );
    gt_time = ground_truth( 3, : );
    
    % save computations
    [BVP, check_data] = POS_WANG_BVP([vidFolder '/vid.avi'], 30, roiDetAlg);
    path = ['2018_12_UBFC_Dataset/bvp_pos_' roiDetAlg '/'];
    mkdir(path);
    save([path dirs{i} '.mat'], 'BVP');
    save([path dirs{i} '_check' '.mat'], 'check_data');
end
