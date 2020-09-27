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
root        =   '~/Desktop/tmp';

% get folder list
dirs     = {
    '00106'
 };

%%
for i=1:size(dirs)
    fprintf('Processing subject %s \n', dirs{i});
    
    vidFolder   =   [root '/' dirs{i}];    
    
    % roiDetector to use - UpperBody or FrontalFaceCART
%     roiDetAlg = 'UpperBody';
    roiDetAlg = 'FrontalFaceCART';
%     roiDetAlg = 'NoDet';

    useFGTransform = false;
    
    method = 'pos';
    
    % load ground truth
%     ground_truth = dlmread( [vidFolder '/ground_truth.txt' ] );
%     gt_trace = ground_truth( 1, : );
%     gt_HR = ground_truth( 2, : );
%     gt_time = ground_truth( 3, : );
    
    % save computations
    if strcmp(method, 'pos')
        [BVP, check_data] = POS_WANG_BVP([vidFolder '/video.avi'], 30, roiDetAlg, useFGTransform);
    elseif strcmp(method, 'chrom')
        [BVP, check_data] = CHROM_DEHAAN_BVP([vidFolder '/vid.avi'], 30, roiDetAlg, useFGTransform);
    end
    
    
    path = ['2018_12_UBFC_Dataset/bvp_' method '_' roiDetAlg '/'];
    if strcmp(method, 'pos') && useFGTransform
        path = ['2018_12_UBFC_Dataset/bvp_' method '_FGTransform_' roiDetAlg '/'];
    end
    mkdir(path);
    save([path dirs{i} '.mat'], 'BVP');
    save([path dirs{i} '_check' '.mat'], 'check_data');
end
