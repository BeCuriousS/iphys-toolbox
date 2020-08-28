clearvars
%% set paths to BVP signals
% DeepPerfusion BVP
dp_path = 'D:\GitLab\deepperfusion_ippg\data\eval\20200821-061926_2\2018_12_UBFC_Dataset\measurements\';
% dp_path = 'D:\GitLab\deepperfusion_ippg\data\eval\20200820-070753\2018_12_UBFC_Dataset\measurements\';
% POS BVP
pos_path = 'D:\GitHub\iphys-toolbox\2018_12_UBFC_Dataset\bvp_pos_FrontalFaceCART\';
% pos_path = 'D:\GitHub\iphys-toolbox\2018_12_UBFC_Dataset\bvp_pos_UpperBody\';
% pos_path = 'D:\GitHub\iphys-toolbox\2018_12_UBFC_Dataset\bvp_pos_FGTransformFrontalFaceCART\';
% pos_path = 'D:\GitHub\iphys-toolbox\2018_12_UBFC_Dataset\bvp_pos_NoDet\';
% ground truth data
gt_path = 'D:\GitHub\iphys-toolbox\2018_12_UBFC_Dataset\ground_truth\';

% read data into memory
dp_folders = dir(dp_path);
dp_data = struct();
for i = 3:length(dp_folders)
    data = load(fullfile(dp_folders(i).folder, dp_folders(i).name, 'ippg_deepPerfusion.mat'));
    % filter values appearing twice (because of overlapping windows
    % convert to double as the network output is single precision
    [~, idx_unq, ~] = unique(data.ts_seq);
    dp_data(i-2).BVP_orig = double(data.ippg_seq.');
    dp_data(i-2).BVP = dp_data(i-2).BVP_orig(idx_unq);
    dp_data(i-2).TS_orig = double(data.ts_seq*1e-6);
    dp_data(i-2).TS = dp_data(i-2).TS_orig(idx_unq);
    dp_data(i-2).subj = dp_folders(i).name;
end

gt_folders = dir(gt_path);
gt_data = struct();
cnt = 0;
for i = 3:length(gt_folders)
    is_in = false;
    for j = 1:size(dp_data, 2)
        if strcmp(gt_folders(i).name, dp_data(j).subj)
            is_in = true;
            idx_dp = j;
            cnt = cnt + 1;
            break
        end
    end
    if is_in
        data = dlmread(fullfile(gt_folders(i).folder, gt_folders(i).name, 'ground_truth.txt'));

%         mask = (data(3, :) >= dp_data(idx).TS(1)) & (data(3, :) <= dp_data(idx).TS(end));
%         mask = ismember(data(3,:), dp_data(idx).TS); % not working cause
%         of small rounding errors (?)

        mask = zeros(1, length(data(3,:)));
        for j = 1:length(dp_data(idx_dp).TS)
            idx = find(abs((data(3, :) - dp_data(idx_dp).TS(j))) < 0.0001);
            if ~isempty(idx)
                mask(1, idx) = true;
            end
        end
        mask = logical(mask);
        [~, idx_unq, ~] = unique(data(3, :));
        mask_unq = zeros(size(mask));
        mask_unq(idx_unq) = 1;
        mask = mask & logical(mask_unq);
        
        gt_data(cnt).BVP_orig = data(1, :);
        gt_data(cnt).BVP = data(1, mask);
        gt_data(cnt).TS_orig = data(3, :);
        gt_data(cnt).TS = data(3, mask);
        gt_data(cnt).HR_orig = data(2, :);
        gt_data(cnt).HR = data(2, mask);
        gt_data(cnt).MASK = mask;
        subj = split(gt_folders(i).name, '.');
        gt_data(cnt).subj = subj{1};
    end
end

pos_files = dir(pos_path);
pos_data = struct();
cnt = 0;
for i = 3:length(pos_files)
    subj = split(pos_files(i).name, '.');
    subj = subj{1};
    for j = 1:size(gt_data, 2)
        if strcmp(subj, gt_data(j).subj)
            idx = j;
        end
    end
    if length(split(subj, '_check')) == 1
        cnt = cnt + 1;
        data = load(fullfile(pos_files(i).folder, pos_files(i).name));
        pos_data(cnt).BVP_orig = data.BVP;
        pos_data(cnt).BVP = data.BVP(gt_data(idx).MASK);
        pos_data(cnt).TS_orig = gt_data(idx).TS_orig;
        pos_data(cnt).TS = gt_data(idx).TS;
        pos_data(cnt).subj = subj;
    end
end

%% calculate metrics
% set parameters
spectrum_window_size = 10; % in sec
window_stride = 1; % in sec
ref_HR = 'PR'; % 'PR' (from BVP), 'PR_fromHR', 'PR_fromPeaks'
resample_fs = 125; % Hz

% Parameters of the bandpass filter
lower_cutOff = 0.6;
upper_cutOff = 8.2;
filter_order = 2;
interp_method = 'pchip';

% Parameters for SNR calculation
snr_freq_bounds = [0.5 8.2]; % Hz
snr_freq_bounds_gt = [0.5 8.2]; % Hz

% resample BVP to specified sample frequency and apply the bandpass-filter
d = designfilt('bandpassiir','FilterOrder',filter_order, ...
    'HalfPowerFrequency1',lower_cutOff,'HalfPowerFrequency2',upper_cutOff, ...
    'SampleRate',resample_fs);
for i = 1:size(gt_data, 2)
    x = gt_data(i).TS;
    xq = x(1):1/resample_fs:x(end);
    
    dp_data(i).BVP_pp = filtfilt(d, interp1(x, dp_data(i).BVP, xq, interp_method));
    pos_data(i).BVP_pp = filtfilt(d, interp1(x, pos_data(i).BVP, xq, interp_method));
    gt_data(i).BVP_pp = filtfilt(d, interp1(x, gt_data(i).BVP, xq, interp_method));
    
    gt_data(i).HR_pp = interp1(x, gt_data(i).HR, xq, interp_method);
    
    dp_data(i).TS_pp = xq;
    gt_data(i).TS_pp = xq;
    pos_data(i).TS_pp = xq;
end

% compute windowed pulse rates from BVP
for i = 1:size(gt_data,2)
    % compute sample freq
    duration = (gt_data(i).TS_pp(end) - gt_data(i).TS_pp(1));
    FS = resample_fs;
    
    % iterate over windows with chosen stride
    n_windows = floor(duration - spectrum_window_size);
    t = 0;
    for j = 1:n_windows
        mask = (gt_data(i).TS_pp >= t) & (gt_data(i).TS_pp < t + spectrum_window_size);
        
        dp_data(i).PR(j) = prpsd(dp_data(i).BVP_pp(mask), FS, 40, 240, false);
        pos_data(i).PR(j) = prpsd(pos_data(i).BVP_pp(mask), FS, 40, 240, false);
        gt_data(i).PR(j) = prpsd(gt_data(i).BVP_pp(mask), FS, 40, 240, false);
        gt_data(i).PR_fromHR(j) = mean(gt_data(i).HR_pp(mask));
        
        [pks, locs] = findpeaks(gt_data(i).BVP_pp(mask), 'MinPeakDistance', ceil(resample_fs * 1/3), 'MinPeakProminence', 1);
        
        gt_data(i).PR_fromPeaks(j) = 1/median(diff(gt_data(i).TS_pp(locs))) * 60;
        
        t = t + window_stride;
    end
end

% compute windowed SNR
for i = 1:size(gt_data,2)
    % compute sample freq
    duration = (gt_data(i).TS_pp(end) - gt_data(i).TS_pp(1));
    FS = resample_fs;
    
    % iterate over windows with chosen stride
    n_windows = floor(duration - spectrum_window_size);
    t = 0;
    for j = 1:n_windows
        mask = (gt_data(i).TS_pp >= t) & (gt_data(i).TS_pp < t + spectrum_window_size);
        
        dp_data(i).SNR(j) = bvpsnr(dp_data(i).BVP_pp(mask), FS, gt_data(i).(ref_HR)(j), false, snr_freq_bounds);
        pos_data(i).SNR(j) = bvpsnr(pos_data(i).BVP_pp(mask), FS, gt_data(i).(ref_HR)(j), false, snr_freq_bounds);
        try
            gt_data(i).SNR(j) = bvpsnr(gt_data(i).BVP_pp(mask), FS, gt_data(i).(ref_HR)(j), false, snr_freq_bounds_gt);
        catch
            gt_data(i).SNR(j) = nan;
        end
        
        t = t + window_stride;
    end
end

% compute metrics
mae_dp = zeros(1, size(gt_data,2));
mae_pos = zeros(1, size(gt_data,2));
rmse_dp = zeros(1, size(gt_data,2));
rmse_pos = zeros(1, size(gt_data,2));
snr_dp = zeros(1, size(gt_data,2));
snr_pos = zeros(1, size(gt_data,2));

mae_dp_signal = zeros(1, size(gt_data,2));
rmse_dp_signal = zeros(1, size(gt_data,2));

pearson_r_dp = zeros(1, size(gt_data,2));
pearson_r_pos = zeros(1, size(gt_data,2));
for i = 1:size(gt_data,2)
    mae_dp(i) = mean(abs(dp_data(i).PR - gt_data(i).(ref_HR)));
    mae_pos(i) = mean(abs(pos_data(i).PR - gt_data(i).(ref_HR)));
    
    rmse_dp(i) = mean((dp_data(i).PR - gt_data(i).(ref_HR)).^2)^0.5;
    rmse_pos(i) = mean((pos_data(i).PR - gt_data(i).(ref_HR)).^2)^0.5;
    
    snr_dp(i) = mean(dp_data(i).SNR);
    snr_pos(i) = mean(pos_data(i).SNR);
    
    gt_BVP_scaled = (gt_data(i).BVP_pp - mean(gt_data(i).BVP_pp))/std(gt_data(i).BVP_pp);
    mae_dp_signal(i) = mean(abs(dp_data(i).BVP_pp - gt_BVP_scaled));
    rmse_dp_signal(i) = mean((dp_data(i).BVP_pp - gt_BVP_scaled).^2)^0.5;
    
    ccfmat = corrcoef([dp_data(i).PR.' gt_data(i).(ref_HR).']);
    pearson_r_dp(i) = ccfmat(1,2);
    ccfmat = corrcoef([pos_data(i).PR.' gt_data(i).(ref_HR).']);
    pearson_r_pos(i) = ccfmat(1,2);
end

fprintf('\n>>>>>>>>>    ERRORS OF DEEPPERFUSION BVP SIGNAL TO Ground truth\n');
fprintf('Mean absolute error:      %2.3f \n', mean(mae_dp_signal));
fprintf('Root Mean squared error:  %2.3f \n', mean(rmse_dp_signal));

fprintf('\n>>>>>>>>>    MEAN ABSOLUTE ERRORS TO Ground truth\n');
fprintf('Mean absolute error deepPerfusion: %2.2f BPM \n', mean(mae_dp));
fprintf('Mean absolute error pos:           %2.2f BPM \n', mean(mae_pos));

fprintf('\n>>>>>>>>>    ROOT MEAN SQUARED ERRORS TO Ground truth\n');
fprintf('Root mean squared error deepPerfusion: %2.2f BPM \n', mean(rmse_dp));
fprintf('Root mean squared error pos:           %2.2f BPM \n', mean(rmse_pos));

fprintf('\n>>>>>>>>>    PEARSON CORRELATION TO Ground truth\n');
fprintf('r deepPerfusion: %1.3f BPM \n', mean(pearson_r_dp));
fprintf('r pos:           %1.3f BPM \n', mean(pearson_r_pos));

fprintf('\n>>>>>>>>>    SNR\n');
fprintf('SNR deepPerfusion: %2.2f dB \n', mean(snr_dp));
fprintf('SNR pos:           %2.2f dB \n', mean(snr_pos));

%% compute metrics with confidence intervals
% compute confidence intervals
pci = paramci(fitdist(mae_dp_signal.', 'normal'));
mae_dp_sig_ci = pci(:, 1);
pci = paramci(fitdist(rmse_dp_signal.', 'normal'));
rmse_dp_sig_ci = pci(:, 1);

pci = paramci(fitdist(mae_dp.', 'normal'));
mae_dp_ci = pci(:, 1);
pci = paramci(fitdist(mae_pos.', 'normal'));
mae_pos_ci = pci(:, 1);

pci = paramci(fitdist(rmse_dp.', 'normal'));
rmse_dp_ci = pci(:, 1);
pci = paramci(fitdist(rmse_pos.', 'normal'));
rmse_pos_ci = pci(:, 1);

pci = paramci(fitdist(pearson_r_dp.', 'normal'));
pearson_r_dp_ci = pci(:, 1);
pci = paramci(fitdist(pearson_r_pos.', 'normal'));
pearson_r_pos_ci = pci(:, 1);

pci = paramci(fitdist(snr_dp.', 'normal'));
snr_dp_ci = pci(:, 1);
pci = paramci(fitdist(snr_pos.', 'normal'));
snr_pos_ci = pci(:, 1);

fprintf('\n>>>>>>>>>    ERRORS OF DEEPPERFUSION BVP SIGNAL TO Ground truth\n');
fprintf('Mean absolute error:      %2.3f (%2.2f/%2.2f) (STD: %2.2f)\n', mean(mae_dp_signal), mae_dp_sig_ci(1), mae_dp_sig_ci(2), std(mae_dp_signal));
fprintf('Root mean squared error:  %2.3f (%2.2f/%2.2f) (STD: %2.2f)\n', mean(rmse_dp_signal), rmse_dp_sig_ci(1), rmse_dp_sig_ci(2), std(rmse_dp_signal));

fprintf('\n>>>>>>>>>    MEAN ABSOLUTE ERRORS TO Ground truth\n');
fprintf('Mean absolute error deepPerfusion: %2.2f BPM (%2.2f/%2.2f) (STD: %2.2f)\n', mean(mae_dp), mae_dp_ci(1), mae_dp_ci(2), std(mae_dp));
fprintf('Mean absolute error pos:           %2.2f BPM (%2.2f/%2.2f) (STD: %2.2f)\n', mean(mae_pos), mae_pos_ci(1), mae_pos_ci(2), std(mae_pos));

fprintf('\n>>>>>>>>>    ROOT MEAN SQUARED ERRORS TO Ground truth\n');
fprintf('Root mean squared error deepPerfusion: %2.2f BPM (%2.2f/%2.2f) (STD: %2.2f)\n', mean(rmse_dp), rmse_dp_ci(1), rmse_dp_ci(2), std(rmse_dp));
fprintf('Root mean squared error pos:           %2.2f BPM (%2.2f/%2.2f) (STD: %2.2f)\n', mean(rmse_pos), rmse_pos_ci(1), rmse_pos_ci(2), std(rmse_pos));

fprintf('\n>>>>>>>>>    PEARSON CORRELATION TO Ground truth\n');
fprintf('r deepPerfusion: %1.3f BPM (%2.2f/%2.2f) (STD: %2.2f)\n', mean(pearson_r_dp), pearson_r_dp_ci(1), pearson_r_dp_ci(2), std(pearson_r_dp));
fprintf('r pos:           %1.3f BPM (%2.2f/%2.2f) (STD: %2.2f)\n', mean(pearson_r_pos), pearson_r_pos_ci(1), pearson_r_pos_ci(2), std(pearson_r_pos));

fprintf('\n>>>>>>>>>    SNR\n');
fprintf('SNR deepPerfusion: %2.2f dB (%2.2f/%2.2f) (STD: %2.2f)\n', mean(snr_dp), snr_dp_ci(1), snr_dp_ci(2), std(snr_dp));
fprintf('SNR pos:           %2.2f dB (%2.2f/%2.2f) (STD: %2.2f)\n', mean(snr_pos), snr_pos_ci(1), snr_pos_ci(2), std(snr_pos));











