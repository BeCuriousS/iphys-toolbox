%% reset legend
lg = {};

%% plot preprocessed BVP signal
idx = 4;
plot_curve = 6;

switch plot_curve
    case 1
        ts1 = timeseries(dp_data(idx).BVP_pp, dp_data(idx).TS_pp);
        lg{end+1} = 'deepPerfusion BVP pp';
        plot(ts1, 'r');
        hold on
        grid on
    case 2
        ts2 = timeseries(pos_data(idx).BVP_pp, pos_data(idx).TS_pp);
        lg{end+1} = 'POS BVP pp';
        plot(ts2, 'g');
        hold on
        grid on
    case 3
        ts3 = timeseries(gt_data(idx).BVP_pp, gt_data(idx).TS_pp);
        lg{end+1} = 'Ground truth BVP pp';
        plot(ts3, 'b');
        hold on
        grid on
    case 4
        ts4 = timeseries(gt_data(idx).BVP, gt_data(idx).TS);
        lg{end+1} = 'Ground truth BVP';
        plot(ts4, 'm');
        hold on
        grid on
    case 5
        ts5 = timeseries((pos_data(idx).BVP_pp - mean(pos_data(idx).BVP_pp))/std(pos_data(idx).BVP_pp), pos_data(idx).TS_pp);
        lg{end+1} = 'POS BVP pp standardized';
        plot(ts5, 'r');
        hold on
        grid on
    case 6
        ts6 = timeseries((dp_data(idx).BVP_pp - mean(dp_data(idx).BVP_pp))/std(dp_data(idx).BVP_pp), dp_data(idx).TS_pp);
        lg{end+1} = 'deepPerfusion BVP pp standardized';
        plot(ts6, 'b');
        hold on
        grid on
end
legend(lg)
% findpeaks(gt_data(i).BVP_pp, 'MinPeakDistance', ceil(resample_fs * 1/3), 'MinPeakProminence', 1)

%% check pos roi data
% detector = 'UpperBody';
detector = 'FrontalFaceCart';
% load subject check data
check_data = struct();
for i = 1:size(gt_data, 2)
    data = load(['D:\GitHub\iphys-toolbox\2018_12_UBFC_Dataset\bvp_pos_' detector '\' gt_data(i).subj '_check.mat']);
    idx_first = 1;
    if isempty(data.check_data(1).frameROI)
        idx_first = 2;
    end
    check_data(i).data = data.check_data(idx_first:end);
end
% 28 leftarrow
% 29 rightarrow
% 30 uparrow
% 31 downarrow
i1 = 1;
i2 = 1;
i1max = size(check_data, 2);
i2max = size(check_data(i1).data, 2);
figure
while true
    k = waitforbuttonpress;
    k = double(get(gcf,'CurrentCharacter'));
    if k == 30
        i1 = i1 + 1;
        if i1 > i1max
            i1 = i1max;
        end
        i2max = size(check_data(i1).data, 2);
    elseif k == 31
        i1 = i1 - 1;
        if i1 == 0
            i1 = 1;
        end
        i2max = size(check_data(i1).data, 2);
    elseif k == 29
        i2 = i2 + 1;
        if i2 > i2max
            i2 = i2max;
        end
    elseif k == 28
        i2 = i2 - 1;
        if i2 == 0
            i2 = 1;
        end
    end
    fprintf('Idx %i - Subject %s - Frame %i \n', i1, gt_data(i1).subj, check_data(i1).data(i2).frameIndex);
    imshow(check_data(i1).data(i2).skinROI);
end
