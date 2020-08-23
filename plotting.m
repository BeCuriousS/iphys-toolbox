%% plot preprocessed BVP signal
idx = 8;

ts1 = timeseries(dp_data(idx).BVP_pp, dp_data(idx).TS_pp);
ts2 = timeseries(pos_data(idx).BVP_pp, pos_data(idx).TS_pp);
ts3 = timeseries(gt_data(idx).BVP_pp, gt_data(idx).TS_pp);
ts4 = timeseries(gt_data(idx).BVP, gt_data(idx).TS);

plot(ts1, 'r');
hold on
grid on
plot(ts2, 'g');
plot(ts3, 'b');
plot(ts4, 'm');
legend({'deepPerfusion BVP pp', 'POS BVP pp', 'Ground truth BVP pp', 'Ground truth BVP'});

% findpeaks(gt_data(i).BVP_pp, 'MinPeakDistance', ceil(resample_fs * 1/3), 'MinPeakProminence', 1)