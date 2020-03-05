% Debug stacks that have frames with fluorescent values lower
% than threshold (within batch_zstacktiff2mat)

%% 1) plot max fluorescence per channel per frame
figure();
plot(maxpertime_g, 'g');
hold on;
plot(maxpertime_r, 'r')

%% 2) display slices to replace (with values below threshold)
find(isinf(M(:, 1)))'
find(isinf(M(:, 2)))'

%% 3) plot zstack along Z axes

% 3.1) all planes
pi = [];
pi.lag = 0.01;
pi.sizY = [size(avgim, 1), size(avgim, 2), size(avgim, 3)];
pi.range

slice3Dmatrix(avgim(:, :, :, 1), pi)
slice3Dmatrix(avgim(:, :, :, 2), pi)

% 3.2) selected planes
planes2plot = 275:279;

slice3Dmatrix(avgim(:, :, planes2plot, 1), pi)
slice3Dmatrix(avgim(:, :, planes2plot, 2), pi)

%% 4) replace planes (using info from 2) and 3))

plane2replace = [1];
channel2use = 2;
avgim(:, :, :, channel2use) = ...
    framegapfill(plane2replace, avgim(:, :, :, channel2use));
