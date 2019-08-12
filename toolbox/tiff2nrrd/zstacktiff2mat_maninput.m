%% Debug batch_zstacktiff2mat
%avgim(:, :, :, 1) = framegapfill(find(isinf(M(:, 1))), avgim(:, :, :, 1));
%avgim(:, :, :, 2) = framegapfill(find(isinf(M(:, 2))), avgim(:, :, :, 2));
pi = []; pi.sizY = [size(avgim, 1), size(avgim, 2), size(avgim, 3)];
%slice3Dmatrix(avgim(:, :, :, 1), pi)
%slice3Dmatrix(avgim(:, :, :, 2), pi)

%% plot max per channel
figure(); plot(MaxPerTime_g, 'g'); hold on; plot(MaxPerTime_r, 'r')

%% display slices to replace
find(isinf(M(:, 1)))'
find(isinf(M(:, 2)))'

%% plot along Z axes
pi = []; pi.lag = 0.01; pi.sizY = [size(avgim, 1), size(avgim, 2), size(avgim, 3)];
channel = 2;
%slice3Dmatrix(avgim(:, :, :, 1), pi)
slice3Dmatrix(avgim(:, :, :, channel), pi)

%% plot some planes
channel = 1; planes2plot = 275:279;
pi.lag = 1;
pi.sizY = [size(avgim, 1), size(avgim, 2),  numel(planes2plot)];
slice3Dmatrix(avgim(:, :, planes2plot, channel), pi)

%% replace planes
plane2replace = [1:4];
channel = 1;
avgim(:, :, :, channel) = framegapfill(plane2replace, avgim(:, :, :, channel));
