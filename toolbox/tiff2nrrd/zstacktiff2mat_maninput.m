% Debug stacks that have frames with fluorescent values lower
% than threshold (within batch_zstacktiff2mat)
%
% Usage: batch_zstacktiff2mat stops with an error when avbrain finds Inf
%   (below threshold) or saturated planes. To fix them, rerun it with
%   "dbstop if error" so the workspace of avbrain is available, then run
%   this script section by section (it uses the variables maxpertime_r,
%   maxpertime_g, M and avgim of avbrain) and return the edited avgim.
%   maxpertime_g exists only if two channels were collected.

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
% (parameters for slice3Dmatrix; do not call this variable 'pi', that
%   would hide the builtin constant)
plot_pars = [];
plot_pars.lag = 0.01;
plot_pars.sizY = [size(avgim, 1), size(avgim, 2), size(avgim, 3)];
plot_pars.range = [0 300];

slice3Dmatrix(avgim(:, :, :, 1), plot_pars)

plot_pars.range = [0 150];
slice3Dmatrix(avgim(:, :, :, 2), plot_pars)

% 3.2) selected planes
planes2plot = 130:140;
plot_pars.lag = 0.5;
plot_pars.lag = 0.01;
plot_pars.sizY = [size(avgim, 1), size(avgim, 2), length(planes2plot)];
%plot_pars.range = [0 130];
plot_pars.range = [0 150];

slice3Dmatrix(avgim(:, :, planes2plot, 1), plot_pars)

%plot_pars.range = [0 350];
plot_pars.range = [0 500];
slice3Dmatrix(avgim(:, :, planes2plot, 2), plot_pars)

%% 4) replace planes (using info from 2) and 3))

%plane2replace = planes2plot([]);
plane2replace = find(isinf(M(:, 1)))';
channel2use = 1;
avgim(:, :, :, channel2use) = ...
    framegapfill(plane2replace, avgim(:, :, :, channel2use));

%plane2replace = planes2plot([]);
plane2replace = find(isinf(M(:, 2)))';
channel2use = 2;
avgim(:, :, :, channel2use) = ...
    framegapfill(plane2replace, avgim(:, :, :, channel2use));
