function batch_plot_compare_fictrac_results(...
    ifiles_1, ifiles_2, fignames, iparams)
% batch_plot_compare_fictrac_results: compare basic stats
%   from fictrac tracking between selected files
%
% Usage:
%   batch_plot_compare_fictrac_results(...
%       ifiles_1, ifiles_2, fignames, iparams)
%
% Args:
%   ifiles_1: name pattern of files to use
%   ifiles_2: name pattern of files to use
%   fignames: figure names
%   iparams: parameters to update
%       (fmetsuffix: suffix of metadata file)
%           (default, '_metadata.mat')
%       (range: range for plotting)
%           (default, [0 1])
%       (refcha: reference channel)
%           (default, 1)
%       (dir_depth: depth of directory search)
%           (default, 0)
%       (oDir: target directory to save figures)
%           (default, [pwd, filesep, 'fictrac'])
%       (ball_radious: depth of directory search)
%           (default, 0)
%       (sig: std of gaussian kernel to smooth motor variable)
%           (default, 3)
%       (siz: size of kernel  to smooth motor variable)
%           (default, 10)
%       (hbins_speed: depth of directory search)
%           (default, -20:1:60)
%       (framerate: frame rate)
%           (default, 100)
% Notes

% default params
metpars = [];
metpars.fmetsuffix = '_metadata.mat';
metpars.fi2reject = 'Zstack';
metpars.range = [0 1];
metpars.refcha = 1;
metpars.dir_depth = 0;
metpars.oDir = [pwd, filesep, 'fictrac_compare'];
metpars.ball_radious = 4.5;
metpars.sig = 3;
metpars.siz = 10;
metpars.hbins_speed = -20:1:60;
metpars.framerate = 100;
metpars.conditiostr = {'condition-1', 'condition-2'};

if ~exist(metpars.oDir, 'dir')
   mkdir(metpars.oDir)
end

% update variables
if ~exist('iparams', 'var'); iparams = []; end
metpars = loparam_updater(metpars, iparams);

if length(ifiles_1) ~= length(ifiles_2)
   fprintf('Input files differ in number\n') 
   return
end

fprintf(['Generating plots for # ', ...
    num2str(length(ifiles_1)), ' pairs\n'])

% setup figure (overall results)
figH = figure('Position', genfigpos(1, 'center', [400 400]));
axH(1) = subplot(1, 1, 1);

colorvect = colorGradient([1 1 1], ...
    [0 0 0], numel(ifiles_1) + 3);
colorvect = colorvect(4:end, :);

for i = 1:numel(ifiles_1)
    
    fprintf('\n')
    
    out = run_per_pair(ifiles_1{i}, ifiles_2{i}, ...
        fignames{i}, metpars);
    
    % plot mean speed across files
    lineH(i) = plot([1:(numel(out{1}) + numel(out{2}))], ...
        [out{1} out{2}], '.-', 'Color', colorvect(i, :), ...
        'Parent', axH(1));
    hold(axH(1), 'on')
    
end

% edit figure
axH(1).XTick = 1:numel(metpars.conditiostr);
axH(1).XTickLabels = metpars.conditiostr;
axH(1).XTickLabelRotation = 45;
axH(1).XLabel.String = 'Manipulation';
axH(1).YLabel.String = 'Average speed (mm/s)';
axH(1).XLim = [0.5 (numel(metpars.conditiostr) + 0.5)];
legend(axH(1), lineH, strrep(ifiles_1, '_', '-'), ...
    'Location', 'northeast')

savefig_int(figH, metpars.oDir, ...
    'summary_pair_locomotion_stats', ...
    [1 0 0 0 0 0 0 0 1 0 1])

close(figH)

fprintf('... Done\n')

end

function out = run_per_pair(ifile_1, ifile_2, ...
    figurename, metpars)
% run_per_pair: function that runs each pair
%
% Usage:
%   out = run_per_pair(ifile_1, ifile_2, ...
%       figurename, metpars)
%
% Args:
%   ifile_1: file name
%   ifile_2: file name
%   metpars: input parameters

% get all possible files
if metpars.dir_depth == 0
    f2run = rdir(['.', filesep, '*', ...
        metpars.fmetsuffix]);
elseif metpars.dir_depth == 1
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', metpars.fmetsuffix]);
else
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', filesep, '*', ...
        metpars.fmetsuffix]);
end

f2run = str2rm(metpars.fi2reject, f2run);

% define ifile_1 to use
filename_1 = str2match(ifile_1, f2run);
filename_1 = {filename_1.name}';
[filename_1, iDir_1] = split_path(filename_1);
filename_1 = strrep(filename_1, metpars.fmetsuffix, '');

fprintf(['Using # ', ...
    num2str(length(filename_1)), ' as group one\n'])

% define ifile_1 to use
filename_2 = str2match(ifile_2, f2run);
filename_2 = {filename_2.name}';
[filename_2, iDir_2] = split_path(filename_2);
filename_2 = strrep(filename_2, metpars.fmetsuffix, '');

fprintf(['Using # ', ...
    num2str(length(filename_2)), ' as group two\n'])

% setup figure-1
figH = figure('Position', genfigpos(1, 'center', [500 900]));
axH(1) = subplot(4, 2, 1);
axH(2) = subplot(4, 2, 2);
axH(3) = subplot(4, 2, [3 4]);
axH(4) = subplot(4, 2, [5 6]);
axH(5) = subplot(4, 2, [7 8]);

% setup color
colorvect_1 = colorGradient([1 1 1], ...
    [0 0 1], numel(filename_1) + 3);
colorvect_1 = colorvect_1(4:end, :);
colorvect_2 = colorGradient([1 1 1], ...
    [1 0 0], numel(filename_2) + 3);
colorvect_2 = colorvect_2(4:end, :);

k = 1;
for i = 1:length(filename_1)
   
    [xy_range_1(i, 1), lineH(k), out{1}(i)] = ...
        plot_pair_pair(filename_1{i}, iDir_1{i}, ...
        metpars, colorvect_1(i, :), axH);
    legend_str{k, 1} = 'group-1';
    k = k + 1;
    
end

for i = 1:length(filename_2)
    
    [xy_range_2(i, 1), lineH(k), out{2}(i)] = ...
        plot_pair_pair(filename_2{i}, iDir_2{i}, ...
        metpars, colorvect_2(i, :), axH);
    legend_str{k, 1} = 'group-2';
    k = k + 1;
    
end

xy_range = [xy_range_1; xy_range_2];

axH(1).XLabel.String = 'X distance (mm)';
axH(1).YLabel.String = 'Y distance (mm)';
axH(1).XLim = [-max(xy_range) max(xy_range)];
axH(1).YLim = [-max(xy_range) max(xy_range)];

axH(2).XLabel.String = 'X distance (mm)';
axH(2).YLabel.String = 'Y distance (mm)';
axH(2).XLim = [-max(xy_range_2) max(xy_range_2)];
axH(2).YLim = [-max(xy_range_2) max(xy_range_2)];

axH(3).XLabel.String = 'speed (mm/s)';
axH(3).YLabel.String = 'Probability';

axH(4).XLabel.String = 'f-speed (mm/s)';
axH(4).YLabel.String = 'Probability';

axH(5).XLabel.String = 'l-speed (mm/s)';
axH(5).YLabel.String = 'Probability';

% save figure
legend(axH(5), lineH, legend_str, 'Location', 'northeast')

savefig_int(figH, metpars.oDir, ...
    [figurename, '_pair_locomotion_stats'], ...
    [1 0 0 0 0 0 0 0 1 0 1])

close(figH)

end

function [xy_range, lineH, mean_speed] = ...
    plot_pair_pair(filename, iDir, ...
    metpars, colorvect, axH)
% run_per_pair: function that runs each pair
%
% Usage:
%   [xy_range, lineH, mean_speed] = ...
%       run_per_pair(filename, iDir, ...
%       metpars, colorvect, axH)
%
% Args:
%   filename: file name
%   iDir: directory name
%   metpars: input parameters
%   colorvect: color
%   axH: axis handle

load([iDir, filesep, filename, ...
        metpars.fmetsuffix], 'wDat');

if exist('wDat', 'var') && ...
        isfield(wDat, 'vid') && ~isempty(wDat.vid.var)

    [xy_range, lineH, mean_speed] = ...
        plot_vid_results(wDat, metpars, ...
        colorvect, axH);

end

clear wDat

end

function [xy_range, lineH, mean_speed] = ...
    plot_vid_results(wDat, metpars, ...
    colorvect, axH)
% run_per_pair: function that runs each pair
%
% Usage:
%   [xy_range, lineH, mean_speed] = ...
%       run_per_pair(wDat, metpars, ...
%       colorvect, axH)
%
% Args:
%   filename: file name
%   iDir: directory name
%   metpars: input parameters
%   colorvect: color
%   axH: axis handle

[Y, ~, ~] = load_edit_fictrac_motor_var(wDat, metpars.ball_radious, ...
    {'speed', 'fV', 'lV'}, metpars.sig, metpars.siz, 0);
x_mm = wDat.vid.var{1}(:, 15)*metpars.ball_radious;
y_mm = wDat.vid.var{1}(:, 16)*metpars.ball_radious;
mean_speed = mean(Y(1, :));

% plot XY
plot(x_mm, y_mm, ...
    'Color', colorvect, ...
    'Parent', axH(1))

hold(axH(1), 'on')

% plot XY
plot(x_mm, y_mm, ...
    'Color', colorvect, ...
    'Parent', axH(2))

hold(axH(2), 'on')

xy_range(1, 1) = max(abs([x_mm; y_mm]));

% plot instantaneous speed
[y_, ~] = hist(Y(1, :), metpars.hbins_speed);

plot(metpars.hbins_speed, y_/sum(y_), ...
    'Color', colorvect, ...
    'Parent', axH(3))

hold(axH(3), 'on')

% plot instantaneous forward speed
[y_, ~] = hist(Y(2, :), metpars.hbins_speed);

lineH = plot(metpars.hbins_speed, y_/sum(y_), ...
    'Color', colorvect, ...
    'Parent', axH(4));

hold(axH(4), 'on')

% plot instantaneous forward speed
[y_, ~] = hist(Y(3, :), metpars.hbins_speed);

lineH = plot(metpars.hbins_speed, y_/sum(y_), ...
    'Color', colorvect, ...
    'Parent', axH(5));

hold(axH(5), 'on')

end
