function batch_plot_mot_var_stats(FolderName, FileName, iparams)
% batch_plot_mot_var_stats: plot motor variable stats
%  (variance explained by mean) and autocorrelation.
%
% Usage:
%   batch_plot_mot_var_stats(FolderName, FileName, iparams)
%
% Args:
%   FolderName: folders to use
%   FileName: files to use
%   iparams: parameters
%       (cDir: current directory)
%       (fo2reject: folders to reject)
%       (fi2reject: files to reject)
%       (fsuffix: suffix of files with roi varaiable)
%           (default, '_prosroi')
%       (metsuffix: suffix of files with wDat variable)
%           (default, '_prosmetadata')
%       %%%%%%%%%%%% filter related %%%%%%%%%%%%
%       (mot2use: motor variables to use)
%           (default, {'speed', 'fV', 'lV', 'yaw', 'pitch', 'roll'})
%       %%%%%%%%%%%% parpool & server related %%%%%%%%%%%%
%       (serId: server id)
%           (default, 'int')
%       (corenum: number of cores)
%           (default, 4)
%       (chunksiz: number of chunks for parpool)
%           (default, 80)
%       %%%%%%%%%%%% motor variable related %%%%%%%%%%%%
%       (ball_radious: radious of ball (mm))
%           (default, 9)
%       (sig: std of gaussian kernel to smooth motor variable)
%       (siz: size of kernel  to smooth motor variable)
%       %%%%%%%%%%%% save summary plots %%%%%%%%%%%%
%       (oDir: output directory to save summary results)
%           (default, [pwd, filesep, 'motvar_stats'])
%
% Notes:

% Default params
pSLM = []; 
pSLM.cDir = pwd;
pSLM.fo2reject = {'.', '..', 'preprocessed', 'BData'}; 
pSLM.fi2reject = {'Zstack'};
pSLM.metsuffix = '_metadata';
pSLM.mot2use = {'speed', 'fV', 'lV', 'yaw', 'pitch', 'roll'};
pSLM.serverid = 'int';
pSLM.corenum = 4;
pSLM.chunksiz = 80;
pSLM.ball_radious = 9;
pSLM.sig = 3;
pSLM.siz = 10;

% update variables
if ~exist('FileName', 'var'); FileName = []; end
if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('iparams', 'var'); iparams = []; end
pSLM = loparam_updater(pSLM, iparams);

% start pararell pool if not ready yet
ppobj = setup_parpool(pSLM.serverid, pSLM.corenum);

% Selecting folders
f2run = dir;
f2run = str2match(FolderName, f2run);
f2run = str2rm(pSLM.fo2reject, f2run);
f2run = {f2run.name};
fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']); 
    
    cd(f2run{i});
    
    pSLM.oDir = [pwd, filesep, 'motvar_stats'];
    if ~exist(pSLM.oDir, 'dir')
        mkdir(pSLM.oDir)
    end
    
    runperfolder(FileName, pSLM);
    cd(pSLM.cDir);
        
    fprintf('\n')
    
end

delete_parpool(ppobj);

fprintf('... Done\n')

end

function runperfolder(filename, pSLM)
% runperfolder: run all files per folder
%
% Usage:
%   runperfolder(filename)
%
% Args:
%   filename: file name pattern
%   pSLM: internal plotting variable

% Files to load
f2run = rdir(['.', filesep, '*', pSLM.metsuffix, '*.mat']);
f2run = str2match(filename, f2run); 
f2run = {f2run.name};

fprintf('Plotting motor variables\n')
fprintf(['Running n-files : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running file : ', strrep(f2run{i}, filesep, ' '), '\n']);
    try
        plotperfile(f2run{i}, pSLM);
    catch
        fprintf(['Failed : ', strrep(f2run{i}, filesep, ' '), '\n'])
    end
    
end

fprintf('****** Done ******\n')

end

function plotperfile(filename, pSLM)
% plotperfile: run stimuli modulation for each file
%
% Usage:
%   plotperfile(filename, pSLM)
%
% Args:
%   filename: file name
%   pSLM: internal plotting variable

t0 = stic;

% 1) load required variables
%   compatible with old format
load(filename, 'wDat')

% 2) load motor variables directly from fictrac (fV, lV, V, yaw, pitch, roll)
motor_fictrac = load_edit_fictrac_motor_var(wDat, pSLM.ball_radious, ...
    pSLM.mot2use, pSLM.sig, pSLM.siz, 0);

% 3) load additional motor variables SVD of cropped video
[motor_SVD, imtime, vidtime, SVD_Dat] = ...
    load_edit_video_SVD(wDat, strrep(filename, '_metadata.mat', '_proc.mat'), ...
    pSLM.sig, pSLM.siz, 0);

% figure settings
color_1 = {'b', 'c', rgb('RoyalBlue'), 'r', rgb('Crimson'), rgb('DarkRed')};
color_2 = [rgb('OrangeRed'); rgb('Orange'); rgb('Green'); rgb('YellowGreen'); ...
    colorGradient([0 0 0], [0.8 0.8 0.8], size(motor_SVD, 1))];

% 4) histograms
vel_bins = -200:0.1:200;
deg_bins = -300:0.1:300;

for i = 1:3
    [y_vel(i, :), ~] = hist(motor_fictrac(i, :), vel_bins);
    [y_deg(i, :), ~] = hist(motor_fictrac(i+3, :), deg_bins);
end

sd_bins = -8:0.1:8;
for i = 1:size(motor_SVD, 1)
    [y_sd(i, :), ~] = hist(motor_SVD(i, :)/std(motor_SVD(i, :)), sd_bins);
end

% plot histograms
plot_hist(vel_bins, deg_bins, sd_bins, ...
    y_vel, y_deg, y_sd, color_1, color_2, filename, pSLM.oDir)
clear vel_bins deg_bins sd_bins ...
    y_vel y_deg y_sd

% 5) calculate crosscorrelations
motor_all = [motor_fictrac; motor_SVD];
lag_t = 300;
auto_corr = mat2cell(motor_all, ones(1, size(motor_all, 1)), size(motor_all, 2));
auto_corr = cellfun(@(x) xcorr(x, lag_t, 'normalized'), auto_corr, 'uniformoutput', false);
auto_corr = cell2mat(auto_corr);

% plot crosscorrelations
plot_croscorr(auto_corr, lag_t, color_1, color_2, ...
    filename, pSLM.oDir)
clear auto_corr

% 6) resample to imaging resolution
motor_all = interp1(vidtime, motor_all', imtime, 'linear');
motor_all = framegapfill(find(isnan(motor_all(:, 1))), motor_all');
motor_all(1, :) = sqrt(motor_all(2, :).^2 + motor_all(3, :).^2);

% 7) get trial structure and calculate trial average
% Get number of stimuli presented
stim_ntype = 1;

if ~contains(wDat.datatype, 'opto') || contains(wDat.datatype, 'prv')
    
    % get stimuli name
    stim_name = cellfun(@(x, y) [x '_AMP_' num2str(y)], ...
        wDat.sPars.name(1, :), chunk2cell(wDat.sPars.int(1, :), 1), ...
        'UniformOutput', false);
    [stim_name, ~, stim_u_idx] = unique(stim_name, 'stable');
     
    % get number of reps per stim
    stim_all_idx = stim_u_idx(wDat.sPars.order(1, :));
    stim_all_idx = stim_all_idx(1:size(wDat.sTime, 1), 1);
    stim_count = hist(stim_all_idx, 1:1:max(stim_all_idx(:)));
    stim_ntype = numel(unique(stim_name));
    
    % get stimulus post baseline time to define filter length
    for i = 1:max(stim_u_idx)
        stim_post_lag(i, 1) = min(wDat.sPars.basPost(stim_u_idx == i));
    end
    clear stim_u_idx i
       
end

[motor_per_trial, time_rel, ~, ~, ~, ~, ~, ~, stim_idx, ~] = trace2trials(wDat, ...
    motor_all, [-10 10], stim_all_idx, 1, []);

% plot mean per trial
stim2plot = unique(stim_all_idx)';

plot_trial_med(motor_per_trial, time_rel, stim_name, ...
    stim2plot, stim_idx, color_1, color_2, filename, pSLM.oDir)

% plot pixel components
% plot mean image, mean energy, first 4 components and traces

plot_SVD_spatial_temporal(SVD_Dat, imtime, ...
    motor_all(7:10, :), filename, pSLM.oDir, color_2)

plot_fictrac_temporal(motor_all(1:6, :), ...
    imtime, filename, pSLM.oDir, color_1)

end

function plot_fictrac_temporal(fictrac_temp_res, ...
    fTime, filename, oDir, color_1)
% plot_fictrac_temporal: plot fictract variables
%
% Usage:
%   plot_fictrac_temporal(fictrac_temp_res, ...
%       fTime, filename, oDir, color_1)
%
% Args:
%   fictrac_temp_res: fictrac variables resampled to fTime
%   fTime: frame time (imaging resolution)
%   filename: filename
%   oDir: outpud directory
%   color_2: colormap

[figH, axH] = makefigs(2, 3, [1000 450], 'center');
axH(1) = subplot(2, 3, [1 2]);
axH(2) = subplot(2, 3, 3);
axH(3) = subplot(2, 3, [4 5]);
axH(4) = subplot(2, 3, 6);

for i = 1:3
    plot(fTime, fictrac_temp_res(i, :), ...
        'Color', color_1{i}, 'Parent', axH(1));
    hold(axH(1), 'on')
end
axH(1).XLim = fTime([1 end]);

for i = 1:3
    plot(fTime, fictrac_temp_res(i, :), ...
        'Color', color_1{i}, 'Parent', axH(2));
    hold(axH(2), 'on')
end
axH(2).XLim = [-10 50];

for i = 4:6
    plot(fTime, fictrac_temp_res(i, :), ...
        'Color', color_1{i}, 'Parent', axH(3));
    hold(axH(3), 'on')
end
axH(3).XLim = fTime([1 end]);

for i = 4:6
    plot(fTime, fictrac_temp_res(i, :), ...
        'Color', color_1{i}, 'Parent', axH(4));
    hold(axH(4), 'on')
end
axH(4).XLim = [-10 50];

axH(1).YLabel.String = 'ground velocities (mm/s)';
axH(3).YLabel.String = 'angular velocities (rad/s)';
axH(3).XLabel.String = 'Time (s)';
axH(4).XLabel.String = 'Time (s)';

figformat = [1 0 0 0 0 0 0 0 1 0 1];
fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;

save_edit_fig_int(axH, figH, oDir, ...
    [strrep(filename, '_metadata.mat', ''), ...
    '_motvar_fictrac_temp'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH)

end

function plot_SVD_spatial_temporal(SVD_Dat, ...
    fTime, SVD_temp_res, filename, oDir, color_2)
% plot_SVD_spatial_temporal: plot SVD components
%
% Usage:
%   plot_SVD_spatial_temporal(SVD_Dat, ...
%       fTime, SVD_temp_res, filename, oDir, color_2)
%
% Args:
%   SVD_Dat: parameter variable from SVD analysis
%   fTime: frame time (imaging resolution)
%   SVD_temp_res: SVD variables resampled to fTime
%   filename: filename
%   oDir: outpud directory
%   color_2: colormap

np = [0 SVD_Dat.npix];
np = cumsum(np);
spa_c2plot = 4;
k = 1;

[figH, axH] = makefigs(4, 3, [1000 900], 'center');
axH(7) = subplot(4, 3, [7 8]);
axH(8) = subplot(4, 3, 9);
axH(9) = subplot(4, 3, [10 11]);
axH(10) = subplot(4, 3, 12);

% calculate bounding box
[x_, y_] = ind2sub(size(SVD_Dat.wpix{k}), find(SVD_Dat.wpix{k} ~= 0));
box = [min(x_) min(y_); max(x_) max(y_)];

% plot average image
im = reshape(SVD_Dat.avgframe, size(SVD_Dat.wpix{1}));
im = im(box(1):box(2), box(3):box(4));

imagesc(im, 'alphadata', ~isnan(im), 'Parent', axH(1));
axis(axH(1), 'image');
colormap(axH(1), 'gray'); 
colorbar(axH(1))
axH(1).XTick = [];
axH(1).YTick = [];

% plot energy image
im = reshape(SVD_Dat.avgmotion, size(SVD_Dat.wpix{1}));
im = im(box(1):box(2), box(3):box(4));

imagesc(im, 'alphadata', ~isnan(im), 'Parent', axH(4));
hold(axH(4), 'on')
contour(SVD_Dat.wpix{1}(box(1):box(2), box(3):box(4)), 'k', 'Parent', axH(4))
axis(axH(4), 'image');
colormap(axH(4), colorGradient([1 0 0], [1 1 1])); 
colorbar(axH(4))
caxis(axH(4), [0 5])
axH(4).XTick = [];
axH(4).YTick = [];

ax2use = [2 3 5 6];
for i = 1:spa_c2plot
    i1 = SVD_Dat.uMotMask{1}(np(k)+[1:SVD_Dat.npix(k)], i);
    ib = NaN*zeros(floor(SVD_Dat.nY{k}/SVD_Dat.sc), floor(SVD_Dat.nX{k}/SVD_Dat.sc));
    ib(SVD_Dat.wpix{k}) = i1;
    ib = ib(box(1):box(2),box(3):box(4));
    imagesc(ib, 'alphadata', ~isnan(ib), 'Parent', axH(ax2use(i)));
    axis(axH(ax2use(i)), 'image');
    caxis(axH(ax2use(i)), [-0.1 0.1])
    colormap(axH(ax2use(i)), 'redblue')
    colorbar(axH(ax2use(i)))
    axH(ax2use(i)).XTick = [];
    axH(ax2use(i)).YTick = [];
    axH(ax2use(i)).XColor = color_2(i, :);
    axH(ax2use(i)).YColor = color_2(i, :);
    axH(ax2use(i)).LineWidth = 2;
end

for i = 1:size(SVD_temp_res, 1)
    if i ~= 1
        plot(fTime, SVD_temp_res(i, :), ...
            'Color', color_2(i, :), 'Parent', axH(7));
        hold(axH(7), 'on')
    end
end
plot(fTime, SVD_temp_res(1, :), ...
    'Color', color_2(1, :), 'Parent', axH(7));
axH(7).XLim = fTime([1 end]);

for i = 1:size(SVD_temp_res, 1)
    if i ~= 1
        plot(fTime, SVD_temp_res(i, :), ...
            'Color', color_2(i, :), 'Parent', axH(8));
        hold(axH(8), 'on')
    end
end
plot(fTime, SVD_temp_res(1, :), ...
    'Color', color_2(1, :), 'Parent', axH(8));
axH(8).XLim = [-10 50];

% plot szcore
SVD_temp_res = zscorebigmem(SVD_temp_res);
for i = 1:size(SVD_temp_res, 1)
    if i ~= 1
        plot(fTime, SVD_temp_res(i, :), ...
            'Color', color_2(i, :), 'Parent', axH(9));
        hold(axH(9), 'on')
    end
end
plot(fTime, SVD_temp_res(1, :), ...
    'Color', color_2(1, :), 'Parent', axH(9));
axH(9).XLim = fTime([1 end]);

for i = 1:size(SVD_temp_res, 1)
    if i ~= 1
        plot(fTime, SVD_temp_res(i, :), ...
            'Color', color_2(i, :), 'Parent', axH(10));
        hold(axH(10), 'on')
    end
end
plot(fTime, SVD_temp_res(1, :), ...
    'Color', color_2(1, :), 'Parent', axH(10));
axH(10).XLim = [-10 50];

axH(7).YLabel.String = 'temporal SVD (a.u)';
axH(9).YLabel.String = 'temporal SVD (SD)';
axH(9).XLabel.String = 'Time (s)';
axH(10).XLabel.String = 'Time (s)';

axH(1).Title.String = 'average image';
axH(2).Title.String = 'SVD-1';
axH(3).Title.String = 'SVD-2';
axH(4).Title.String = 'average motion energy';
axH(5).Title.String = 'SVD-3';
axH(6).Title.String = 'SVD-4';

figformat = [1 0 0 0 0 0 0 0 1 0 1];
fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;

save_edit_fig_int(axH([1 4 7:10]), figH, oDir, ...
    [strrep(filename, '_metadata.mat', ''), ...
    '_motvar_SVD_temp_spa'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH)

end

function plot_trial_med(motor_per_trial, time_rel, stim_name, ...
    stim2plot, stim_idx, color_1, color_2, filename, oDir)
% plot_trial_med: plot median of each motor variable across trials
%
% Usage:
%   plot_trial_med(motor_per_trial, time_rel, ...
%   	stim2plot, stim_idx, color_1, color_2, filename, oDir)
%
% Args:
%   motor_per_trial: motor variable per trial
%   time_rel: frame time (imaging resolution)
%   stim_name: stimuli name
%   stim2plot: stimuli to plot
%   stim_idx: indeces of each stimuli
%   color_1: colormap
%   color_2: colormap
%   filename: filename
%   oDir: outpud directory

figH = figure('Position', genfigpos(1, 'center', [400 900]));
for i = 1:numel(stim2plot)
    axH(i) = subplot(numel(stim2plot), 1, i);
end

for i = 1:numel(stim2plot)
    
    % normalize
    motor_per_trial_mean = cf(@(x) x(:, stim_idx == stim2plot(i)), motor_per_trial);
    motor_per_trial_mean = cf(@(x, y) reshape(zscorebigmem(x(:)'), y), motor_per_trial_mean, ...
        cf(@(x) size(x), motor_per_trial_mean));
    motor_per_trial_mean = cell2mat(cf(@(x) median(x, 1), motor_per_trial_mean));
    
    for j = 8:size(motor_per_trial_mean, 1)
        plot(time_rel(stim_idx == 1), motor_per_trial_mean(j, :), ...
            'Color', color_2(j - 6, :), 'Parent', axH(i))
        hold(axH(i), 'on')
    end
    
    plot(time_rel(stim_idx == 1), motor_per_trial_mean(7, :), ...
        'Color', color_2(7 - 6, :), 'Parent', axH(i))
    
    for j = 1:6
        plot(time_rel(stim_idx == 1), motor_per_trial_mean(j, :), ...
            'Color', color_1{j}, 'Parent', axH(i))
    end
   
    axH(i).Title.String = stim_name{stim2plot(i)};
    
end

axH(round(numel(stim2plot)/2)).YLabel.String = 'motor variable change (SD)';
axH(end).XLabel.String = 'Time (s)';

figformat = [1 0 0 0 0 0 0 0 1 0 1];
fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;

save_edit_fig_int(axH, figH, oDir, ...
    [strrep(filename, '_metadata.mat', ''), '_motvar_trial_med'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH)

end

function plot_croscorr(temp_, lag_t, ...
    color_1, color_2, filename, oDir)
% plot_croscorr: plot autocorrelation across variables
%
% Usage:
%   plot_croscorr(temp_, lag_t, ...
%      color_1, color_2, filename, oDir)
%
% Args:
%   temp_: motor variables
%   lag_t: time lags to use for autocorrelation
%   color_1: colormap
%   color_2: colormap
%   filename: filename
%   oDir: outpud directory

figH = figure('Position', genfigpos(1, 'center', [300 300]));
axH(1) = subplot(1, 1, 1);

for i = 5:(size(temp_, 1) - 6)
    plot((-lag_t:lag_t)*0.01, temp_(i + 6, :), ...
        'Color', color_2(i, :), 'Parent', axH(1))
    hold(axH(1), 'on')
end

for i = 1:4
    plot((-lag_t:lag_t)*0.01, temp_(i + 6, :), ...
        'Color', color_2(i, :), 'Parent', axH(1))
end

for i = 1:6
    plot((-lag_t:lag_t)*0.01, temp_(i, :), ...
        'Color', color_1{i}, 'Parent', axH(1))
    hold(axH(1), 'on')
end

axH(1).XLim = [-3 3];
axH(1).XLabel.String = 'lag (s)';
axH(1).YLabel.String = 'autocorrelation';

figformat = [1 0 0 0 0 0 0 0 1 0 1];
fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;

save_edit_fig_int(axH, figH, oDir, ...
    [strrep(filename, '_metadata.mat', ''), '_motvar_autocorr'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH)

end

function plot_hist(vel_bins, deg_bins, sd_bins, ...
    y_vel, y_deg, y_sd, color_1, color_2, filename, oDir)
% plot_hist: plot histograms across variables
%
% Usage:
%   plot_hist(vel_bins, deg_bins, sd_bins, ...
%       y_vel, y_deg, y_sd, color_1, color_2, filename, oDir)
%
% Args:
%   vel_bins: velocity bins (mm/s)
%   deg_bins: angular bins (rads/s)
%   sd_bins: SD bins
%   y_vel: velocity histogram
%   y_deg: angular histogram
%   y_sd: SD histogram
%   color_1: colormap
%   color_2: colormap
%   filename: filename
%   oDir: outpud directory

figH = figure('Position', genfigpos(1, 'center', [1200 300]));

axH(1) = subplot(1, 3, 1);
axH(2) = subplot(1, 3, 2);
axH(3) = subplot(1, 3, 3);

for i = 1:size(y_sd, 1)
    
    if i > 4
        plot(sd_bins, y_sd(i, :)/sum(y_sd(i, :)), ...
            'Color', color_2(i, :), 'Parent', axH(3))
    end
    hold(axH(3), 'on')
    
end

for i = 1:4
    
    linH_3(i) = plot(sd_bins, y_sd(i, :)/sum(y_sd(i, :)), ...
        'Color', color_2(i, :), 'Parent', axH(3));        
    
end

for i = 1:3
    
    linH_1(i) = plot(vel_bins, y_vel(i, :)/sum(y_vel(i, :)), ...
        'Color', color_1{i}, 'Parent', axH(1));
    hold(axH(1), 'on')
    
    linH_2(i) = plot(deg_bins, y_deg(i, :)/sum(y_deg(i, :)), ...
        'Color', color_1{i + 3}, 'Parent', axH(2));
    hold(axH(2), 'on')
    
end

legend(axH(1), linH_1, {'speed', 'fV', 'lV'})
legend(axH(2), linH_2, {'yaw', 'pitch', 'roll'})
legend(axH(3), linH_3, {'SVD-1', 'SVD-2', 'SVD-3', 'SVD-4'})

axH(1).XLim = [-15 15];
axH(2).XLim = [-4 4];
axH(3).XLim = [-5 5];
axH(1).XLabel.String = 'velocity mm/s';
axH(2).XLabel.String = 'velocity rads/s';
axH(3).XLabel.String = 'SD';
axH(1).YLabel.String = 'Probability';
axH(2).YLabel.String = 'Probability';
axH(3).YLabel.String = 'Probability';
axH(1).Title.String = 'ground velocity';
axH(2).Title.String = 'angular velocity';
axH(3).Title.String = 'pixel SVDs';

figformat = [1 0 0 0 0 0 0 0 1 0 1];
fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;

save_edit_fig_int(axH, figH, oDir, ...
    [strrep(filename, '_metadata.mat', ''), '_motvar_hist'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH)

end
