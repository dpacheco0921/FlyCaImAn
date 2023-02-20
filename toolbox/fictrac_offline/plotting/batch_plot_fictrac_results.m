function batch_plot_fictrac_results(FolderName, Filename, iparams)
% batch_plot_fictrac_results: plot basic stats from raw fictrac tracking, and
%   and from locomotion parsed according to stimuli presentation.
%
% Usage:
%   batch_plot_fictrac_results(FolderName, Filename, iparams)
%
% Args:
%   Foldername: name pattern of directories to use
%       (default, [])
%   Filename: name pattern of files to use
%       (default, [])
%   iparams: parameters to update
%       (fmetsuffix: suffix of metadata file)
%           (default, '_metadata.mat')
%       (range: range for plotting)
%           (default, [0 1])
%       (dir_depth: depth of directory search)
%           (default, 0)
%       (oDir: target directory to save figures)
%           (default, [pwd, filesep, 'fictrac'])
%       (ball_radious: depth of directory search)
%           (default, 4.5)
%       (hbins_speed: bins to use for histogram)
%           (default, [-20:0.1:60])
%       (framerate: frame rate)
%           (default, 100)
%       (time_buffer: time before and after stimulus)
%           (default, [-5 +5])
%       (timeres: temporal resolution in sec)
%           (default, 0.1)
%
% Notes
%   this setup x = yaw, y = pitch, z = roll (this needs to be tested for every setup)
%   this setup y = forward and x = lateral (this needs to be tested for every setup)

% default params
metpars = [];
metpars.fmetsuffix = '_metadata.mat';
metpars.fi2reject = 'Zstack';
metpars.range = [0 1];
metpars.dir_depth = 0;
metpars.oDir = [pwd, filesep, 'fictrac'];
metpars.ball_radious = 4.5;
metpars.hbins_speed = -20:0.1:60;
metpars.framerate = 100;
metpars.time_buffer = [-5 +5];
metpars.timeres = 0.1;

if ~exist('Filename', 'var')
    Filename = [];
end

if ~exist('FolderName', 'var')
    FolderName = [];
end

% update variables
if ~exist('iparams', 'var'); iparams = []; end
metpars = loparam_updater(metpars, iparams);

if ~exist(metpars.oDir, 'dir')
   mkdir(metpars.oDir)
end

% define files to use
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
f2run = {f2run.name}';
[filename, iDir] = split_path(f2run);
filename = strrep(filename, metpars.fmetsuffix, '');

if ~isempty(Filename)
    f2run = find(contains(filename, Filename));
    filename = filename(f2run);
    iDir = iDir(f2run);
end

if ~isempty(FolderName)
    f2run = find(contains(iDir, FolderName));
    filename = filename(f2run);
    iDir = iDir(f2run);
end

fprintf(['Generating plots for ', ...
    num2str(numel(filename)), ' files\n'])

% plot fictrac + stimuli results
for i = 1:numel(filename)

    load([iDir{i}, filesep, filename{i}, ...
        metpars.fmetsuffix], 'wDat');

    if exist('wDat', 'var') && ...
            isfield(wDat, 'vid') && ~isempty(wDat.vid.var)
        
        plot_vid_results(wDat, filename{i}, metpars)
        
        if isfield(wDat, 'sPars') && ~isempty(wDat.sPars.name)
            plot_vid_stim_results(wDat, filename{i}, metpars)
        end
        
    end
    
    clear wDat
    
end

fprintf('... Done\n')

end

function plot_vid_results(wDat, fname, metpars)
% plot_vid_results: plot basic locomotion stats
%
% Usage:
%   plot_vid_results(wDat, fname, metpars)
%
% Args:
%   wDat: name pattern of directories to use
%       (default, [])
%   fname: figure name
%       (default, [])
%   metpars: parameters to update
%       (oDir: target directory to save figures)
%           (default, [pwd, filesep, 'fictrac'])
%       (ball_radious: depth of directory search)
%           (default, 4.5)
%       (hbins_speed: bins to use for histogram)
%           (default, [-20:0.1:60])
%       (framerate: frame rate)
%           (default, 100)
%
% Notes
%   this setup x = yaw, y = pitch, z = roll (this needs to be tested for every setup)
%   this setup y = forward and x = lateral (this needs to be tested for every setup)

% plot basic motion stats
figH = figure('Position', genfigpos(1, 'center', [900 600]));
axH(1) = subplot(2, 2, 1);
axH(2) = subplot(2, 2, 3:4);
axH(3) = subplot(2, 2, 2);

for i = 1:numel(wDat.vid.var)
    
    x_mm = wDat.vid.var{i}(:, 15)*metpars.ball_radious;
    y_mm = wDat.vid.var{i}(:, 16)*metpars.ball_radious;
    speed_mm_s = wDat.vid.var{i}(:, 19)*metpars.ball_radious*metpars.framerate;
    fv_mm_s = diff(y_mm)*metpars.framerate;
    lv_mm_s = diff(x_mm)*metpars.framerate;
    time_ = wDat.vid.fstEn{i}(:, 1);
    
    % plot XY
    plot(x_mm, y_mm, 'Color', 'k', ...
        'Parent', axH(1))
    
    xy_range(i, 1) = max(abs([x_mm; y_mm]));
     
    % plot instantaneous speed
    plot(time_, speed_mm_s, 'Color', 'b', ...
        'Parent', axH(2))
    hold(axH(2), 'on')
    
    % plot instantaneous forward speed
    plot(time_(2:end), fv_mm_s, 'Color', 'r', ...
        'Parent', axH(2))
    
    % plot instantaneous forward speed
    plot(time_(2:end), lv_mm_s, 'Color', 'm', ...
        'Parent', axH(2))
    
    t_range(i, :) = [min(time_), max(time_)];
    
    [y_, ~] = hist(speed_mm_s, metpars.hbins_speed);
    
    lineH(1) = plot(metpars.hbins_speed, y_/sum(y_), ...
        'Color', 'b', ...
        'Parent', axH(3));
            
    hold(axH(3), 'on')

    [y_, ~] = hist(fv_mm_s, metpars.hbins_speed);
    
    lineH(2) = plot(metpars.hbins_speed, y_/sum(y_), ...
        'Color', 'r', ...
        'Parent', axH(3));
    
    [y_, ~] = hist(lv_mm_s, metpars.hbins_speed);
    
    lineH(3) = plot(metpars.hbins_speed, y_/sum(y_), ...
        'Color', 'm', ...
        'Parent', axH(3));
    
end


axH(1).Title.String = 'Distance traveled';
axH(1).XLabel.String = 'X distance (mm)';
axH(1).YLabel.String = 'Y distance (mm)';
axH(1).XLim = [-max(xy_range) max(xy_range)];
axH(1).YLim = [-max(xy_range) max(xy_range)];

axH(2).Title.String = 'Speed/velocity over time';
axH(2).XLabel.String = 'Time (s)';
axH(2).YLabel.String = 'instantaneous speed (mm/s)';
axH(2).YLim = [-70 70];
axH(2).XLim = [min(t_range(:, 1)) max(t_range(:, 2))];

axH(3).Title.String = 'Distribution of Speed/velocity';
axH(3).XLabel.String = 'speed (mm/s)';
axH(3).YLabel.String = 'Probability';

legend(axH(3), lineH, {'speed', 'f-velocity', 'l-velocity'}, 'Location', 'northeast')

savefig_int(figH, metpars.oDir, [fname, '_locomotion_stats'], ...
    [1 0 0 0 0 0 0 0 1])
close(figH)

end

function plot_vid_stim_results(wDat, fname, metpars)
% plot_vid_stim_results: plot basic stim vs locomotion stats
%
% Usage:
%   plot_vid_stim_results(wDat, fname, metpars)
%
% Args:
%   wDat: name pattern of directories to use
%       (default, [])
%   fname: figure name
%       (default, [])
%   metpars: parameters to update
%       (oDir: target directory to save figures)
%           (default, [pwd, filesep, 'fictrac'])
%       (ball_radious: depth of directory search)
%           (default, 4.5)
%       (hbins_speed: bins to use for histogram)
%           (default, [-20:0.1:60])
%       (framerate: frame rate)
%           (default, 100)
%       (timeres: temporal resolution in sec)
%           (default, 0.1)
%
% Notes
%   It resamples speed/velocity to metpars.timeres
%   this setup x = yaw, y = pitch, z = roll (this needs to be tested for every setup)
%   this setup y = forward and x = lateral (this needs to be tested for every setup)

% assumes each metadata only has one set of stimuli presented
stimuli_name = cellfun(@(x, y) [x '_AMP_' num2str(y)], ...
    wDat.sPars.name(1, :), chunk2cell(wDat.sPars.int(1, :), 1), ...
    'UniformOutput', false);

[sName, ~, stim_u_idx] = unique(stimuli_name, 'stable');
stim_all_idx = stim_u_idx(wDat.sPars.order(1, :));

if size(stim_all_idx, 1) == 1
    stim_all_idx = stim_all_idx';
end

stim_all_idx = stim_all_idx(1:size(wDat.sTime, 1), 1);

% plot basic motion stats
figH = figure('Position', genfigpos(1, 'center', [1200 900]));

% time to resample
time_ = wDat.vid.fstEn{1}(:, 1);   
time_res = ceil(time_(1)*10)/10:metpars.timeres:floor(time_(end)*10)/10;   
time_res = round(time_res*100)/100;

intmethod = 'linear';

pre_x_mm = interp1(time_, wDat.vid.var{1}(:, 15), ...
     time_res, intmethod);
pre_y_mm = interp1(time_, wDat.vid.var{1}(:, 16), ...
     time_res, intmethod);
pre_speed_mm_s = interp1(time_, wDat.vid.var{1}(:, 19), ...
     time_res, intmethod);
   
x_mm = pre_x_mm*metpars.ball_radious;
y_mm = pre_y_mm*metpars.ball_radious;
speed_mm_s = pre_speed_mm_s*metpars.ball_radious*metpars.framerate;
fv_mm_s = diff(y_mm)*metpars.framerate;
lv_mm_s = diff(x_mm)*metpars.framerate;

% round stimuli time
wDat.sTime = round(wDat.sTime*100)/100;

for ii = 1:size(wDat.sTime)

    idx2use = find(time_res >= metpars.time_buffer(1) + wDat.sTime(ii, 1) ...
        &  time_res <= metpars.time_buffer(2) + wDat.sTime(ii, 2));

    if ~isempty(idx2use)
        
        time_per_stim{ii, 1} = time_res(idx2use) - wDat.sTime(ii, 1);
        speed_mm_per_stim{ii, 1} = speed_mm_s(idx2use);
        fspeed_mm_per_stim{ii, 1} = fv_mm_s(idx2use + 1);
        lspeed_mm_per_stim{ii, 1} = lv_mm_s(idx2use + 1);
        
        % normalize to baseline
        index2use = time_per_stim{ii, 1} < 0;
        speed_mm_bs = nanmean(speed_mm_per_stim{ii, 1}(:, index2use), 2);
        fspeed_mm_bs = nanmean(fspeed_mm_per_stim{ii, 1}(:, index2use),2);
        lspeed_mm_bs = nanmean(lspeed_mm_per_stim{ii, 1}(:, index2use), 2);
        
        speed_mm_sd = nanstd(speed_mm_per_stim{ii, 1}(:, index2use), [], 2);
        fspeed_mm_sd = nanstd(fspeed_mm_per_stim{ii, 1}(:, index2use), [], 2);
        lspeed_mm_sd = nanstd(lspeed_mm_per_stim{ii, 1}(:, index2use), [], 2);
        
        % z-score to baseline
        speed_mm_per_stim_c{ii, 1} = speed_mm_per_stim{ii, 1} - speed_mm_bs;
        fspeed_mm_per_stim_c{ii, 1} = fspeed_mm_per_stim{ii, 1} - fspeed_mm_bs;
        lspeed_mm_per_stim_c{ii, 1} = lspeed_mm_per_stim{ii, 1} - lspeed_mm_bs;
        
        speed_mm_per_stim{ii, 1} = bsxfun(@rdivide, speed_mm_per_stim_c{ii, 1}, speed_mm_sd);
        fspeed_mm_per_stim{ii, 1} = bsxfun(@rdivide, fspeed_mm_per_stim_c{ii, 1}, fspeed_mm_sd);
        lspeed_mm_per_stim{ii, 1} = bsxfun(@rdivide, lspeed_mm_per_stim_c{ii, 1}, lspeed_mm_sd);       
        
    end
    
    clear idx2use

end

try 
    stim_all_idx = stim_all_idx(1:numel(time_per_stim));
    fprintf([fname, ' stim # ', ...
        num2str(numel(time_per_stim)), '\n'])
catch
    fprintf([fname, '\n'])
    return
    keyboard
end

for i = 1:numel(sName)
    
    axH(i*3 - 2) = subplot(numel(sName), 3, i*3 - 2);
    axH(i*3 -1) = subplot(numel(sName), 3, i*3 - 1);
    axH(i*3) = subplot(numel(sName), 3, i*3);

    stim2plot = find(stim_all_idx == i)';
    colorvect_1 = colorGradient([1 1 1], ...
        [0 0 1], numel(stim2plot) + 3);
    colorvect_1 = colorvect_1(4:end, :);
    colorvect_2 = colorGradient([1 1 1], ...
        [1 0 0], numel(stim2plot) + 3);
    colorvect_2 = colorvect_2(4:end, :);
    colorvect_3 = colorGradient([1 1 1], ...
        [1 0 1], numel(stim2plot) + 3);
    colorvect_3 = colorvect_3(4:end, :);
    
    colormap = parula(numel(stim2plot));
    
    for ii = stim2plot

        % plot speed/velocities
        plot(time_per_stim{ii, 1}, ...
            speed_mm_per_stim_c{ii, 1}, ...
            'Color', colormap(ii == stim2plot, :), ...
            'Parent', axH(i*3))

        hold(axH(i*3), 'on')

        plot(time_per_stim{ii, 1}, ...
            fspeed_mm_per_stim_c{ii, 1}, ...
            'Color', colormap(ii == stim2plot, :), ...
            'Parent', axH(i*3 - 1))

        hold(axH(i*3 - 1), 'on')

        plot(time_per_stim{ii, 1}, ...
            lspeed_mm_per_stim_c{ii, 1}, ...
            'Color', colormap(ii == stim2plot, :), ...
            'Parent', axH(i*3 - 2))

        hold(axH(i*3 - 2), 'on')

    end
    
    % overlay mean
    lineH = plot(time_per_stim{ii, 1}, ...
        mean(cell2mat(speed_mm_per_stim_c(stim2plot, 1)), 1), ...
        'Color', 'k', 'Linewidth', 1, ...
        'Parent', axH(i*3));
    
    plot(time_per_stim{ii, 1}, ...
        mean(cell2mat(fspeed_mm_per_stim_c(stim2plot, 1)), 1), ...
        'Color', 'k', 'Linewidth', 1, ...
        'Parent', axH(i*3 - 1))
    
    plot(time_per_stim{ii, 1}, ...
        mean(cell2mat(lspeed_mm_per_stim_c(stim2plot, 1)), 1), ...
        'Color', 'k', 'Linewidth', 1, ...
        'Parent', axH(i*3 - 2))
    
    axH(i*3 - 2).Title.String = {'l-velocity', strrep(sName{i}, '_', '-')};
    axH(i*3 - 1).Title.String = {'f-velocity', ''};
    axH(i*3).Title.String = {'speed', ''};
    axH(i*3 - 2).XLabel.String = 'Time (s)';
    axH(i*3 - 1).XLabel.String = 'Time (s)';
    axH(i*3).XLabel.String = 'Time (s)';
    axH(i*3 - 2).YLabel.String = 'centered speed (mm/s)';
    axH(i*3 - 2).XLim = [metpars.time_buffer(1) + wDat.sTime(1, 1), ...
        metpars.time_buffer(2) + wDat.sTime(1, 2)];

    legend(axH(i*3), lineH, 'mean')
    
end

savefig_int(figH, metpars.oDir, [fname, '_locomotion_vs_stim_centered'], ...
    [1 0 0 0 0 0 0 0 1])
close(figH)

figH_ = figure('Position', genfigpos(1, 'center', [1200 900]));

for i = 1:numel(sName)
    
    axH_(i*3 - 2) = subplot(numel(sName), 3, i*3 - 2);
    axH_(i*3 -1) = subplot(numel(sName), 3, i*3 - 1);
    axH_(i*3) = subplot(numel(sName), 3, i*3);

    stim2plot = find(stim_all_idx == i)';
    colorvect_1 = colorGradient([1 1 1], ...
        [0 0 1], numel(stim2plot) + 3);
    colorvect_1 = colorvect_1(4:end, :);
    colorvect_2 = colorGradient([1 1 1], ...
        [1 0 0], numel(stim2plot) + 3);
    colorvect_2 = colorvect_2(4:end, :);
    colorvect_3 = colorGradient([1 1 1], ...
        [1 0 1], numel(stim2plot) + 3);
    colorvect_3 = colorvect_3(4:end, :);
    
    colormap = parula(numel(stim2plot));
    
    for ii = stim2plot

        % plot speed/velocities
        plot(time_per_stim{ii, 1}, ...
            speed_mm_per_stim{ii, 1}, ...
            'Color', colormap(ii == stim2plot, :), ...
            'Parent', axH_(i*3))

        hold(axH_(i*3), 'on')

        plot(time_per_stim{ii, 1}, ...
            fspeed_mm_per_stim{ii, 1}, ...
            'Color', colormap(ii == stim2plot, :), ...
            'Parent', axH_(i*3 - 1))

        hold(axH_(i*3 - 1), 'on')

        plot(time_per_stim{ii, 1}, ...
            lspeed_mm_per_stim{ii, 1}, ...
            'Color', colormap(ii == stim2plot, :), ...
            'Parent', axH_(i*3 - 2))

        hold(axH_(i*3 - 2), 'on')

    end
    
    % overlay mean
    lineH = plot(time_per_stim{ii, 1}, ...
        mean(cell2mat(speed_mm_per_stim(stim2plot, 1)), 1), ...
        'Color', 'k', 'Linewidth', 1, ...
        'Parent', axH_(i*3));
    
    plot(time_per_stim{ii, 1}, ...
        mean(cell2mat(fspeed_mm_per_stim(stim2plot, 1)), 1), ...
        'Color', 'k', 'Linewidth', 1, ...
        'Parent', axH_(i*3 - 1))
    
    plot(time_per_stim{ii, 1}, ...
        mean(cell2mat(lspeed_mm_per_stim(stim2plot, 1)), 1), ...
        'Color', 'k', 'Linewidth', 1, ...
        'Parent', axH_(i*3 - 2))
    
    axH_(i*3 - 2).Title.String = {'l-velocity', strrep(sName{i}, '_', '-')};
    axH_(i*3 - 1).Title.String = {'f-velocity', ''};
    axH_(i*3).Title.String = {'speed', ''};
    axH_(i*3 - 2).XLabel.String = 'Time (s)';
    axH_(i*3 - 1).XLabel.String = 'Time (s)';
    axH_(i*3).XLabel.String = 'Time (s)';
    axH_(i*3 - 2).YLabel.String = 'z-scored speed (SD)';
    axH_(i*3 - 2).XLim = [metpars.time_buffer(1) + wDat.sTime(1, 1), ...
        metpars.time_buffer(2) + wDat.sTime(1, 2)];

    legend(axH_(i*3), lineH, 'mean')
    
end

savefig_int(figH_, metpars.oDir, [fname, '_locomotion_vs_stim_zscored'], ...
    [1 0 0 0 0 0 0 0 1])
close(figH_)

end
