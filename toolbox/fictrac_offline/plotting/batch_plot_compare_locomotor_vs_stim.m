function batch_plot_compare_locomotor_vs_stim(...
    ifiles_1, ifiles_2, figname, iparams)
% batch_plot_compare_locomotor_vs_stim: plot stimulus 
%   triggered changes in locomotion
%
% Usage:
%   batch_compare_fictrac_results(...
%       ifiles_1, ifiles_2, figname, iparams)
%
% Args:
%   ifiles_1: name pattern of files to use
%   ifiles_2: name pattern of files to use
%   figname: figure name
%   iparams: parameters to update
%       (fmetsuffix: suffix of metadata file)
%           (default, '_metadata.mat')
%       (refcha: reference channel)
%           (default, 1)
%       (dir_depth: depth of directory search)
%           (default, 0)
%       (oDir: target directory to save figures)
%           (default, [pwd, filesep, 'fictrac'])
%       (ball_radious: depth of directory search)
%           (default, 4.5)
%       (sig: std of gaussian kernel to smooth motor variable)
%           (default, 3)
%       (siz: size of kernel  to smooth motor variable)
%           (default, 10)
%       (framerate: frame rate)
%           (default, 100)
%       (time_buffer: time before and after stimulus)
%           (default, [-5 +5])
% Notes

% default params
metpars = [];
metpars.fmetsuffix = '_metadata.mat';
metpars.fi2reject = 'Zstack';
metpars.refcha = 1;
metpars.dir_depth = 0;
metpars.oDir = [pwd, filesep, 'fictrac_compare'];
metpars.ball_radious = 4.5;
metpars.sig = 3;
metpars.siz = 10;
metpars.framerate = 100;
metpars.time_buffer = [-10 +10];

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

for i = 1:numel(ifiles_1)
    wDat{i} = run_each(ifiles_1{i}, metpars);
end

for i = 1:numel(ifiles_2)
    wDat_{i} = run_each(ifiles_2{i}, metpars);
end

% collect stimuli name
sName_ = cf(@(x) x.sName', wDat);
sName_ = cat(1, sName_{:});
stimuli_name = unique(sName_, 'stable');

wDat = homogenize_stim_name(wDat, stimuli_name);
wDat_ = homogenize_stim_name(wDat_, stimuli_name);

plot_per_wDat(wDat, stimuli_name, ...
    metpars, [figname, '_a'], ifiles_1)
plot_per_wDat(wDat_, stimuli_name, ...
    metpars, [figname, '_b'], ifiles_2)

fprintf('... Done\n')

end

function plot_per_wDat(wDat, stimuli_name, ...
    metpars, figurename, filenames)
% plot_per_wDat: plot average speed per stimuli
%
% Usage:
%   plot_per_wDat(wDat, stimuli_name, ...
%       metpars, figurename)
%
% Args:
%   wDat: file name
%   stimuli_name: input parameters

% setup figure-1
figH = figure('Position', genfigpos(1, 'center', [300 900]));

for i = 1:numel(stimuli_name)
    axH(i) = subplot(numel(stimuli_name), 1, i);
end

color_1 = colorGradient(rgb('Purple'), [1 1 1], 4);
colorvect = colorGradient([0.2 0.2 0.2], [1 1 1], numel(wDat) + 3);

for i = 1:numel(stimuli_name)

    k = 1;
    
    for j = 1:numel(wDat)
        
        stim2plot = wDat{j}.sIdx == i;
        
        if sum(stim2plot) ~= 0
            lineH(j) = plot(wDat{j}.time_per_stim, ...
                wDat{j}.speed_mm_per_stim{stim2plot}, ...
                'Color', colorvect(j, :), ...
                'Parent', axH(i));
            hold(axH(i), 'on')
            
            x_(k, :) = wDat{j}.speed_mm_per_stim{stim2plot};
            k = k + 1;
        end
                    
    end
    
    plot(wDat{1}.time_per_stim, ...
        median(x_, 1), 'Color', color_1(1, :), ...
        'Linewidth', 2, 'Parent', axH(i))
    
    clear x_
    
end

for i = 1:numel(stimuli_name)
    axH(i).Title.String = stimuli_name{i};
    axH(i).YLim = [-2 2];
    axH(i).XLabel.String = 'time (s)';
    axH(i).YLabel.String = 'speed (SD)';
end

legend(axH(1), lineH, strrep(filenames, '_', '-'), ...
    'Location', 'northeast')

savefig_int(figH, metpars.oDir, ...
    [figurename, '_locomotion_vs_stim'], ...
    [1 0 0 0 0 0 0 0 1 0 1])

close(figH)

end

function wDat = homogenize_stim_name(wDat, stimuli_name)
% homogenize_stim_name: function that homogenizes stimuli name across
%   different recordings
%
% Usage:
%   wDat = homogenize_stim_name(wDat, stimuli_name)
%
% Args:
%   ifile: file name
%   metpars: input parameters

for i = 1:numel(wDat)
    
    fprintf('*')
    [stim_name_u, ~, sti_idx] = ...
        unique(wDat{i}.sName(wDat{i}.sIdx), 'stable');
    [~, idx_lo] = ismember(stimuli_name, stim_name_u);
    idx_target = 1:numel(stimuli_name);

    % remove idx (sName) not found
    idx2del = find(idx_lo == 0);
    idx_lo(idx2del) = [];
    idx_target(idx2del) = [];
    sti_idx = changem(sti_idx, idx_target(:), idx_lo(:));

    % update stimuli name: sName
    wDat{i}.sName = stimuli_name;

    % update stimuli order: sIdx
    wDat{i}.sIdx = sti_idx;

    clear ib sti_idx idx_target idx_lo idx_lo

end

end

function wDat = run_each(ifile, metpars)
% run_each: function that runs each pair
%
% Usage:
%   wDat = run_per_pair(ifile, metpars)
%
% Args:
%   ifile: file name
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
filename_1 = str2match(ifile, f2run);
filename_1 = {filename_1.name}';
[filename_1, iDir_1] = split_path(filename_1);
filename_1 = strrep(filename_1, metpars.fmetsuffix, '');

if numel(filename_1) == 1
    wDat = load_gen_vars2plot([iDir_1{1}, ...
        filesep, filename_1{1}, metpars.fmetsuffix], metpars);
else
    fprintf('more than one possible file\n')
end

end

function wDat = load_gen_vars2plot(...
    filename_full, metpars)
% load_gen_vars2plot: load wDat, generate 
%   stimulus name and indeces and all variables to plot
%
% Usage:
%   wDat = load_gen_vars2plot(...
%      filename_full, metpars)
%
% Args:
%   name_full: full path of file to load
%   metpars: internal parameters

load(filename_full, 'wDat');

if ~exist('wDat', 'var') || ...
        ~isfield(wDat, 'vid') || isempty(wDat.vid.var)
    return
end

% generate stimulus name and indeces
stimuli_name = cellfun(@(x, y) [x '_AMP_' num2str(y)], ...
    wDat.sPars.name(1, :), chunk2cell(wDat.sPars.int(1, :), 1), ...
    'UniformOutput', false);

[sName, ~, stim_u_idx] = unique(stimuli_name, 'stable');
stim_all_idx = stim_u_idx(wDat.sPars.order(1, :));

if size(stim_all_idx, 1) == 1
    stim_all_idx = stim_all_idx';
end

stim_all_idx = stim_all_idx(1:size(wDat.sTime, 1), 1);

wDat.sIdx = stim_all_idx;
wDat.sName = sName;

% time to resample
time_ = wDat.vid.fstEn{1}(:, 1);   
time_res = ceil(time_(1)*10)/10:0.1:floor(time_(end)*10)/10;   
time_res = round(time_res*100)/100;
intmethod = 'linear';

[Y, ~, ~] = load_edit_fictrac_motor_var(wDat, metpars.ball_radious, ...
    {'speed', 'fV', 'lV'}, metpars.sig, metpars.siz, 0);

Y = interp1(time_, Y', time_res, intmethod);
Y = framegapfill(find(isnan(Y(:, 1))), Y');

% round stimuli time
wDat.sTime = round(wDat.sTime*100)/100;

% get traces per stimuli and do baseline substraction
for ii = 1:size(wDat.sTime)

    idx2use = find(time_res >= metpars.time_buffer(1) + wDat.sTime(ii, 1) ...
        &  time_res <= metpars.time_buffer(2) + wDat.sTime(ii, 2));

    if ~isempty(idx2use)
        time_per_stim{ii, 1} = time_res(idx2use) - wDat.sTime(ii, 1);
        
        speed_mm_per_stim{ii, 1} = Y(1, idx2use);
        fspeed_mm_per_stim{ii, 1} = Y(2, idx2use);
        lspeed_mm_per_stim{ii, 1} = Y(3, idx2use);
        
        speed_mm_per_stim{ii, 1} = (speed_mm_per_stim{ii, 1} - ...
            nanmean(speed_mm_per_stim{ii, 1}(time_per_stim{ii, 1} < 0)))/...
            std(speed_mm_per_stim{ii, 1}(time_per_stim{ii, 1} < 0));
        fspeed_mm_per_stim{ii, 1} = (fspeed_mm_per_stim{ii, 1} - ...
            nanmean(fspeed_mm_per_stim{ii, 1}(time_per_stim{ii, 1} < 0)))/...
            std(fspeed_mm_per_stim{ii, 1}(time_per_stim{ii, 1} < 0));
        lspeed_mm_per_stim{ii, 1} = (lspeed_mm_per_stim{ii, 1} - ...
            nanmean(lspeed_mm_per_stim{ii, 1}(time_per_stim{ii, 1} < 0)))/...
            std(lspeed_mm_per_stim{ii, 1}(time_per_stim{ii, 1} < 0));

    end
    
    clear idx2use

end

wDat.time_per_stim = time_per_stim{1};

% reduce stim time to match video frames
wDat.sTime = wDat.sTime(1:numel(time_per_stim), :);
wDat.sIdx = wDat.sIdx(1:numel(time_per_stim));

% calculate mean per stimulus
idx_unique = unique(wDat.sIdx);

for i = 1:numel(idx_unique)
    
    trial2avg = idx_unique(i) == wDat.sIdx;
    
    wDat.speed_mm_per_stim{i} = ...
        nanmedian(cell2mat(speed_mm_per_stim(trial2avg, 1)), 1);
    wDat.fspeed_mm_per_stim{i} = ...
        nanmedian(cell2mat(fspeed_mm_per_stim(trial2avg, 1)), 1);
    wDat.lspeed_mm_per_stim{i} = ...
        nanmedian(cell2mat(lspeed_mm_per_stim(trial2avg, 1)), 1);
    
end

wDat.sIdx = idx_unique;

end
