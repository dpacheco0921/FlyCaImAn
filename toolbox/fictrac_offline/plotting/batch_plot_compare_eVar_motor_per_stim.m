function batch_plot_compare_eVar_motor_per_stim(...
    ifiles_1, figname, iparams)
% batch_plot_compare_eVar_motor_per_stim: plot explained
%   variance by mean/median locomotor variable
%
% Usage:
%   batch_plot_compare_eVar_motor_per_stim(...
%       ifiles_1, figname, iparams)
%
% Args:
%   ifiles_1: name pattern of files to use
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
metpars.dir_depth = 0;
metpars.oDir = [pwd, filesep, 'fictrac_compare'];

if ~exist(metpars.oDir, 'dir')
   mkdir(metpars.oDir)
end

% update variables
if ~exist('iparams', 'var'); iparams = []; end
metpars = loparam_updater(metpars, iparams);

fprintf(['Generating plots for # ', ...
    num2str(length(ifiles_1)), ' pairs\n'])

for i = 1:numel(ifiles_1)
    wDat{i} = run_each(ifiles_1{i}, metpars);
end

% collect stimuli name
sName_ = cf(@(x) x.stim_name', wDat);
sName_ = cat(1, sName_{:});
stimuli_name = unique(sName_, 'stable');

plot_per_wDat(wDat, stimuli_name, ...
    metpars, figname)

fprintf('... Done\n')

end

function plot_per_wDat(wDat, stimuli_name, ...
    metpars, figurename)
% plot_per_wDat: plot average speed per stimuli
%
% Usage:
%   plot_per_wDat(wDat, stimuli_name, ...
%       metpars, figurename)
%
% Args:
%   wDat: file name
%   stimuli_name: input parameters

figH = figure('Position', genfigpos(1, 'center', [225 900]));
for i = 1:numel(stimuli_name)
    axH(i) = subplot(numel(stimuli_name), 1, i);
end

for i = 1:numel(stimuli_name)
    
    for j = 1:numel(wDat)
        
        stim2plot = contains(wDat{j}.stim_name, stimuli_name(i));
        trace2plot(j, :) = wDat{j}.eVar_mean_per_stim(1:7, stim2plot)*100;
        plot((1:7) + rand(1, 7)*0.1, trace2plot(j, :), ...
            '.', 'Color', [0.5 0.5 0.5], 'markersize', 20, 'Parent', axH(i));
        hold(axH(i), 'on')        
        
    end
    
    plot((1:7), mean(trace2plot, 1), ...
        '.k', 'markersize', 20, 'Parent', axH(i))
    clear trace2plot
    
end

for i = 1:numel(stimuli_name)
    axH(i).Title.String = stimuli_name{i};
    axH(i).YLabel.String = 'explained variance (%)';
    axH(i).YLim = [0 20];
end

axH(end).XTick = 1:1:7;
axH(end).XTickLabel = {'speed', 'fV', 'lV', ...
    'yaw', 'pitch', 'roll', 'SVD-1'};
axH(end).XTickLabelRotation = 45;

savefig_int(figH, metpars.oDir, ...
    [figurename, '_evar_motor_per_stim'], ...
    [1 0 0 0 0 0 0 0 1 0 1])

close(figH)

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
filename = str2match(ifile, f2run);
filename = {filename.name}';
[filename, iDir_1] = split_path(filename);
filename = strrep(filename, metpars.fmetsuffix, '');

if numel(filename) == 1
    
    load([iDir_1{1}, filesep, filename{1}, metpars.fmetsuffix], 'wDat', 'sti_ln_mot')
    wDat.stim_name = sti_ln_mot.stim_name;
    wDat.eVar_mean_per_stim = sti_ln_mot.eVar_mean_per_stim;
    wDat.eVar_med_per_stim = sti_ln_mot.eVar_med_per_stim;
    
else
    fprintf('more than one possible file\n')
end

end
