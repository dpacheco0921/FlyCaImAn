function batch_plot_roiseg_results(Filename, oDir, iparams)
% batch_plot_roiseg_results: plots results from ROI segmentation
%
% Usage:
%   batch_plot_roiseg_results(FileName, oDir, iparams)
%
% Args:
%   Filename: name pattern of files to use
%       (default, [])
%   oDir: output directory.
%   iparams: parameters to update
%       (fsuffix: suffix of roi file - contains roi var)
%           (default, '_rawdata.mat')
%       (fmetsuffix: suffix of metadata file - contains wDat var)
%           (default, '_metadata.mat')
%       (dir_depth: depth of directory search)
%           (default, 0)
%       (roi2sel: indeces of ROIs to use)
%           (default, [])
%       (channel2use: define the channel to use)
%           (default, 'wDat.GreenChaMean')
%       (com_flag: flag to overlay ROI number using its center of mass)
%           (default, 0)
%       (indroi_color_flag: flag to use different color for each ROI)
%           (default, 0)
%
% Notes:
%   when using indroi_color_flag, be aware that for this you are generating
%       a full volume per ROI to be plotted, so it is fine for plotting tens of
%       ROIs but not for hundreds (but this depends on how big your volume is)

% default params
motpar = [];
motpar.fsuffix = '_prosroi.mat';
motpar.fmetsuffix = '_prosmetadata.mat';
motpar.dir_depth = 0;
motpar.roi2sel = [];
motpar.cha2use = 'wDat.GreenChaMean';
motpar.com_flag = 0;
motpar.indroi_color_flag = 0;
    
if ~exist('Filename', 'var')
    Filename = [];
end

if ~exist('oDir', 'var') || isempty(oDir)
    oDir = [pwd, filesep, 'smodrel'];
end

% update variables
if ~exist('iparams', 'var'); iparams = []; end
motpar = loparam_updater(motpar, iparams);

if ~isempty(oDir) && ~exist('oDir', 'var')
    mkdir(oDir)
end

% define files to use
if motpar.dir_depth == 0
    f2run = rdir(['.', filesep, '*', ...
        motpar.fsuffix]);
elseif motpar.dir_depth == 1
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', motpar.fsuffix]);
else
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', filesep, '*', ...
        motpar.fsuffix]);
end

f2run = {f2run.name}';
[filename, iDir] = split_path(f2run);
filename = strrep(filename, motpar.fsuffix, '');

if ~isempty(Filename)
    f2run = find(contains(filename, Filename));
    filename = filename(f2run);
    iDir = iDir(f2run);
end

fprintf(['Generating plots for ', ...
    num2str(numel(filename)), ' files\n'])

% plot ROI segmentation results
for i = 1:numel(filename)

    plot_roi_coverage_int([iDir{i}, filesep, filename{i}], ...
        oDir, motpar.fmetsuffix, motpar.roi2sel, ...
        motpar.cha2use, motpar.com_flag, motpar.indroi_color_flag);

end

fprintf('... Done\n')

end

function plot_roi_coverage_int(...
    filename, oDir, filesuffix, roi2sel, ...
    cha2use, com_flag, indroi_color_flag)
% plot_roi_coverage_int: plot roi coverage
%   it plots both sum of max normalized weights and binarized voxels
%
% Usage:
%   plot_roi_coverage_int(...
%       filename, oDir, filesuffix, roi2sel, ...
%       cha2use, com_flag, indroi_color_flag)
%
% Args:
%   filename: file name
%   oDir: output directory
%   filesuffix: metadata suffix
%   roi2sel: indeces of ROIs to use
%   channel2use: define the channel to use
%   com_flag: flag to overlay ROI number using its center of mass
%   indroi_color_flag: flag to use different color for each ROI

% load metadata file
load([filename, filesuffix], 'wDat')

% load ROI file
load([filename, '_prosroi'], 'roi')

% plot videos
plot_roi_coverage(filename, ...
    [1 1 1], wDat, roi, oDir, roi2sel, ...
    cha2use, com_flag, indroi_color_flag)
% plot_roi_coverage(filename, ...
%     vid2plot, wDat, roi, oDir, roi2sel, ...
%     cha2use, com_flag, indroi_color_flag)

end
