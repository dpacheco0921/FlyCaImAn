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
%       (fsuffix: suffix of raw data)
%           (default, '_rawdata.mat')
%       (fmetsuffix: suffix of metadata file)
%           (default, '_metadata.mat')
%       (hbins: bins for histogram)
%           (default, -1:0.01:1)
%       (hbinsp: bins for histogram pvals)
%           (default, 0:0.01:1)
%       (prct2use: bins for histogram)
%           (default, 30)
%       (fdr: bins for histogram)
%           (default, .01)
%       (mccor_method: multiple comparison 
%           correction to use: dep, pdep, bh)
%           see calculate_pval.m
%            (default, 'dep')
%       (dir_depth: depth of directory search)
%           (default, 0)
%
% Notes:
%
% ToDo:
%   add legends to pull plots and stitch ones

% default params
motpar = [];
motpar.fsuffix = '_prosroi.mat';
motpar.fmetsuffix = '_prosmetadata.mat';
motpar.hbins = -1:0.01:1;
motpar.hbinsp = 0:0.01:1.1;
motpar.prct2use = 30;
motpar.fdr = 0.01;
motpar.mccor_method = 'dep';
motpar.dir_depth = 0;

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

% plot sMod results
for i = 1:numel(filename)

    plotperfly(filename{i}, iDir{i}, oDir, motpar);

end

fprintf('... Done\n')

end

function plot_roi_coverage_int(...
    filename, oDir, filesuffix)
% plot_roi_coverage_int: plot roi coverage
%   it plot both sum of max norm weights and binarized voxels
%
% Usage:
%   plot_roi_coverage_int(...
%       filename, oDir, filesuffix)
%
% Args:
%   filename: file name
%   oDir: output directory
%   filesuffix: metadata suffix

% load metadata file
load([filename, filesuffix], 'wDat')
% load ROI file
load([filename, '_prosroi'], 'roi')

% plot videos
plot_roi_coverage(filename, ...
    [1 1 1], wDat, roi, oDir)

end
