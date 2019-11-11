function batch_plot_brainside_MIP(...
    Filename, oDir, iparams)
% batch_plot_brainside_MIP: plot MIP of segment to 
%   see match with brain side
%
% Usage:
%   batch_plot_brainside_MIP(...
%       Filename, oDir, iparams)
%
% Args:
%   Filename: name pattern of files to use
%       (default, [])
%   oDir: output directory.
%   iparams: parameters to update
%       (fmetsuffix: suffix of metadata file)
%           (default, '_metadata.mat')
%       (dir_depth: depth of directory search)
%           (default, 0)
% Notes

% default params
metpars = [];
metpars.fmetsuffix = '_prosmetadata.mat';
metpars.dir_depth = 0;

if ~exist('Filename', 'var')
    Filename = [];
end

if isempty(oDir)
    oDir = [pwd, filesep, 'MIP'];
end

% update variables
if ~exist('iparams', 'var'); iparams = []; end
metpars = loparam_updater(metpars, iparams);

if ~isempty(oDir) && ~exist(oDir, 'dir')
    mkdir(oDir)
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

f2run = {f2run.name}';
[filename, iDir] = split_path(f2run);
filename = strrep(filename, metpars.fmetsuffix, '');

if ~isempty(Filename)
    f2run = find(contains(filename, Filename));
    filename = filename(f2run);
    iDir = iDir(f2run);
end

fprintf(['Generating plots for ', ...
    num2str(numel(filename)), ' files\n'])

% plot sMod results
for i = 1:numel(filename)

    load([iDir{i}, filesep, filename{i}, ...
        metpars.fmetsuffix], 'wDat');

    plot_bside_MIP(wDat, filename{i}, ...
        oDir)

end

fprintf('... Done\n')

end

function plot_bside_MIP(wDat, filename, ...
        oDir)
% plot_bside_MIP: plot MIP of red and green channel
%
% Usage:
%   plot_bside_MIP(wDat, filename, ...
%       oDir)
%
% Args:
%   wDat: input structure
%   filename: figure name
%   oDir: output directory to save figures or videos
    
figH = figure('Position', genfigpos(1, 'center', [600 400])); 
axH(1) = subplot(1, 2, 1);
axH(2) = subplot(1, 2, 2);

try imagesc(max(wDat.RedChaMean, [], 3), 'Parent', axH(1)); end
try imagesc(max(wDat.GreenChaMean, [], 3), 'Parent', axH(2)); end

axH(1).Title.String = strrep(wDat.bSide, '_', '-');
axH(2).Title.String = strrep(wDat.bSide, '_', '-');

figEdit(axH, figH, [], [], [], 'off')
savefig_int(figH, oDir, [filename, '_MIP'], ...
    [0 0 0 0 0 0 0 0 1])
close(figH)

end
