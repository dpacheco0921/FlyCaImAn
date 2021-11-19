function batch_plot_man_ROI_centers(Filename)
% batch_plot_man_ROI_centers: tool to display 3D matrices, 
%   Y could be a 3D matrix or a flatten arragement of a 3D matrix
%
% Usage:
%   slice3Dmatrix(Y, iparams)
%
% Args:
%   Filename: files to load
%   iparams: parameters to update


% default params
pi.figpos = genfigpos(1, 'center', [900 900]);
pi.cmap = parula;
pi.Y3text = [];
pi.Y3textcmap = 'cyan';
pi.Y3fontsiz = 12;
pi.vquality = 100;
pi.frate = 10;
pi.metadata_suffix = '_prosmetadata.mat';
pi.dir_depth = 0;
pi.oDir = [];
pi.cha2use = 'RedChaMean';

if ~exist('Filename', 'var')
    Filename = [];
end

if ~exist('oDir', 'var') || isempty(pi.oDir)
    pi.oDir = [pwd, filesep, 'man_roi'];
end

if ~exist(pi.oDir, 'dir')
   mkdir(pi.oDir)
end

% define files to use
if pi.dir_depth == 0
    f2run = rdir(['.', filesep, '*', ...
        pi.metadata_suffix]);
elseif pi.dir_depth == 1
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', pi.metadata_suffix]);
else
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', filesep, '*', ...
        pi.metadata_suffix]);
end

f2run = {f2run.name}';
[filename, iDir] = split_path(f2run);

if ~isempty(Filename)
    f2run = find(contains(filename, Filename));
    filename = filename(f2run);
    iDir = iDir(f2run);
end

filename = strrep(filename, pi.metadata_suffix, '');

fprintf(['Generating plots for ', ...
    num2str(numel(filename)), ' files\n'])

% plot manual segmentation results
for i = 1:numel(filename)
    
    runperfile(filename{i}, iDir{i}, pi)
    
end

end

function runperfile(filename, iDir, pi)

load([iDir, filesep, filename, pi.metadata_suffix], 'wDat')

eval(['im2use = wDat.', pi.cha2use, ';'])
figH = figure();
set(figH, 'color', 'w', ...
    'Position', pi.figpos);
axH = subplot(1, 1, 1);
colormap(axH, pi.cmap)

vHandle = [];

vHandle = VideoWriter([pi.oDir, filesep, filename, '_0']);
vHandle.Quality = pi.vquality;
vHandle.FrameRate = pi.frate;
open(vHandle);

for z = 1:size(im2use, 3)
    
    try
        plot_plane(im2use(:, :, z), wDat.manROI, figH, axH, z, pi)
        writeVideo(vHandle, getframe(axH));
    catch
        keyboard
    end
end

close(vHandle);

vHandle = [];

vHandle = VideoWriter([pi.oDir, filesep, filename, '_1']);
vHandle.Quality = pi.vquality;
vHandle.FrameRate = pi.frate;
open(vHandle);

% get centers from 3D image
new_centers = find(wDat.ROI_center_matrix(:) ~= 0);
[new_centers(:, 2) new_centers(:, 1) new_centers(:, 3)] = ...
    ind2sub(wDat.vSize, new_centers);
wDat.manROI.roi_center = new_centers;

for z = 1:size(im2use, 3)
    
    try
        plot_plane(im2use(:, :, z), wDat.manROI, figH, axH, z, pi)
        writeVideo(vHandle, getframe(axH));
    catch
        keyboard
    end
end

close(vHandle);

close(gcf)

end

function plot_plane(im2use, manROI, figH, axH, z, pi)
% function that plots each plane

% plot plane
imagesc(im2use, 'parent', axH)
axH.XTick = [];
axH.YTick = [];

% show channel name
axH.Title.String = [' z-' num2str(z)];

% plot ROIs in this plane
if ~isempty(manROI.roi_center)
    
    roi2use = find(manROI.roi_center(:, 3) == z);
    
    roi_center = manROI.roi_center(roi2use, :);
    roi_idx = manROI.roi_idx(roi2use, :);
    roi_radious = manROI.roi_radious(roi2use, :);
    
    for i = 1:numel(roi_idx)
        
        text(roi_center(i, 1), roi_center(i, 2), ...
            num2str(roi_idx(i)), 'color', pi.Y3textcmap,...
            'FontSize', max(roi_radious(i)/4, 14), 'Parent', axH)

    end
    
end

end
