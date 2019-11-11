function batch_fill_frames_vols(FolderName, FileName, iparams)
% batch_fill_frames_vols: batch file to fill 
%   frames or volumes of specific files using framegapfill
%
% Usage:
%   Y = batch_fill_frames_vols(FolderName, FileName, iparams)
%
% Args:
%   FolderName: folders to use
%   FileName: files to use
%   iparams: input parameters
%       (ftype, type of replacing ('mean', default, or other (direct replacement))
%       (fside: side to use (right or left to the timepoint to replace))
%
% Notes
% see framegapfill

ipars = [];
ipars.ftype = 'mean';
ipars.fside = 'r';

if ~exist('FolderName','var')
    FolderName = [];
end

if ~exist('FileName','var')
    FileName = [];
end

if ~exist('iparams', 'var')
    iparams = [];
end

ipars = loparam_updater(ipars, iparams);

% Selecting folders
cDir = pwd;
f2reject = {'.', '..', 'preprocessed', 'BData'};
f2run = dir;
f2run = str2match(FolderName, f2run);
f2run = str2rm(f2reject, f2run);
f2run = {f2run.name};
fprintf('Filling missing Data\n')
fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']); 
    cd(f2run{i});
    runperfolder(FileName, ipars);
    cd(cDir)
    
end

fprintf('... Done\n')

end

function runperfolder(fname, ipars)
% runperfolder: run all files per folder
%
% Usage:
%   runperfolder(fname, ipars)
%
% Args:
%   fname: file name pattern
%   ipars: parameter variable

if ~exist('fname','var') || ...
        isempty(fname)
    fname = [];
end

f2fill = rdir(['.', filesep, '*_rawdata.mat']);

if ~isempty(fname)
    
    if iscell(fname)
        
        for i = 1:numel(fname)
            fname{i} = [fname{i}, '_rawdata.mat'];
        end
        
    else
        
        fname = [fname, '_rawdata.mat'];
        
    end
    
end

f2fill = str2match(fname, f2fill);
f2fill = {f2fill.name};
fprintf(['Filling n-flies : ', ...
    num2str(numel(f2fill)), '\n'])

for fly_i = 1:numel(f2fill)
    
    fprintf(['Filling file: ', ...
        strrep(f2fill{fly_i}, '\', '*'), '\n'])
    runperfile(f2fill{fly_i}, ipars)
    
end

end

function runperfile(f2fill, ipars)
% runperfile: fill volumes
%
% Usage:
%   runperfile(f2fill, ipars)
%
% Args:
%   f2fill: file name
%   ipars: parameter variable

% Load Data
load(f2fill, 'Data')
load(strrep(f2fill, '_rawdata', '_metadata'), 'iDat')

% Show all planes
figH = figure('Name', strrep(f2fill, '_', ' '));

aXh(1) = subplot(3, 2, 1);
aXh(2) = subplot(3, 2, 2);
aXh(3) = subplot(3, 2, 3);
aXh(4) = subplot(3, 2, 4);
aXh(5) = subplot(3, 2, 5);
aXh(6) = subplot(3, 2, 6);

Y2display = squeeze(max(max(Data, [], 1), [], 2));

imagesc(Y2display(:, :, 1), 'Parent', aXh(1))
title(aXh(1), 'Green')

if size(Y2display, 3) > 1
    
    imagesc(Y2display(:, :, 2), 'Parent', aXh(2))
    title(aXh(2), 'Red')
    
end

% Decide which planes / volumes to delete
framestdel = input('Frames | Volumes to fill : ');

% update iDat.DelFrames
iDat.DelFrames = framestdel;

% plot frames to delete
imagesc(Y2display(:, framestdel, 1), 'Parent', aXh(3))
title(aXh(3), 'Green2Fill')

if size(Y2display, 3) > 1
    
    imagesc(Y2display(:, framestdel, 2), 'Parent', aXh(4))
    title(aXh(4), 'Red2Fill')
    
end

% Fill and save Data
Data = framegapfill(iDat.DelFrames, ...
    Data, ipars.ftype, ipars.fside);

Y2display = squeeze(max(max(Data, [], 1), [], 2));

% plot new Data
imagesc(Y2display(:, :, 1), 'Parent', aXh(5))
title(aXh(5), 'GreenFilled')

if size(Y2display, 3) > 1
    
    imagesc(Y2display(:, :, 2), 'Parent', aXh(6))
    title(aXh(6), 'RedFilled')
    
end

fprintf('Saving ... ')
save(f2fill, 'Data', '-v7.3')
save(strrep(f2fill, '_rawdata', '_metadata'), ...
    'iDat', '-append')
fprintf('Done\n')

end
