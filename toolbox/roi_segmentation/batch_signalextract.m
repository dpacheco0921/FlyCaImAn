function batch_signalextract(FolderName, FileName, iparams)
% batch_signalextract: re-extract the ROI signals (raw, filtered, dfof)
%   from existing ROI segmentation results (*_prosroi.mat)
%
% Usage:
%   batch_signalextract(FolderName, FileName, iparams)
%
% Args:
%   FolderName: name of folders to run
%       (default, [], all folders)
%   FileName: name pattern of files to run
%       (default, [], all files)
%   iparams: parameters to update (fields of the global 'p' below), e.g.
%       (redo: 1 to extract again even if signals already exist)
%           (default, 0)
%       (fisuffix, osuffix: suffix of the ROI file, '_prosroi')
%       (orawfsuffix: suffix of the green channel data, '_prosdata')
%       (orefsuffix: suffix of the reference (red) channel, '_prosref')
%       (serverid, corenum: parpool settings, 'int' and 1)
%
% Returns:
%   nothing, it appends the updated 'roi' variable to each
%   *_prosroi.mat file
%
% Notes:
%   uses the global variable p to share settings with its sub-functions
global p
p = [];
p.cDir = pwd;
p.sep = filesep;
p.fo2reject = {'.', '..', 'preprocessed', 'BData'};
p.fi2reject = {'Zstack'};
p.FolderName = []; 
p.FileName = [];
p.corenum = 1; % numbre of cores to request
p.fisuffix = '_prosroi'; 
p.orawfsuffix = '_prosdata';
p.orefsuffix = '_prosref';
p.osuffix = '_prosroi';
p.redo = 0;
p.serverid = 'int';
p.corenum = 1;

if ~exist('iparams', 'var'); iparams = []; end
p = loparam_updater(p, iparams);

if exist('FolderName','var'); p.FolderName = FolderName; end
if exist('FileName','var'); p.FileName = FileName; end

% set up parpool
setup_parpool(p.serverid, p.corenum);

% Selecting folders
f2run = dir;
f2run = str2match(p.FolderName, f2run);
f2run = str2rm(p.fo2reject, f2run);
f2run = {f2run.name};
fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']); 
    cd(f2run{i})
    runperfolder(p.FileName)
    cd(p.cDir)
    fprintf('\n')
    
end

delete(gcp)
fprintf('... Done\n')

end

function runperfolder(fname)

global p

% Run all flies per folder
f2run = rdir(['.', p.sep, '*', p.fisuffix, '.mat']);
f2run = str2rm(p.fi2reject, f2run);
f2run = str2match(fname, f2run);
f2run = {f2run.name};
f2run = strrep(f2run, ['.', p.sep], '');
f2run = strrep(f2run, [p.fisuffix, '.mat'], '');
f2run = unique(f2run);

if isempty(f2run)
    warning('batch_signalextract: no %s files found in %s', [p.fisuffix, '.mat'], pwd)
end

fprintf(['Extracting signal from ', num2str(numel(f2run)),' flies\n'])

for fly_i = 1:numel(f2run)
    fprintf(['Extracting signal from fly: ', f2run{fly_i} ,'\n'])
    batch_getCatrace(f2run{fly_i});
end

end

function batch_getCatrace(filename)
% batch_getCatrace: extract the signals of the ROIs already segmented in
%   filename_prosroi.mat from the raw data and save them back in that file
%   (it does not segment ROIs again)

% parameters
global p
obj = CaROISeg;
% load roi data
load([filename, p.osuffix], 'roi')
% populate fields
populate_obj(obj, roi)
% add Green-Data
obj.data = matfile([filename, p.orawfsuffix, '.mat'], 'Writable', true);
userparams = roi.userparams;

% run signal extraction
if ~isfield(roi, 'filtered_') || ~isfield(roi, 'filtered') ...
        || isempty(roi.filtered_) || p.redo
    
    % extract signal from raw data
    rawsignal(obj, [filename, p.orefsuffix])
    % extract signal from raw and ref data
    %rawsignal_red(obj, [filename, p.orefsuffix])
    bas_estimate(obj)
    roi = snapshot(obj);
    roi.userparams = userparams;
    save([filename, p.osuffix], 'roi', '-append')
    
else
    
    fprintf(['Already extracted signal of fly: ', filename, '\n']);
    display([pwd, p.sep, filename, p.osuffix, '.mat'])
    
end

end
