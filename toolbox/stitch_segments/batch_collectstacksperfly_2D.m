function batch_collectstacksperfly_2D(FolderName, FileName, iparams)
% batch_collectstacksperfly_2D: Second part of preprocessing for data
%   that does not undergo compiling
%   1) reads metadata generated by batch_SpaTemp_ResFilt_2D 
%   2) saves variables in a ROIseg-compatible way
%
% Usage:
%   batch_collectstacksperfly_a(FolderName, FileName, iparams)
%
% Args:
%   FolderName: name of folders to load
%   FileName: name of files to load
%   iparams: parameters to update
%       (cDir: current directory)
%       (fo2reject: folders to reject)
%       (fi2reject: files to reject)
%       (fsuffix: suffix of files to load)
%           (default, '_rawdata')
%       (oDir: temporary folder to copy data to)
%           (default, [])
%       %%%%%%%%%%%% shift fluorescence distribution %%%%%%%%%%%%
%       (bkgate: flag for background substraction)
%           (default, 0)
%       (blowcap: fluorescence below which it is zerored)
%           (default, 0)
%       (fshift: shift distribution of F to the positive side)
%           (default, 6)
%       %%%%%%%%%%%% parpool & server related %%%%%%%%%%%%
%       (serId: server id)
%           (default, 'int')
%       (corenum: number of cores)
%           (default, 4)
%       (idp_run_flag: flag to run each selected file independently)
%           (default, 0)
%       (green_field2save: flag to save Data variable as [Y Yr], Yr is Y in 2D)
%           (default, [1 1])
%       (red_field2save:  flag to save Data variable as [Y Yr], Yr is Y in 2D)
%           (default, [0 1])

cspf_2d = [];
cspf_2d.cDir = pwd;
cspf_2d.fo2reject = {'.', '..', 'preprocessed', 'BData'};
cspf_2d.fi2reject = {'Zstack'};
cspf_2d.fisuffix = '_rawdata';
cspf_2d.oDir = [];
cspf_2d.bkgate = 0;
cspf_2d.blowcap = 0;
cspf_2d.fshift = 6;
cspf_2d.serId = 'int';
cspf_2d.corenum = 4;
cspf_2d.idp_run_flag = 0;
cspf_2d.green_field2save = [1 1];
cspf_2d.red_field2save = [0 1];

if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('FileName', 'var'); FileName = []; end
if ~exist('iparams', 'var'); iparams = []; end
cspf_2d = loparam_updater(cspf_2d, iparams);

% start pararell pool if not ready yet
ppobj = setup_parpool(cspf_2d.serId, cspf_2d.corenum);

if ~isempty(cspf_2d.oDir)
    
    if ~exist(cspf_2d.oDir, 'dir')
        mkdir(cspf_2d.oDir);
    end    
    fprintf(['Copying output files at : ', ...
        strrep(cspf_2d.oDir, filesep, ' '), '\n'])
    
end

% Selecting folders
f2run = dir;
f2run = str2match(FolderName, f2run);
f2run = str2rm(cspf_2d.fo2reject, f2run);
f2run = {f2run.name};
fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']); 
    cd(f2run{i});
    runperfolder(FileName, cspf_2d);
    cd(cspf_2d.cDir);
    fprintf('\n')
    
end

delete_parpool(ppobj);

fprintf('... Done\n')

end

function runperfolder(fname, cspf_2d)
% runperfolder: run all files per folder
%
% Usage:
%   runperfolder(fname)
%
% Args:
%   fname: file name pattern
%   cspf_2d: parameter variable

% determine if fname narrows down to just one file
[f2plot, ~, rep2plot] = rdir_namesplit(...
    fname, '.mat', cspf_2d.fisuffix, ...
    cspf_2d.fi2reject, [], 1);
fprintf('Reformatting stacks from');

if cspf_2d.idp_run_flag == 1
    
    for i = 1:numel(f2plot)
        % run single file
        f2plot{i} = [f2plot{i}, '_', num2str(rep2plot(i))];
    end
    fprintf(' each seg independently ') 
    
else
    
    if numel(f2plot) == 1 && numel(rep2plot) == 1
        % run single file
        f2plot{1} = [f2plot{1}, '_', num2str(rep2plot)];
        fprintf(' one fly_seg ')
    else
        % Run all flies per folder
        f2plot = unique(f2plot);
    end
    
end

fprintf([num2str(numel(f2plot)), ' flies\n'])

for file_i = 1:numel(f2plot)
    
    % Compiling all stacks per fly
    [filename, ~, repnum] = ...
        rdir_namesplit(f2plot{file_i}, ...
        '.mat', cspf_2d.fisuffix, cspf_2d.fi2reject, [], 1);
    
    % run per fly
    for rep_i = 1:numel(repnum)
        
        if numel(filename) == 1
            
            fprintf(['Reformatting fly: ', ...
                filename{1}, '_', num2str(repnum(rep_i)),'\n'])
            fcompiler([filename{1}, '_', num2str(repnum(rep_i))], cspf_2d);
            
        else
            fprintf('error')
        end
        
    end
    
end

end

function fcompiler(fname, cspf_2d)
% fcompiler: for each filename compile all sub-stacks in the right order
%
% Usage:
%   fcompiler(fname)
%
% Args:
%   fname: file name
%   cspf_2d: parameter variable

% load previously generated wDat
load([fname, '_metadata.mat'], 'wDat');
wDat.cDir = pwd;

% add local folder
if isfield(wDat, 'min_f')
    wDat = rmfield(wDat, 'min_f');
end

% correct each channel separately
vList = whos('-file', [fname, '_rawdata']);

if ~contains([vList.name], 'Y')
    
    % saves Y and Yr
    savedata_per_cha(fname, 1, wDat, ...
        cspf_2d.green_field2save, ...
        cspf_2d.bkgate, cspf_2d.fshift, cspf_2d.blowcap);
    
    % only saves Yr for red channel
    try
        savedata_per_cha(fname, 2, wDat, ...
        cspf_2d.red_field2save, ...
            cspf_2d.bkgate, cspf_2d.fshift, cspf_2d.blowcap);
    end
    
end

% save metadata
dataObj = matfile([fname, '_rawdata.mat'], 'Writable', true);

% perform 2D neighcorr
wDat.lc2D = neighcorr_2D(dataObj);

% compile field
wDat.cspf = 1;
save([fname, '_metadata.mat'], 'wDat', '-append');

if ~isempty(cspf_2d.oDir)
    
    copyfile([fname, '_metadata.mat'], ...
        [cspf_2d.oDir, filesep, fname, '_metadata.mat']);
    copyfile([fname, '_rawdata.mat'], ...
        [cspf_2d.oDir, filesep, fname, '_rawdata.mat']);
    
    try
        copyfile([fname, '_refdata.mat'], ...
            [cspf_2d.oDir, filesep, fname, '_refdata.mat']);
    end
    
end

end

function savedata_per_cha(fname, cha2use, wDat, ...
    format2save, bkgate_, fshift_, blowcap_)

tinit = tic;

% Load files
load([fname, '_metadata.mat'], 'iDat');

if cha2use == 1
    
    % green channel
    load([fname, '_rawdata.mat'], 'Data');
    delete([fname, '_rawdata.mat'])
    dataObj = matfile([fname, '_rawdata.mat'], 'Writable', true);
    
else
    
    % red channel
    load([fname, '_refdata.mat'], 'Data');
    delete([fname, '_refdata.mat'])
    dataObj = matfile([fname, '_refdata.mat'], 'Writable', true);
    
end

% background substract F
if bkgate_
    if cha2use == 1
        Data = Data - iDat.bs(end); 
    else
        Data = Data - iDat.bs(1);
    end
end
Data = Data + fshift_;
Data(Data < blowcap_) = blowcap_;

% Prune Data
Data = pruneIm(Data, wDat.mask);

% saving Data
if cha2use == 1
    wDat.GreenTrend(1, :) = stacktrend(Data, wDat.bMask);
else
    wDat.RedTrend(1, :) = stacktrend(Data, wDat.bMask);
end

% saving data using relative indexes
if format2save(1)
    dataObj.Y(1:wDat.fSize(1), 1:wDat.fSize(2), 1:wDat.Tn) = Data;
end

% indexes are sequential from plane to plane
if format2save(2)
    pixelindex = 1:prod(wDat.vSize);
    dataObj.Yr(pixelindex, 1:wDat.Tn) = ...
        reshape(Data, [prod(wDat.vSize), wDat.Tn]);
end

dataObj.nY = min(Data(:));
dataObj.sizY = [wDat.vSize(1:2), wDat.Tn];
clear iDat Data

fprintf([num2str(toc(tinit)), ' seconds\n'])

end

function fmed = stacktrend(Y, mask)
% stacktrend: get median trend per stack for just neural tissue (using mask)
%
% Usage:
%   stacktrend(fname)
%
% Args:
%   Y: 2DxT image
%   mask: 2D mask

dDim = size(Y);
mask = mask(:);
Y = reshape(Y, [prod(dDim(1:2)), dDim(3)]);
Y = Y(mask ~= 0, :);
fmed = median(Y, 1);

end
