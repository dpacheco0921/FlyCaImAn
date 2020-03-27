function batch_brainmaskgen(FolderName, FileName, iparams)
% batch_brainmaskgen: function to generate a binary mask of pixels that belong
%   to neural tissue (based on pixel fluorescence) to use for ROI segmentation.
%
% Usage:
%   batch_brainmaskgen(FolderName, FileName, iparams)
%
% Args:
%   FolderName: name of folders to load
%   FileName: name of files to load
%   iparams: parameters to update
%       (cDir: current directory)
%       (redo: redo gate)
%       (debug: debug gate)
%       (fsuffix: suffix of files to load)
%           ('_prosmetadata' or '_metadata')
%       (fo2reject: folders to reject)
%       (fi2reject: files to reject)
%       (tsuffix: name type)
%           (deaulf, 0 [day _ fly])
%           (1 [day _ fly _ rep])
%       %%%%%%%%%%%% editing of mean image before generating histograms %%%%%%%%%%%%
%       (minsize: minumun size of connected components when generating binary mask)
%           (default, 5^3)
%       (blurflag: flag to blur image)
%           (default, flag to blur image)
%       (manual: flag to require manual confirmation of binary mask)
%           (default, 1)
%       (all: flag to use all pixels, instead of mask)
%           (default, 0)
%
% Notes:
%   generates wDat.fths_zbin, wDat.fths_zval and wDat.min_f which are used by
%   batch_formatstacks or batch_stitch_format_stacks_*.m
%
% ToDo:

% Deafult params
pbm = [];
pbm.cDir = pwd;
pbm.redo = 0;
pbm.fsuffix = '_prosmetadata'; 
pbm.fo2reject = {'.', '..', 'preprocessed', 'BData'};
pbm.fi2reject = {'Zstack'};
pbm.tsuffix = 0;
pbm.minsize = 5^3;
pbm.blurflag = 1;
pbm.manual = 1;
pbm.all = 0;

% update variables
if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('FileName', 'var'); FileName = []; end
if ~exist('iparams', 'var'); iparams = []; end
pbm = loparam_updater(pbm, iparams);

% Selecting folders
f2run = dir;
f2run = str2match(FolderName, f2run);
f2run = str2rm(pbm.fo2reject, f2run);
f2run = {f2run.name};

fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']); 
    cd(f2run{i});
    runperfolder(FileName, pbm);
    cd(pbm.cDir)
    fprintf('\n')
    
end

fprintf('... Done\n')

end

function runperfolder(fname, pbm)
% runperfolder: run all files per folder
%
% Usage:
%   runperfolder(fname, pbm)
%
% Args:
%   fname: file name pattern
%   pbm: internal parameters structure

% Run all flies per folder
if pbm.tsuffix == 0
    nametype = 1;
else
    nametype = 2;
end

[f2plot, ~, ~] = ...
    rdir_namesplit(fname, '.mat', ...
    pbm.fsuffix, pbm.fi2reject, [], nametype);
f2plot = unique(f2plot);

fprintf(['Generating brain mask from ', ...
    num2str(numel(f2plot)),' flies\n'])

for fly_i = 1:numel(f2plot)
    processfunc(f2plot{fly_i}, pbm);
end

end

function processfunc(f2run, pbm)
% processfunc: for each filename compile 
%   all sub-stacks in the right order
%
% Usage:
%   processfunc(f2run, pbm)
%
% Args:
%   fname: file name
%   pbm: internal parameters structure

% Load metadata
load([f2run, '_', strrep(pbm.fsuffix, '_', ''), ...
    '.mat'], 'wDat')

if pbm.redo || (~isfield(wDat, 'bSide') ...
        && ~isfield(wDat, 'bMask')) ...
        && isfield(wDat, 'vSize')
    
    % Input brain side used
    if pbm.manual
        
        % plot z-projection
        htemp = figure();
        try
            imagesc(max(wDat.RedChaMean, [], 3));
        catch
            imagesc(max(wDat.GreenChaMean, [], 3));
        end
        
        % input brain side details
        fprintf(['Define hemisphere ', ...
            '(L (left), R(right), or C(center))'])
        fprintf('Define Dorso-ventral axis (d, dorsal, v: ventral) : \n')
        
        wDat.bSide = input('brain side : ');
        close(htemp)
        
        fprintf(['Brain side chosen : ', ...
            wDat.bSide, '\n'])
        
    else
        
        wDat.bSide = 'L_d';
        
    end
    
    if pbm.manual
        
        % smooth mean images (XYZ)
        bGreen = [];

        if pbm.blurflag
            bGreen = imblur(wDat.GreenChaMean, ...
                [2 2 2], [5 5 3], 3);
        end

        % plot fluorescence histograms per sub-segment
        figH = figure();
        axH(1) = subplot(2, 1, 1);
        axH(2) = subplot(2, 1, 2);
        wDat_F_hist_per_stack(wDat, 1, ...
            [], bGreen, figH, axH(1))

        % manually input the final intensity threshold
        lnbis = wDat_get_stackedges(wDat);
        wDat = fitgauss_to_fluohist(wDat, lnbis, ...
            'gauss5', [2 5], 1, bGreen, ...
            figH, axH(2), 1);
        
        pi.figpos = genfigpos(1, 'nw', [1156 510]);
        
        plot4Dproj(wDat.RedChaMean, [], wDat.vSize, pi); 
        plot4Dproj(wDat.GreenChaMean, [], wDat.vSize, pi)
        
        % always manually inspect mask
        edit brainmaskgen_maninput;
        keyboard;
        
        close all

        % generate binary mask of brain pixels
        wDat.bMask = wDat_generatemask(wDat, pbm.minsize, bGreen);

    else
        
        % use all pixels (overwrite previous masking)
        fprintf('Using all pixels\n');
        wDat.bF_ths = [];
        wDat.bF_zbin = [];
        wDat.bMask = true(size(wDat.GreenChaMean));
        
    end
    
    save([f2run, '_', strrep(pbm.fsuffix, '_', ''), '.mat'], ...
        'wDat', '-append')
    
elseif isfield(wDat, 'bSide') && ~isfield(wDat, 'bMask')
    
    fprintf('No imaging data found\n')
    
elseif isfield(wDat, 'vSize')

    fprintf('no \n')
    
end

end
