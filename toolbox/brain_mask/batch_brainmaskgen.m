function batch_brainmaskgen(FolderName, FileName, iparams)
% batch_brainmaskgen(FolderName, FileName, iparams)
% This function basically reads all flyfiles ('*prosmetadata') and displays
% distribution of the F per pixel and their gaussian fits so to manually decide on thresholds
% for mask and min_f
% generates wDat.fths_zbin, wDat.fths_zval and wDat.min_f which is used by
% batch_stackperfly_b

% Deafult params
global pbm

pbm = [];
pbm.cDir = pwd;
pbm.redo = 0;
pbm.fo2reject = {'.', '..', 'preprocessed', 'BData'};
pbm.fi2reject = {'Zstack'};
pbm.FolderName = [];
pbm.FileName = [];
pbm.fsuffix = '_prosmetadata'; 
pbm.tsuffix = 0; % default name type 0 [day _ fly], if 1 [day _ fly _ rep]
% editing of mean image before generating histograms
pbm.minsize = 5^3; % delete connected components that are smaller than this
pbm.blurgate = 1; % blur image
pbm.manual = 1; % require manual confirmation that masks are ok
pbm.all = 0; % use all pixels

if exist('FolderName','var'); pbm.FolderName = FolderName; end
if exist('FileName','var'); pbm.FileName = FileName; end
if ~exist('iparams', 'var'); iparams = []; end
pbm = loparam_updater(pbm, iparams);

% Selecting folders
f2run = dir;
f2run = GS_str2match(pbm.FolderName, f2run);
f2run = GS_str2rm(pbm.fo2reject, f2run);
f2run = {f2run.name};

fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']); 
    cd(f2run{i});
    runperfolder(pbm.FileName, pbm);
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
% processfunc: for each filename compile all sub-stacks in the right order
%
% Usage:
%   processfunc(f2run, pbm)
%
% Args:
%   fname: file name
%   cspf: internal parameters structure

% Load metadata
load([f2run, '_', strrep(pbm.fsuffix, '_', ''), '.mat'], 'wDat')

if pbm.redo || (~isfield(wDat, 'bSide') && ~isfield(wDat, 'bMask'))
    
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

        if pbm.blurgate
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
            pbm.blurgate, figH, ...
            axH(2), 1);
        
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
    
else
    
    fprintf('bMask already generated\n')
    
end

end
