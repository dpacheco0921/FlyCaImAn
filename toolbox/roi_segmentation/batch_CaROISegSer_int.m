function batch_CaROISegSer_int(FolderName, ...
    FileName, iparams, roiparams)
% batch_CaROISegSer_int: runs ROI segmentation per mat file, 
%   it requires the mat file to be in the right mat format
%
% Usage:
%   batch_CaROISegSer_int(FolderName, ...
%       FileName, iparams, roiparams)
%
% Args:
%   FolderName: name of folders to load
%   FileName: name of files to load
%   iparams: roi file selection and job related params
%   roiparams: roi segmentation paramters (see roiseg_params)
% 
% Notes:
% unlike ROI segmentation code for big stacks (CaROISegSer),
%   this one does jobpart-1 and -2 in the same run

global roipars
roipars = [];
roipars.cDir = pwd;
roipars.redo = 0;
% params dir related
roipars.fo2reject = {'.', '..', 'preprocessed', 'BData'};
roipars.fi2reject = {'Zstack'};
% metadata suffix
roipars.fimetsuffix = '_metadata';
% suffix of input data (it could be raw or reformated to (Y, Yr, sizY, nY))
roipars.fisuffix = '_rawdata';
roipars.fosuffix = '_prosdata';
% hardcoded suffix of data in the right format (Y, Yr, sizY, nY)
roipars.rfisuffix = '_prosref';
roipars.rfosuffix = '_prosref'; 
roipars.c_ths = []; % threshold used for C during segmentation
% params parpool related
roipars.serId = 'int';
roipars.corenum = 1;

if ~exist('iparams', 'var'); iparams = []; end
roipars = loparam_updater(roipars, iparams);

if ~exist('FolderName','var'); FolderName = []; end
if ~exist('FileName', 'var'); FileName = []; end

% set up parpool
ppobj = setup_parpool(roipars.serId, roipars.corenum);

% Selecting folders
f2run = dir;
f2run = str2match(FolderName, f2run);
f2run = str2rm(roipars.fo2reject, f2run);
f2run = {f2run.name};

timeint = tic;
fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']); 
    cd(f2run{i});
    runperfolder(FileName, roiparams, roipars);
    cd(roipars.cDir)
    fprintf('\n')
    
end

if ~ispc || ~ismac
    delete_parpool(ppobj);
end

fprintf('Done\n');
toc(timeint)

end

function runperfolder(fname, roiparams, roipars)
% runperfolder: run all files per folder
%
% Usage:
%   runperfolder(fname, roiparams, roipars)
%
% Args:
%   fname: file name pattern
%   roiparams: roi segmentation paramters (see roiseg_params)
%   roipars: roi file selection and job related params

% determine if fname narrows down to just one file
[f2plot, ~, rep2plot] = ...
    rdir_namesplit(fname, '.mat', ...
    roipars.fisuffix, roipars.fi2reject, [], 1);

fprintf('ROI segmentation of ');

if numel(f2plot) == 1 && numel(rep2plot) == 1
    
    % run single file
    f2plot{1} = [f2plot{1}, '_', num2str(rep2plot)];
    fprintf(' one file_seg ')
    
else
    
    % Run all files per folder
    f2plot = unique(f2plot);
    
end

fprintf([num2str(numel(f2plot)), ' files\n'])

for fly_i = 1:numel(f2plot)
    
    [filename, ~, repnum] = ...
        rdir_namesplit(f2plot{fly_i}, ...
        '.mat', roipars.fisuffix, ...
        roipars.fi2reject, [], 1);
    filename = unique(filename);
    repnum = sort(repnum);
    
    if numel(filename) == 1
        
        for rep_i = 1:numel(repnum)
            
            fprintf(['ROI segmenting file: ', ...
                filename{fly_i}, '_', ...
                num2str(repnum(rep_i)),'\n'])
            
            roiseg_int([filename{fly_i}, '_', ...
                num2str(repnum(rep_i))], ...
                roiparams, roipars);
            
        end
        
    else
        
        fprintf('error')
        
    end
    
end

end

function roiseg_int(filename, roiparams, roipars)
% roiseg_int: segment ROIs from 2DxT or
%   3DxT data (not big data though)
%
% Usage:
%   roiseg_int(fname, roiparams, roipars)
%
% Args:
%   filename: file name pattern
%   roiparams: roi segmentation paramters (see roiseg_params)
%   roipars: roi file selection and job related params
%
% Notes:
% runs ROI segmentation per mat file
%   it requires the mat file to be in the right mat format


% make object and update parameters
loDir = pwd;

if isfield(roiparams, 'tDir') && ...
        ~isempty(roiparams.tDir)
    loDir = roiparams.tDir;
    fprintf(['Target directory : ', loDir, '.\n']);
end

obj = CaROISeg;

% load metadata
load([filename, roipars.fimetsuffix], 'wDat')
dsize = wDat.vSize;

% load mat file (load and generate temp mat file)
loaddata(obj, [filename, roipars.fosuffix]) 

% make patches (generate patches using set parameters)
makepatch(obj, roiparams.patch_size, roiparams.overlap)

% update parameters
updateParams(obj, ...
    'search_method', roiparams.search_method, ...
    'dist', roiparams.dist, ...
    'se', roiparams.se, ...
    'deconv_method', 'constrained_foopsi', ...
    'temporal_iter', roiparams.temporal_iter, ...
    'ssub', roiparams.ssub, ...
    'tsub', roiparams.tsub, ...
    'fudge_factor', roiparams.fudge_factor, ...
    'merge_thr', roiparams.merge_thr, ...
    'gSig', roiparams.tau, ...
    'gSiz', roiparams.gSiz, ...
    'init_method', roiparams.init_method, ...
    'cluster_pixels', roiparams.cluster_pixels, ...
    'noise_norm', roiparams.noise_norm, ...
    'noise_range', roiparams.noise_range, ...
    'nrgthr', roiparams.nrgthr, ...
    'thr_method', roiparams.thr_method, ...
    'rem_prct', roiparams.rem_prct, ...
    'df_prctile', roiparams.df_prctile, ...
    'spatial_method', roiparams.spatial_method);

% add extra options 
obj.options.sr = roiparams.sr;
obj.options.freq_th = roiparams.freq_th;
obj.options.brainmask = wDat.bMask;

% add extra options for detending
obj.options.dtype = roiparams.dtype; % 0 = lowpass filter, 1 = percentile filter
obj.options.d_prct = roiparams.d_prct;
obj.options.d_window = roiparams.d_window;
obj.options.d_shift = roiparams.d_shift;
obj.options.se_parallel = 1;
obj.options.Asize_ths = roiparams.Asize_ths;
obj.options.t_init = 2;

% add c_ths for data with LED bleedthrough
if isfield(wDat.sPars, 'led_delta') && ...
        ~isempty(wDat.sPars.led_delta)
    
    if isfield(roipars, 'c_ths') && ...
            ~isempty(roipars.c_ths)
        obj.options.l_df = roipars.c_ths;
    else
        obj.options.l_df = max(wDat.sPars.led_delta);
    end
    
else
    
    obj.options.l_df = 0;
    
end

if isnan(obj.options.l_df)
    obj.options.l_df = 0;
end

% 1) run each file
% 1.1) do roi segmentation for each patch (XZY sub-segment)
% this generates the variable: RESULTS and saves files *_t_*.mat
fprintf(['Running patches: ', num2str(numel(obj.patches)),'\n'])

for patch_idx = 1:length(obj.patches)
    
    try
        
        if ~exist([filename, '_t_', num2str(patch_idx), '.mat'], 'file')
            
            initComponents_part1(obj, roiparams.K, ...
                roiparams.tau, roiparams.p, patch_idx, filename)
            % edit run_CNMF_patches_int_runperseg
            
        else
            
            fprintf(['Already processed patch # ', ...
                num2str(patch_idx), ' out of ',...
            num2str(length(obj.patches)), '.\n']);
            disp([filename, '_t_', num2str(patch_idx), '.mat'])
            
        end
        
    catch error
        
        fprintf(['Failed processing patch # ', ...
            num2str(patch_idx), ' out of ', ...
            num2str(length(obj.patches)), '.\n']);
        error.identifier
        error.message
        error.cause
        
    end
    
end

% 2.2) load *_t_*.mat files, compile patches, save *_prosroi
if ~exist([loDir, filesep, filename, '_prosroi.mat'], 'file') || roipars.redo
    
    initComponents_part2(obj, roiparams.K, ...
        roiparams.tau, roiparams.p, filename)
    
    rawsignal(obj, [filename, roipars.rfosuffix])
    
    bas_estimate(obj);
    roi = snapshot(obj);
    roi.userparams = roiparams;
    
    % saves files locally
    save([loDir, filesep, filename, '_prosroi'], 'roi', '-v7.3')
    
    % plot ROI coverage results
    if ~exist([roipars.cDir, filesep, 'roicov'], 'dir')
        mkdir([roipars.cDir, filesep, 'roicov'])
    end
    
    plot_roi_coverage_int(filename, ...
        [roipars.cDir, filesep, 'roicov'], ...
        roipars.fimetsuffix);
    
else
    
    fprintf(['Already compiled fly: ', filename, '.\n']);
    disp([loDir, filesep, filename, '_prosroi.mat'])
    
end

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
