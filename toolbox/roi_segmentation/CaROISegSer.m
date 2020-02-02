function CaROISegSer(paramfile, serverid, IntID)
% CaROISegSer: runs internal functions for ROI 
%   segmentation of volumes to be partitioned
%   into subvolumes (defined by patches)
%
% Usage:
%   CaROISegSer(paramfile, serverid, IntID)
%
% Args:
%   paramfile: mat file with all parameters
%   serverid: server id
%   IntID: job index manually inserted
%
% Notes:
% This function is run by batch_CaROISegSer
%   cluster compatible

timeint = tic;
fprintf('Running roi-segmentation using CNMF\n');

p = [];
p.cDir = pwd;
p.tDir = pwd;

pwd

% load parameter file
if exist('IntID', 'var')
    
    taskID = IntID;
    load([p.tDir, filesep, paramfile, '_impre.mat'], ...
        'filename', 'patchidx', 'p', 'roiparams')

else
    
    [~, ~, ~, tempfiledir, ~] = ...
        user_defined_directories(serverid);
    p.tDir = [tempfiledir, filesep, 'jobsub', filesep, 'roirel'];
    load([p.tDir, filesep, paramfile, '_impre.mat'], ...
        'filename', 'patchidx', 'p', 'roiparams')
    % p.Cdir is updated to the one saved in the _impre.mat file
    taskID = getclus_taskid(serverid);
    
end

% Display parameters only when running from clusters
if ~ispc && ~ismac
    p
    roiparams
end

% start pararell pool if not ready yet
ppobj = setup_parpool(serverid, p.corenum);

% Selecting file and frames to use
filename = filename{taskID};
if p.jobpart == 1 || ...
        p.jobpart == 3 && ...
        ~isempty(patchidx)
    patchidx = patchidx{taskID};
else
    patchidx = [];
end

% pick task index
taskID = num2str(taskID);

% moving to data folder
fprintf('\nStep (1): loading folder info\n');
fprintf(['File to run: ', filename, '\n']);
cd(p.cDir)

% register lrefIm to refIm (affine)
fprintf('\nStep (2): do roi segmentation\n')
batch_CaROISeg_3D(filename, roiparams, patchidx, p)

if ~ispc || ~ismac
    delete_parpool(ppobj);
end

fprintf('Done\n');
toc(timeint)

end

function batch_CaROISeg_3D(...
    filename, roiparams, patchidx, p)
% CaROISegSer: runs ROI segmentation per mat file, 
%   it requires the mat file to be in the right mat format
%   parameters
%
% Usage:
%   batch_CaROISeg_3D(...
%       filename, roiparams, patchidx, p)
%
% Args:
%   filename: filename (usually *metadata.mat)
%   roiparams: roi segmentation paramters (see roiseg_params)
%   patchidx: index of patch to use (patch of whole original volume)
%   p: roi file selection and job related params
%
% Notes:
% data file usually has a 3DxT volume as Y and a flattened Y variable Yr

obj = CaROISeg;

% load metadata
load([filename, p.fimetsuffix], 'wDat')
dsize = wDat.vSize;

% load and generate temp mat file
loaddata(obj, [filename, p.fisuffix])

% generate patches using set parameters
if roiparams.patchtype == 0
    makepatch(obj, ...
        roiparams.patch_size, roiparams.overlap)
else
    makepatch(obj, ...
        roiparams.patch_size, roiparams.overlap, ...
        roiparams.patchtype, wDat.Zstitch.Zidx)
end

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
obj.options.d_prct = roiparams.d_prct;
obj.options.d_window = roiparams.d_window;
obj.options.d_shift = roiparams.d_shift;
obj.options.se_parallel = 1;
obj.options.Asize_ths = roiparams.Asize_ths;
% skip first 5 seconds
obj.options.t_init = max([sum(wDat.fTime - wDat.fTime(1)  < 5) 1]);

% add c_ths for data with LED bleedthrough
if isfield(wDat.sPars, 'led_delta') && ...
        ~isempty(wDat.sPars.led_delta)
    
    % roiparams.c_ths is userdefined threshold (deprecated)
    
    if isfield(roiparams, 'c_ths') && ...
            ~isempty(roiparams.c_ths)
        obj.options.l_df = roiparams.c_ths;
    else
        obj.options.l_df = max(wDat.sPars.led_delta);
    end
    
else
    
    obj.options.l_df = 0;
    
end

if isnan(obj.options.l_df)
    obj.options.l_df = 0;
end

if p.jobpart == 1
    
    if ~exist([p.cDir, filesep, filename, '_prosroi.mat'], 'file')

        % 1) run each patch and save intermediate results
    
        % do roi segmentation for each patch (XZY sub-segment)
        %   this generates the variable: RESULTS
        fprintf(['Running patches: ', num2str(patchidx),'\n'])

        for i = 1:length(patchidx)

            try

                if ~exist([p.cDir, filesep, filename, '_t_', ...
                        num2str(patchidx(i)), '.mat'], 'file')

                    initComponents_part1(obj, roiparams.K, ...
                        roiparams.tau, roiparams.p, patchidx(i), filename)
                    % edit run_CNMF_patches_int_runperseg

                else

                    fprintf(['Already processed patch # ', ...
                        num2str(patchidx(i)), ' out of ',...
                        num2str(length(obj.patches)), '.\n']);
                    display([p.cDir, filesep, filename, '_t_', ...
                        num2str(patchidx(i)), '.mat'])

                end

            catch error

                fprintf(['Failed processing patch # ', ...
                    num2str(patchidx(i)), ' out of ',...
                    num2str(length(obj.patches)), '.\n']);
                error.identifier
                error.message
                error.cause

            end

        end
    
    else
        
        fprintf(['Already compiled fly: ', filename, '.\n']);
        display([p.cDir, filesep, filename, '_prosroi.mat'])
        
    end
        
elseif p.jobpart == 2
    
    % 2) Compile intermediate results
    
    if ~exist([p.cDir, filesep, filename, '_prosroi.mat'], 'file')
        
        initComponents_part2(obj, roiparams.K, ...
            roiparams.tau, roiparams.p, filename)
        rawsignal(obj, [filename, p.rfisuffix])
        %bas_estimate(obj)
        roi = snapshot(obj);
        roi.userparams = roiparams;
        
        % saves files locally
        save([filename, '_prosroi'], 'roi', '-v7.3')
        
        % plot ROI coverage results
        if ~exist([p.cDir, filesep, 'roicov'], 'dir')
            mkdir([p.cDir, filesep, 'roicov'])
        end
        
        plot_roi_coverage_int(filename, ...
            [p.cDir, filesep, 'roicov'], ...
            p.fimetsuffix);
        
    else
        
        fprintf(['Already compiled fly: ', filename, '.\n']);
        display([p.cDir, filesep, filename, '_prosroi.mat'])
        
    end
    
elseif p.jobpart == 3

    % 3) run each patch and save results directly
    %   it does: 1), and saves all ROI segmentation 
    %   parameters as in 2)
    
    % uses internally:
    % roi = snapshot(obj);
    % roi.userparams = roiparams;

    fprintf(['Running patches: ', num2str(patchidx),'\n'])
    
    for i = 1:length(patchidx)
        
        try
            
            if ~exist([p.cDir, filesep, filename, '_t_', ...
                    num2str(patchidx(i)), '.mat'], 'file')
                
                initComponents_part3(obj, roiparams.K, ...
                    roiparams.tau, roiparams.p, patchidx(i), ...
                    filename, p.rfisuffix)
                % edit run_CNMF_patches_int_runperseg
                
            else
                
                fprintf(['Already processed patch # ', ...
                    num2str(patchidx(i)), ' out of ',...
                    num2str(length(obj.patches)), '.\n']);
                display([p.cDir, filesep, filename, '_t_', ...
                    num2str(patchidx(i)), '.mat'])
                
            end
            
        catch error
            
            fprintf(['Failed processing patch # ', ...
                num2str(patchidx(i)), ' out of ',...
                num2str(length(obj.patches)), '.\n']);
            error.identifier
            error.message
            error.cause
            
        end
        
    end
    
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
