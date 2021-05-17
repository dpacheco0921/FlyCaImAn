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
        roiparams.patch_size, roiparams.overlap, ...
        [], [], roiparams.patch_minsiz)
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

% add user defined ROI centers
if isfield(wDat, 'ROI_center_matrix') && ~isempty(wDat.ROI_center_matrix)
    
    obj.options.ROI_center_matrix = wDat.ROI_center_matrix;
    
    % overwrite K
    roiparams.K = numel(find(obj.options.ROI_center_matrix ~= 0));
    
end

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

% update target output directory
if isfield(p, 'oDir') && ~isempty(p.oDir)
    target_dir = p.oDir;
else
    target_dir = p.cDir;
end

fprintf(['Files are saved at: ', target_dir])

if p.jobpart == 1
    
    % move to target folder
    cd(target_dir)
    
    if ~exist([target_dir, filesep, filename, '_prosroi.mat'], 'file')

        % 1) run each patch and save intermediate results
    
        % do roi segmentation for each patch (XZY sub-segment)
        %   this generates the variable: RESULTS
        fprintf(['Running patches: ', num2str(patchidx),'\n'])

        for i = 1:length(patchidx)

            try

                if ~exist([target_dir, filesep, filename, '_t_', ...
                        num2str(patchidx(i)), '.mat'], 'file')

                    initComponents_part1(obj, roiparams.K, ...
                        roiparams.tau, roiparams.p, patchidx(i), filename)
                    % edit run_CNMF_patches_int_runperseg

                else

                    fprintf(['Already processed patch # ', ...
                        num2str(patchidx(i)), ' out of ',...
                        num2str(length(obj.patches)), '.\n']);
                    display([target_dir, filesep, filename, '_t_', ...
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
        display([target_dir, filesep, filename, '_prosroi.mat'])
        
    end
   
    % move back to origin folder
    cd(p.cDir)
    
elseif p.jobpart == 2
    
    % 2) Compile intermediate results
    
    if ~exist([target_dir, filesep, filename, '_prosroi.mat'], 'file')
        
        % move to target folder
        cd(target_dir)
    
        initComponents_part2(obj, roiparams.K, ...
            roiparams.tau, roiparams.p, filename)
        % edit run_CNMF_patches_int_compile
        
        % move back to origin folder
        cd(p.cDir)
        
        rawsignal(obj, [filename, p.rfisuffix])
        % bas_estimate(obj)
        roi = snapshot(obj);
        roi.userparams = roiparams;
        
        % saves files locally
        save([target_dir, filesep, filename, '_prosroi'], 'roi', '-v7.3')
        
        % plot ROI coverage results
        if ~exist([target_dir, filesep, 'roicov'], 'dir')
            mkdir([target_dir, filesep, 'roicov'])
        end
        
        plot_roi_coverage_int(...
            [p.cDir, filesep, filename, p.fimetsuffix], ...
            [target_dir, filesep, filename, '_prosroi'], ...
            [target_dir, filesep, 'roicov']);
        
    else
        
        fprintf(['Already compiled fly: ', filename, '.\n']);
        display([target_dir, filesep, filename, '_prosroi.mat'])
        
    end
    
elseif p.jobpart == 3

    % move to target folder
    cd(target_dir)
    
    % 3) run each patch and save results directly
    %   it does: 1), and saves all ROI segmentation 
    %   parameters as in 2)
    
    % uses internally:
    % roi = snapshot(obj);
    % roi.userparams = roiparams;

    fprintf(['Running patches: ', num2str(patchidx),'\n'])
    
    for i = 1:length(patchidx)
        
        try
            
            if ~exist([target_dir, filesep, filename, '_t_', ...
                    num2str(patchidx(i)), '.mat'], 'file')
                
                initComponents_part3(obj, roiparams.K, ...
                    roiparams.tau, roiparams.p, patchidx(i), ...
                    filename, p.rfisuffix)
                % edit run_CNMF_patches_selected_segment
                
            else
                
                fprintf(['Already processed patch # ', ...
                    num2str(patchidx(i)), ' out of ',...
                    num2str(length(obj.patches)), '.\n']);
                display([target_dir, filesep, filename, '_t_', ...
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
    
    % move back to origin folder
    cd(p.cDir)
    
elseif p.jobpart == 4
    
    % 4) Compile intermediate results without 
    %   merging across planes, it also removes 
    %   any weight from out of plane ROI pixels
    
    if ~exist([target_dir, filesep, filename, '_prosroi.mat'], 'file')
        
        % move to target folder
        cd(target_dir)
        
        initComponents_part4(obj, roiparams.K, ...
            roiparams.tau, roiparams.p, filename)
        % edit run_CNMF_patches_int_compile_2
        
        % move back to origin folder
        cd(p.cDir)
        
        rawsignal(obj, [filename, p.rfisuffix])
        % bas_estimate(obj)
        roi = snapshot(obj);
        roi.userparams = roiparams;
        
        % saves files locally
        save([target_dir, filesep, filename, '_prosroi'], 'roi', '-v7.3')
        
        % plot ROI coverage results
        if ~exist([target_dir, filesep, 'roicov'], 'dir')
            mkdir([target_dir, filesep, 'roicov'])
        end
        
        plot_roi_coverage_int(...
            [p.cDir, filesep, filename, p.fimetsuffix], ...
            [target_dir, filesep, filename, '_prosroi'], ...
            [target_dir, filesep, 'roicov']);
        
    else
        
        fprintf(['Already compiled fly: ', filename, '.\n']);
        display([target_dir, filesep, filename, '_prosroi.mat'])
        
    end
    
end

end

function plot_roi_coverage_int(...
    metadata_fullpath, roi_fullpath, oDir)
% plot_roi_coverage_int: plot roi coverage
%   it plot both sum of max norm weights and binarized voxels
%
% Usage:
%   plot_roi_coverage_int(...
%       metadata_fullpath, roi_fullpath, oDir)
%
% Args:
%   metadata_fullpath: full file path of *_metadata.mat
%   roi_fullpath: full file path of *_prosroi.mat
%   oDir: output directory

% load metadata file
load(metadata_fullpath, 'wDat')
% load ROI file
load(roi_fullpath, 'roi')

[filename, ~] = ...
    split_path(strrep(roi_fullpath, '_prosroi', ''));

% plot videos
plot_roi_coverage(filename, ...
    [1 1 1], wDat, roi, oDir)

end
