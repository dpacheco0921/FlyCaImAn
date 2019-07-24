function batch_NoRMCorre(FolderName, FileName, iparams)
% batch_NoRMCorre: function to perform motion correction to files within a
% folder
%
% Usage:
%   batch_NoRMCorre(FolderName, FileName, iparams)
%
% Args:
%   FolderName: name of folders to load
%   FileName: name of files to load
%   iparams: parameters to update
%       (cDir: current directory)
%       (redo: redo gate)
%       (debug: debug gate)
%       (fsuffix: suffix of files to load)
%       (fo2reject: folders to reject)
%       %%%%%%%%%%%% parpool & server related %%%%%%%%%%%%
%       (serId: server id)
%           (default, 'int')
%       (corenum: number of cores)
%           (default, 4)
%       %%%%%%%%%%%% motion corretion params %%%%%%%%%%%%
%       (stack2del: timepoints to delete (for first frames when the laser warms up))
%           (default, [])
%       (rigidg: ridig motion correction gate)
%           (default, 1)
%       (nrigidg: nonridig motion correction gate)
%           (default, 0)
%       (refcha: reference channel)
%           (red channel, 1, default)
%           (green channel, 2)
%       %%%%%%%%%%%% internal NoRMCorre params %%%%%%%%%%%%
%       (def_maxshift: max shift allowed from frame-to-frame)
%           (default, [5 5 2])
%       (phaseflag: flag for using phase correlation, use this when SNR is high)
%           (default, 1)
%       (shifts_method: method to use for resampling)
%           ('linear', 'fft', 'cubic')
%       (boundary: method of boundary treatment)
%           (default, 'NaN')
%       (grid_size: size of non-overlapping regions)
%           (default, [64, 64, 5])
%       (bin_width: width of each bin
%           (deafult, 50)
%       (mot_uf: degree of patches upsampling)
%           (deafult, [4 4 1])
%       (us_fac: upsampling factor for subpixel registration)
%           (deafult, 10)
%       (overlap_pre: size of overlapping region)
%           (deafult, [16 16 3])
%       (overlap_post: size of overlapping region after upsampling)
%           (deafult, [16 16 3])
%       (use_parallel: for each frame, update patches in parallel)
%           (deafult, true)
%       (min_patch_size: minimum size of patch)
%           (deafult, [64 64 5])
%       %%%%%%%%%%%% extra editing %%%%%%%%%%%%
%       (shift_ths: max cumulative motion or delta from fram-to-frame (assummes smoothness))
%           (default, [18, 0.7])
%       (span: number of frames to smooth)
%           (penalizes fast motion, assumes it is mostly slow motion == drift))
%           (default, 10)
%       (sgate: type of editing)
%           (1 = smooth and zero big displacements)
%           (2 = smooth, 3 = zeroing)
%           (0 = raw)
%       (readextfile_flag: flag to read shifts from another matfile and apply those to this file)
%       (readextfile_dir: input directory)
%
% Notes:
% This function uses NoRMCorre (https://github.com/flatironinstitute/NoRMCorre)
% see NoRMCorre functions: NoRMCorreSetParms, normcorre_batch
%
% ToDo:
% add optical flow
% re-run old correction readextfile_flag/readextfile_dir

% default params
pMC = [];
pMC.cDir = pwd;
pMC.redo = 0;
pMC.debug = 0;
pMC.fsuffix = '_rawdata.mat';
pMC.fo2reject = {'.', '..', 'preprocessed', 'BData'};
pMC.serId = 'int';
pMC.corenum = 4; 
pMC.stack2del = [];
pMC.rigidg = 1;
pMC.nrigidg = 0;
pMC.refcha = 1;
pMC.def_maxshift = [5 5 2];
pMC.phaseflag = 1;
pMC.shifts_method = 'linear';
pMC.boundary = 'NaN';
pMC.grid_size = [64, 64, 5];
pMC.bin_width = 50;
pMC.mot_uf = [4 4 1];
pMC.us_fac = 10;
pMC.overlap_pre = [16 16 3];
pMC.overlap_post = [16 16 3];
pMC.use_parallel = true;
pMC.min_patch_size = [64 64 4];
pMC.shift_ths = [18, 0.7];
pMC.span = 10;
pMC.sgate = 1;
pMC.readextfile_flag = 0;
pMC.readextfile_dir = [];

% update variables
if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('FileName', 'var'); FileName = []; end
if ~exist('iparams', 'var'); iparams = []; end
pMC = loparam_updater(pMC, iparams);

% start pararell pool if not ready yet
ppobj = setup_parpool(pMC.serId, pMC.corenum);

% Selecting folders
f2run = dir;
f2run = str2match(FolderName, f2run);
f2run = str2rm(pMC.fo2reject, f2run);
f2run = {f2run.name};

fprintf('Motion Correction NoRMCorre\n')
fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']);
    cd(f2run{i});
    runperfolder(FileName, pMC);
    cd(pMC.cDir)
    
end

delete_parpool(ppobj);

fprintf('... Done\n')

end

function runperfolder(fname, pMC)
% runperfolder: run all files per folder
%
% Usage:
%   runperfolder(fname)
%
% Args:
%   fname: file name pattern
%   pMC: parameter variable

global RedCha GreenCha

shifts_r = [];
shifts_nr = [];
RedCha = [];
GreenCha = [];

% Files to load
sep = filesep;
f2run = rdir(['.', sep, '*', pMC.fsuffix]);
fname = addsuffix(fname, pMC.fsuffix);
f2run = str2match(fname, f2run);
f2run = {f2run.name};
fprintf('Correcting motion\n')

for F2R_idx = 1:numel(f2run)
    
    % Loading files
    fprintf(['Correcting file: ', ...
        f2run{F2R_idx}((max(strfind(f2run{F2R_idx}, sep))+1):end), '\n'])
    % Loading metadata
    load(strrep(f2run{F2R_idx}, '_rawdata', '_metadata'), ...
        'iDat', 'lStim', 'fDat')
    
    if iDat.MotCorr == 0 || pMC.redo == 1
        
        % Loading Data (uint16) and keeping the format as double
        %   (cos some functions require double precision)
        if pMC.redo == 1 && exist(strrep(f2run{F2R_idx}, '_rawdata', '_refdata'), 'file')
            
            load(strrep(f2run{F2R_idx}, '_rawdata', '_metadata'), 'mcDat')
            green_obj = matfile(f2run{F2R_idx}, 'Writable', true);            
            GreenCha = green_obj.Data; 
            red_obj = matfile(strrep(f2run{F2R_idx}, ...
                '_rawdata', '_refdata'), 'Writable', true);            
            RedCha = red_obj.Data;
            
        else
            
            data_obj = matfile(f2run{F2R_idx}, 'Writable', true);
            % Delete frame during flyback / update frameN / RedChaMean
            iDat = deleteflybackframe(iDat, fDat, data_obj);
            
            if ~exist('lStim', 'var')
                lStim = [];
            end
            
            [lStim, iDat] =  volumeprunner(lStim, iDat, data_obj, pMC.stack2del);
            % update lStim and iDat
            save(strrep(f2run{F2R_idx}, '_rawdata', '_metadata'), ...
                'iDat', 'lStim', '-append')
        
            siz = size(data_obj.Data);
            ChannelSplitter(siz, data_obj)
            
        end
        
        % fill in gaps in RedCha if it is opto data
        %   (when opto stim is overlapping with red PMT)
        fillgaps_redcha(fDat, iDat, pMC.debug)
        
        % Do motion correction
        dgDim = size(GreenCha);

        % motion correction parameters

        if length(dgDim) == 4
            d3 = dgDim(3);
        else
            d3 = 1;
            pMC.overlap_pre = pMC.overlap_pre(1);
            pMC.overlap_post = pMC.overlap_post(1);
            pMC.grid_size(3) = 1;
        end
        
        % rigid correction settings
        options_r = NoRMCorreSetParms(...
            'd1', dgDim(1), 'd2', dgDim(2), 'd3', d3, ...
            'grid_size', [dgDim(1:2) d3], 'bin_width', pMC.bin_width, ...
            'mot_uf', pMC.mot_uf, 'us_fac', pMC.us_fac, ...
            'overlap_pre', pMC.overlap_pre, 'overlap_post', pMC.overlap_post, ...
            'use_parallel', pMC.use_parallel, 'max_shift', pMC.def_maxshift, ...
            'phase_flag', pMC.phaseflag, 'boundary', pMC.boundary, ...
            'shifts_method', pMC.shifts_method);
        options_r.correct_bidir = 0;

        % non-rigid correction settings
        options_nr = NoRMCorreSetParms(...
            'd1', dgDim(1), 'd2', dgDim(2), 'd3', d3, ...
            'grid_size', pMC.grid_size, 'bin_width', pMC.bin_width, ...
            'mot_uf', pMC.mot_uf, 'us_fac', pMC.us_fac, ...
            'overlap_pre', pMC.overlap_pre, 'overlap_post', pMC.overlap_post, ...
            'min_patch_size', pMC.def_maxshift, 'use_parallel', pMC.use_parallel, ...
            'max_shift', pMC.def_maxshift, 'phase_flag', pMC.phaseflag, ...
            'boundary', pMC.boundary, 'shifts_method', pMC.shifts_method);
        options_nr.correct_bidir = 0;
        
        % crispness of the mean        
        if pMC.refcha == 1
            mcDat.crisp(1, 1) = get_crisp_idx(RedCha);
        else
            mcDat.crisp(1, 1) = get_crisp_idx(GreenCha);
        end

        % Run rigid
        if pMC.rigidg
            
            fprintf('RigidMC\n')
            % correct for motion (using selected ref channel)
            % shifts: [height, width, depth]
            if pMC.refcha == 1
                [~, shifts_pre, template_, ~, col_shift] = ...
                    normcorre_batch(RedCha, options_r);          
            else
                [~, shifts_pre, template_, ~, col_shift] = ...
                    normcorre_batch(GreenCha, options_r);           
            end
            
            % correlation with the mean (CM):
            if pMC.refcha == 1
                mcDat.CM(:, 1) = get_CM(template_, RedCha);
            else
                mcDat.CM(:, 1) = get_CM(template_, GreenCha);
            end
            
            % shifts_g = struct('shifts',cell(T,1),'shifts_up',cell(T,1),'diff',cell(T,1));
            
            % Edit shifts
            igate = 0;
            if pMC.sgate == 1 % smooth shift and zero it if the delta is too big
                fprintf('\n smoothing and zeroing if necessary\n')
                [shifts_r, shifts_s, igate] = ...
                    checkshift(shifts_pre, pMC.span, pMC.shift_ths);
            elseif pMC.sgate == 2 % just smooth shift
                fprintf('\n smoothing\n');
                shifts_r = smoothshift(shifts_pre, [], pMC.span);
            elseif pMC.sgate == 3 % just zero shift
                fprintf('\n zeroing\n')
                shifts_r = zeroshift(shifts_pre);
            else % use raw
                fprintf('\n not editing\n');
                shifts_r = shifts_pre;
            end
            
            if ~exist('shifts_s', 'var'); shifts_s = shifts_r; end
            
            if pMC.debug
                
                % plot shifts
                figName = strrep(strrep(strrep(f2run{F2R_idx}, ...
                    ['.', sep], ''), '_rawdata.mat', ''), '_', '-');
                
                if length(dgDim) == 4 
                    figH = plot_NoRMCorr_shitfs(3, ...
                        shifts_pre, shifts_s, shifts_r);
                    %[figH, aH] = plot_NoRMCorr_shitfs(im_dim, varargin)
                else
                    figH = plot_NoRMCorr_shitfs(2, ...
                        shifts_pre, shifts_s, shifts_r);
                end
                
                figH.Name = figName;
                
            end
            
            if pMC.debug
                % Plot video of single plane                
                %keyboard
                %plot_2DoR(1, RedCha, GreenCha, [], 0.001, {'RedCha','GreenCha'})
                %plot_2DoR(5, GreenCha, Ymr, [], 0.001, {'Raw','Corrected'})
                
            end
            
            % apply shifts to Y           
            if pMC.refcha == 1
                RedCha = apply_shifts(RedCha, shifts_r, options_r, [], [], [], col_shift);
            elseif pMC.refcha == 2
                GreenCha = apply_shifts(GreenCha, shifts_r, options_r, [], [], [], col_shift);
            end
            
            if pMC.refcha == 1 && ~isempty(GreenCha)
                GreenCha = apply_shifts(GreenCha, shifts_r, options_r, [], [], [], col_shift);
            elseif pMC.refcha == 2 && ~isempty(RedCha)
                RedCha = apply_shifts(RedCha, shifts_r, options_r, [], [], [], col_shift);
            end
            
        end
        
        % Run non rigid  still in progress
%         if pMC.nrigidg
%             
%             fprintf('NonRigidMC\n')
%             if pMC.refcha == 1
%                 if length(dgDim) == 4
%                     Y = RedCha(:, :, :, perm);
%                 else
%                     Y = RedCha(:, :, perm);
%                 end
%             else
%                 if length(dgDim) == 4
%                     Y = GreenCha(:, :, :, perm);
%                 else
%                     Y = GreenCha(:, :, perm);
%                 end
%             end
%             
%             [Y, shifts_nr, ~] = normcorre_batch(Y, options_nr);
%             
%             if pMC.refcha == 1 
%                 if length(dgDim) == 4
%                     RedCha(:, :, :, perm) = Y;
%                 else
%                     RedCha(:, :, perm) = Y;
%                 end
%             else
%                 if length(dgDim) == 4
%                     GreenCha(:, :, :, perm) = Y;
%                 else
%                     GreenCha(:, :, perm) = Y;
%                 end
%             end
%             
%             shifts_nr(perm)  = shifts_nr;
%             
%             if pMC.refcha == 1 && ~isempty(GreenCha)
%                 GreenCha = apply_shifts(GreenCha, shifts_nr, options_nr);
%             elseif pMC.refcha == 2 && ~isempty(RedCha)
%                 RedCha = apply_shifts(RedCha, shifts_nr, options_nr);
%             end
%             
%             clear Ymr
%             
%         end
        
        % Save metadata, avg image and save
        fprintf('Saving ... ')
        
        % Get mean volumes
        iDat.GreenChaMean = mean(GreenCha, length(dgDim));
        iDat.RedChaMean = mean(RedCha, length(dgDim));
        
        % Get displacements
        iDat.MotCorr = 1;
        shifts = parseshift(shifts_r, shifts_pre, shifts_nr);
        mcDat.axes = {'Y', 'X', 'Z'};
        
        if pMC.redo == 1 && ...
                exist(strrep(f2run{F2R_idx}, '_rawdata', '_refdata'), 'file')
            mcDat.rigid = mcDat.rigid + shifts{1};
            mcDat.nonrigid = mcDat.nonrigid + shifts{3}; 
            mcDat.rigids = mcDat.rigids + shifts{2};
        else
            mcDat.rigid = shifts{1};
            mcDat.nonrigid = shifts{3};
            mcDat.rigids = shifts{2};
        end
        mcDat.ergate = igate;
        
        % remove nan pixels and whole slices
        nan_mask = max(isnan(iDat.GreenChaMean), [], 3);
        floatIm = [];
        
        if length(dgDim) == 4
            nan_plane = sum(reshape(isnan(iDat.GreenChaMean), ...
                [prod(dgDim([1 2])) dgDim(3)])) ~= prod(dgDim([1 2]));
            
            % use non-nan planes
            nan_mask = max(isnan(iDat.GreenChaMean(:, :, nan_plane)), [], 3);
            template_ = template_(:, :, nan_plane);
            
            if pMC.refcha == 1
                floatIm = pruneIm(RedCha(:, :, nan_plane, :), nan_mask);
            else
                floatIm = pruneIm(GreenCha(:, :, nan_plane, :), nan_mask);
            end
            
        else
            
            if pMC.refcha == 1
                floatIm = pruneIm(RedCha, nan_mask);
            else
                floatIm = pruneIm(GreenCha, nan_mask);
            end            
            
        end
        template_ = pruneIm(template_, nan_mask);
        
        % correlation after correction      
        mcDat.CM(:, 2) = get_CM(template_, floatIm);
                
        % crispness of the mean
        mcDat.crisp(1, 2) = get_crisp_idx(floatIm);
        
        clear floatIm
        
        % add optic flow
        % edit opticalFlowFarneback
        
        % saving processed video to .mat file
        save(strrep(f2run{F2R_idx}, '_rawdata', '_metadata'), ...
            'iDat', 'lStim', 'mcDat', '-append')
        
        Data = GreenCha;
        save(f2run{F2R_idx}, 'Data', '-v7.3');
        GreenCha = [];
        
        Data = RedCha;
        if ~isempty(Data)
            save(strrep(f2run{F2R_idx}, '_rawdata', '_refdata'), ...
                'Data', '-v7.3');
        end
        RedCha = [];
        
        Data = [];
        shifts_r = [];
        shifts_nr = [];
        
        clear shifts_s;
        fprintf('Done\n')
        
    else
        
        fprintf(' *already corrected*\n')
        
    end
    
    clear dDim iDat lStim fDat
    
end

fprintf('****** Done ******\n')

end

function CM = get_CM(templateIm, floatIm)
% get_CM: calculate correlation of template to each frame of floating image
%
% Usage:
%   CM = get_crisp_idx(iDat, fDat, data_obj)
%
% Args:
%   templateIm: reference channel to use
%   floatIm: floating image
%
% Returns:
%   CM: correlation

siz_ = size(floatIm);
templateIm = reshape(templateIm, [prod(siz_(1:end-1)) 1]);
floatIm = reshape(floatIm, [prod(siz_(1:end-1)) siz_(end)]);

CM = corr(templateIm, floatIm)';

end

function crisp_idx = get_crisp_idx(Y)
% get_crisp_idx: calculate crispness index as in (http://dx.doi.org/10.1016/j.jneumeth.2017.07.031)
%
% Usage:
%   crisp_idx = get_crisp_idx(Y)
%
% Args:
%   Y: input 3D or 2D matrix
%
% Returns:
%   crisp_idx: crispness

crisp_idx = [];

siz_ = size(Y);
pre_matrix = abs(reshape(squeeze(mean(Y, length(siz_))),...
    [prod(siz_(1:end-1)) 1]));

crisp_idx = norm(pre_matrix, 'fro');

end

function iDat = deleteflybackframe(iDat, fDat, data_obj)
% deleteflybackframe: delete deleting z frame and related time stamps (Piezzo fly back)
%
% Usage:
%   iDat = deleteflybackframe(iDat, fDat, data_obj)
%
% Args:
%   iDat: image metadata variable
%   fDat: metadata variable
%   data_obj: memmap data object
%
% Returns:
%   iDat: updated image metadata variable
%
% Notes:
% choose plane to delete based on data input type 
%   (Also, if already LEDcorr then this step has already happened)

if contains(fDat.DataType, '3DxT') && ~iDat.LEDCorr ...
        && ~contains(fDat.DataType, 'remotefocus')
    
    if contains(fDat.DataType, 'old')
        z2del = size(data_obj.Data, 3);
        % old, Chop Data and update iDat
        iDat = rmfield(iDat, 'sstEn');
    else
        z2del = 1;
        % new, Chop Data and update iDat
        iDat = rmfield(iDat, 'sstEn');
    end
    
    % update frame time / Data / RedChaMean
    fprintf('Deleting flyback frame, final size: ')
    iDat.fstEn(z2del:iDat.FrameN:end, :) = [];
    
    % load data
    data = data_obj.Data;
    
    % delete flyback frame
    data(:, :, z2del, :, :) = [];
    
    % overwrite data
    data_obj.Data = data;
    
else
    
    if iDat.LEDCorr
        fprintf('flyback frame deletion already performed, final size: ')
    elseif contains(fDat.DataType, 'remotefocus')
        fprintf('there is no flyback frame, so it does not require slice removal ')
    else
        fprintf('It is 2DxT data, so it does not require slice removal ')
    end
    
end

if contains(fDat.DataType, '3DxT')
    iDat.FrameN = size(data_obj.Data, 3);
elseif contains(fDat.DataType, '2DxT')
    iDat.FrameN = 1;
end

fprintf([num2str(iDat.FrameN), '\n'])

end

function ChannelSplitter(siz, data_obj)
% ChannelSplitter: splits data_obj.Data variable into channels that
%   comprise it.
%
% Usage:
%   ChannelSplitter(siz, data_obj)
%
% Args:
%   siz: image size
%   data_obj: memmap data object

global RedCha GreenCha

if length(siz) == 4 && siz(end) == 2
    
    RedCha = single(data_obj.Data(:, :, :, 1));
    GreenCha = single(data_obj.Data(:, :, :, 2));
    
elseif length(siz) == 5 && siz(end) == 2
    
    RedCha = single(data_obj.Data(:, :, :, :, 1));
    GreenCha = single(data_obj.Data(:, :, :, :, 2));
    
else
    GreenCha = single(data_obj.Data);
end

end

function [lStim, iDat] = ...
    volumeprunner(lStim, iDat, data_obj, stack2del)
% volumeprunner: prunes the original Data
%   and adjust lStim and iDat metadata accordingly
%
% Usage:
%   [lStim, iDat] =  ...
%       volumeprunner(lStim, iDat, data_obj, stack2del)
%
% Args:
%   lStim: stimuli metadata variable
%   iDat: image metadata variable
%   data_obj: memmap data object
%   stack2del: timepoints to delete
%
% Return:
%   lStim: updated stimuli metadata variable
%   iDat: updated image metadata variable

% Prune Data
if ~isempty(stack2del)
    
    % load data
    data = data_obj.Data;
    
    if iDat.StackN > 1
        data(:, :, :, stack2del, :) = [];
    else
        data(:, :, :, stack2del) = [];
    end
    
    % overwrite data
    data_obj.Data = data;
    
    % Rechape fsTen
    initF = reshape(iDat.fstEn(:, 1), [iDat.FrameN, iDat.StackN]);
    initF(:, stack2del) = [];
    initF = initF(:);
    
    endF = reshape(iDat.fstEn(:, 2), [iDat.FrameN, iDat.StackN]);
    endF(:, stack2del) = [];
    endF = endF(:);
    
    iDat.fstEn = [];
    iDat.fstEn(:, 1) = initF;
    iDat.fstEn(:, 2) = endF;
    
    % Prune lStim
    if ~isempty(lStim) && isfield(lStim, 'trace') ...
            && isfield(lStim, 'lstEn')
        lStim.trace = lStim.trace(iDat.fstEn(1, 1):iDat.fstEn(end, 2));
        lStim.lstEn = lStim.lstEn - iDat.fstEn(1, 1) + 1;
    end
    
    % Zero frame init
    iDat.fstEn = iDat.fstEn - iDat.fstEn(1, 1) + 1;
    
    % Get frame times and prune them too
    iDat.StackN = iDat.StackN - numel(stack2del);
    
    % Update sstEn
    if iDat.StackN > 1
        preInit = min(reshape(iDat.fstEn(:, 1), ...
            [iDat.FrameN, iDat.StackN]), [], 1)';
        preEnd =  max(reshape(iDat.fstEn(:, 1), ...
            [iDat.FrameN, iDat.StackN]), [], 1)';
        iDat.sstEn = [preInit, preEnd];
        clear preInit preEnd
    end
    
    % Update PMT_fscore
    if isfield(iDat, 'PMT_fscore') && ~isempty(iDat.PMT_fscore)
        iDat.PMT_fscore(:, stack2del) = [];
    end
    
end

end

function fillgaps_redcha(fDat, iDat, idebug)
% fillgaps_redcha: fill gaps in red channel for opto data from LEDcontroler
%
% Usage:
%   fillgaps_redcha(fDat, iDat, idebug)
%
% Args:
%   fDat: stimuli metadata variable
%   iDat: image metadata variable
%   idebug: debug gate

global RedCha

if contains(fDat.DataType, 'opto') && ~isempty(RedCha)
    
    tp2fill = [];
    
    if contains(fDat.DataType, '3DxT') && isfield(iDat, 'PMT_fscore')
        tp2fill = find(iDat.PMT_fscore < (iDat.FrameN + 1));
    elseif contains(fDat.DataType, '2DxT') && isfield(iDat, 'PMT_fscore')
        tp2fill = find(iDat.PMT_fscore < 1);
    end
    
    if idebug
        Data4display = squeeze(max(max(RedCha, [], 1), [], 2));
        subplot(1, 2, 1);
        imagesc(Data4display);
        title('Red_input')
    end
    
    RedCha = framegapfill(tp2fill, RedCha);
    
    if idebug
        Data4display = squeeze(max(max(RedCha, [], 1), [], 2));
        subplot(1, 2, 2);
        imagesc(Data4display);
        title('Red_output')
    end
    
end

end

function [shifts_i, shifts_s, igate] = ...
    checkshift(shifts_i, smooth_span, shift_ths)
% checkshift: get all shifts and smooth the change over time
%
% Usage:
%   [shifts_i, shifts_s, igate] = checkshift(shifts_i)
%
% Args:
%   shifts_i: input shifts
%   smooth_span: span for smoothing
%   shift_ths: threshold for the shifts
%
% Returns:
%   shifts_i, shifts_s, igate: zeroed shifts in vector format

igate = 0;
shifts_s = shifts_i;
shifts_r = shifts_editor(shifts_i, 'read');

for i = 1:size(shifts_r, 2)
    
    % Smooth
    shifts_sr(:, i) = smooth(shifts_r(:, i), smooth_span);
    % Round to 2 decimal
    shifts_sr(:, i) = round(shifts_sr(:, i)*100)/100;
    
end

shifts_delta = abs(shifts_sr - shifts_r);

for i = 1:size(shifts_r, 2)
    display([max(shifts_delta(:, i)) sum(shifts_delta(:, i))])
    
    if max(shifts_delta(:, i)) > shift_ths(2) ...
            || sum(shifts_delta(:, i)) > shift_ths(1)
        shifts_sr(:, i) = 0;
        igate = 1;     
    end
    
end

% if a max delta is greater than ths then zero all those shifts
% save smooth trace
shifts_s = shifts_editor(shifts_s, 'write_r', shifts_sr);
shifts_i = shifts_editor(shifts_i, 'write_r', shifts_sr);

end

function shifts_s = smoothshift(shifts_i, ax2use, smooth_span)
% smoothshift: smooth shifts, axes specific
%
% Usage:
%   shifts_s = smoothshift(shifts_i, ax2use, smooth_span)
%
% Args:
%   shifts_i: input shifts
%   ax2use: axes to use
%   smooth_span: span for smoothing
%
% Returns:
%   shifts: smoothed shifts in vector format

shifts_s = shifts_i;
shifts_r = shifts_editor(shifts_i, 'read');
shifts_sr = shifts_r;

if ~exist('ax2use', 'var') || isempty(ax2use)
    ax2use = 1:size(shifts_r, 2);
end

for i = ax2use
    % Smooth
    shifts_sr(:, i) = smooth(shifts_r(:, i), smooth_span);
    % Round to 2 decimal
    shifts_sr(:, i) = round(shifts_sr(:, i)*100)/100;
end

% save smooth trace
shifts_s = shifts_editor(shifts_s, 'write_r', shifts_sr);

end

function shifts_s = zeroshift(shifts_i, ax2use)
% zeroshift: zero shifts, axes specific
%
% Usage:
%   shifts_s = smoothshift(shifts_i, ax2use)
%
% Args:
%   shifts_i: input shifts
%   ax2use: axes to use
%
% Returns:
%   shifts: zeroed shifts in vector format

shifts_s = shifts_i;
shifts_r = shifts_editor(shifts_i, 'read');

if ~exist('ax2use', 'var')
    ax2use = 1:size(shifts_r, 2);
end

shifts_sr = shifts_r;

% Zeroing values
for i = ax2use
    shifts_sr(:, i) = shifts_r(:, i)*0;
end

% save trace
shifts_s = shifts_editor(shifts_s, 'write_r', shifts_sr);

end

function shifts = parseshift(varargin)
% parseshift: get shifts in vector format per dimension
%
% Usage:
%   shifts = parseshift(varargin)
%
% Args:
%   varargin: shifts (in cells)
%
% Returns:
%   shifts: shifts in vector format

shifts = cell(numel(varargin), 1);

for i = 1:numel(varargin)
    
    % get structure shifts_i and plots them
    shifts_i = varargin{i};
    
    if ~isempty(shifts_i)
        shifts_t = shifts_editor(shifts_i, 'read');
        if ~iscell(shifts_t); shifts_t = shifts_t'; end
    else
        shifts_t = [];
    end
    
    shifts{i} = shifts_t;
    clear shifts_i shifts_t shifts_r
    
end

end
