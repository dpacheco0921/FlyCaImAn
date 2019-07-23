function batch_SpaTemp_ResFilt_3D(FolderName, FileName, iparams)
% batch_SpaTemp_ResFilt_3D: function to perform spatial and temporal
%   resampling for all files per folder. For segments of the same fly 
%   (based on flyname) it aligns them all to stimulus start. 
%
% Usage:
%   batch_SpaTemp_ResFilt_3D(FolderName, FileName, iparams)
%
% Args:
%   FolderName: Folder name to load
%   FileName: File name to load
%   iparams: parameters to update
%       (cDir: current directory)
%           (default, pwd)
%       (debug: debug gate)
%           (default, 0)
%       (fo2reject: folders to reject)
%           (default, {'.', '..', 'preprocessed', 'BData'})
%       (fi2reject: files to reject)
%           (default, {'Zstack'})
%       (fsuffix: suffix of files to load)
%           (default, '_rawdata')
%       %%%%%%%%%%%% temporal resampling %%%%%%%%%%%%
%       (newtimeres: final temporal resolution in sec)
%           (default, 0.5)
%       (time: internal variable that defines start and end of recordings per file)
%           (default, [])
%       %%%%%%%%%%%% spatial resampling & smoothing %%%%%%%%%%%%
%       (sigma: timepoints to delete (for first frames when the laser warms up))
%           (default, [])
%       (size: ridig motion correction gate)
%           (default, [])
%       (newres: target voxel size (homogenize width and heigth))
%           (default, [1.2 1.2 1])
%       (fres: max sub-micron resolution, up to 4 digits, it rounds number smaller than this)
%           (default, 10^4)
%       (art_val: negative value to use for nan replacements before spatial resampling)
%           (default, -1*10^3)
%       (wDat_gen: flag to generate wDar)
%           (default, 0)
%       (direction: used for iverting order of planes in the z axis, see notes)
%           (default, 'invert')
%       %%%%%%%%%%%% shift fluorescence distribution %%%%%%%%%%%%
%       (bkgate: flag for background substraction)
%           (default, 0)
%       (blowcap: fluorescence below which it is zerored)
%           (default, 0)
%       (fshift: shift distribution of F to the positive side)
%           (default, 6)
%       (stack2del: remove timepoints)
%           (default, [])
%
% Notes:
% Works specifically for 3DxT
% Assumes original data does not have negative values below "art_val".
% FileName means strictly file name (minus '.mat')
% How to orient consecutive planes and stacks (particularly to match reference orientation: ventral-dorsal)
%   Imaging from the central brain dorsal (need to be inverted to have a ventral-dorsal orientation)
%   Imaging from the VNC from ventral side already has a ventral-dorsal orientation.

% default params
spte = [];
spte.cDir = pwd;
spte.debug = 0;
spte.fo2reject = {'.', '..', 'preprocessed', 'BData'};
spte.fi2reject = {'Zstack'};
spte.fsuffix = '_rawdata';
spte.newtimeres = 0.5;
spte.time = [];
spte.sigma = [];
spte.size = [];
spte.newres = [1.2 1.2 1];
spte.fres = 10^4;
spte.art_val = -1*10^3;
spte.wDat_gen = 0;
spte.direction = 'invert';
spte.bkgate = 0;
spte.blowcap = 0;
spte.fshift = 6;
spte.stack2del = [];

% update variables
if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('FileName', 'var'); FileName = []; end
if ~exist('iparams', 'var'); iparams = []; end
spte = loparam_updater(spte, iparams);

% Selecting folders
f2run = dir;
f2run = str2match(FolderName, f2run);
f2run = str2rm(spte.fo2reject, f2run);
f2run = {f2run.name};

fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']); 
    cd(f2run{i});
    runperfolder(FileName, spte);
    cd(spte.cDir)
    fprintf('\n')
    
end

fprintf('... Done\n')

end

function runperfolder(fname, spte)
% runperfolder: run all files per folder
%
% Usage:
%   runperfolder(fname_trial)
%
% Args:
%   fname: file name pattern

% determine if fname narrows down to just one file
[f2plot, ~, rep2plot] = ...
    rdir_namesplit(fname, '.mat', ...
    spte.fsuffix, spte.fi2reject, [], 1);
fprintf('Temporal resampling & aligment of');

if numel(f2plot) == 1 && numel(rep2plot) == 1
    
    % run single file
    f2plot{1} = [f2plot{1}, '_', num2str(rep2plot)];
    fprintf(' one fly_seg ')  
    
else
    
    % Run all flies per folder
    f2plot = unique(f2plot);
    
end

fprintf([num2str(numel(f2plot)), ' flies\n'])

for file_i = 1:numel(f2plot)
    
    % doing resampling and aligning per fly in case there are repnum > 1
    [filename, ~, repnum] = ...
        rdir_namesplit(f2plot{file_i}, '.mat', ...
        spte.fsuffix, spte.fi2reject, [], 1);
    filename = unique(filename);
    repnum = sort(repnum);
    
    if numel(filename) == 1
        processfunc(filename{1}, repnum, spte);
    else
        fprintf('error')
    end
    
end

end

function processfunc(f2run, repnum, spte)
% processfunc: run all reps per f2run name
%
% Usage:
%   runperfolder(fname_trial)
%
% Args:
%   f2run: file name
%   repnum: reps to load and align

% load all trials or reps per fly to get the best spte.time
spte.time = getInitEnd(f2run, repnum, spte.newtimeres);
spte

% generate new timestamps
lstEn = [];
endT = spte.time(1):spte.newtimeres:spte.time(2);

fprintf(['Running fly : ', strrep(strrep(f2run, ['.', filesep], ''), '_', ' '), '\n'])

for i = 1:numel(repnum)
    
    % load metadata & data
    rep2run = [f2run, '_', num2str(repnum(i)), '_rawdata.mat']; 
    load(rep2run)
    
    if exist(strrep(rep2run, 'raw', 'ref'), 'file')
        Data_ref = load(strrep(rep2run, 'raw', 'ref'), 'Data'); 
        Data_ref = Data_ref.Data;
    end

    load(strrep(rep2run, '_rawdata', '_metadata'), ...
        'iDat', 'lStim', 'fDat')

    if ~isempty(strfind(fDat.DataType, 'opto'))
        load(strrep(rep2run, '_rawdata', '_metadata'), 'cDat')
    end
    
    % remove unwanted stacks
    
    if ~isempty(spte.stack2del)
        [lStim, iDat, Data, Data_ref] =  ...
            volumeprunner(lStim, iDat, spte.stack2del, Data, Data_ref);
    end
    
    % match X and Y resolution
    inpres = round(iDat.MetaData{3}*spte.fres)/spte.fres;
    
    if isempty(spte.newres)
        spte.newres = inpres;
    end
    
    if ~isfield(iDat, 'XYresmatch') || iDat.XYresmatch == 0
        
        % resample functional data
        % (it removes nans prior to resampling and then puts them back)
        fprintf(['iVoxelsize: ', num2str(iDat.MetaData{3}(1:3)), ' '])
        fprintf('XY green ');
        [Data, nan_idx, nan_idx_2_i] = ...
            spa_temp_resamp_int(Data, spte.art_val, ...
            inpres, spte.newres, [], []);
        
        if exist('Data_ref', 'var')
            fprintf('XY red ');
            [Data_ref, ~, ~] = ...
                spa_temp_resamp_int(Data_ref, spte.art_val, ...
                inpres, spte.newres, nan_idx, nan_idx_2_i);
        end
        
        % resample structural Data
        if ~isempty(iDat.RedChaMean)
            fprintf('XY redmean ');
            iDat.RedChaMean = interp3DxT(iDat.RedChaMean, ...
                inpres, spte.newres, 2);
        end
        
        if ~isempty(iDat.GreenChaMean)
            fprintf('XY greenmean ');
            iDat.GreenChaMean = interp3DxT(iDat.GreenChaMean, ...
                inpres, spte.newres, 2);
        end
        
        % update metadata
        iDat.MetaData{3}(1:3) = spte.newres;
        iDat.FrameSize = [size(Data, 1) size(Data, 2)];
        iDat.XYresmatch = 1;
        fprintf([' eVoxelsize: ', num2str(iDat.MetaData{3}(1:3)), ' ']);
        clear nan_idx nan_idx_2_i
        
    else
        fprintf('Already XY resampled ')
    end
    
    % do spatial smoothing
    if ~isempty(spte.sigma) && ~isempty(spte.size) && iDat.sSmooth == 0
        
        fprintf('XY smoothing ')
        
        % functional Data
        Data = imblur(Data, spte.sigma, spte.size, 2);
        % delete edges
        Data = Data(2:end-1, 2:end-1, :, :);
        
        if exist('Data_ref', 'var')
            Data_ref = imblur(Data_ref, spte.sigma, spte.size, 2);
            % delete edges
            Data_ref = Data_ref(2:end-1, 2:end-1, :, :);
        end
        
        % structural Data
        if ~isempty(iDat.RedChaMean)
            iDat.RedChaMean = iDat.RedChaMean(2:end-1, 2:end-1, :);
        end
        
        if ~isempty(iDat.GreenChaMean)
            iDat.GreenChaMean = iDat.GreenChaMean(2:end-1, 2:end-1, :);
        end
        
        iDat.FrameSize = [size(Data, 1) size(Data, 2)];
        iDat.sSmooth = 1;
        iDat.sSmoothpar = [spte.sigma, spte.size];
        fprintf('\n')
        
    else
        
        if isempty(spte.sigma) || isempty(spte.size)
            fprintf('Not XY-smoothing data ')
        else 
            fprintf('Already XY-smoothed ')
        end
        
    end
    
    % Report when whole slices are nan
    nan_data = isnan(Data);
    nan_data = squeeze(sum(sum(nan_data, 1), 2)) ...
        == size(Data, 1)*size(Data, 2);
    
    if sum(nan_data(:) == 1) ~= 0
        
        fprintf([' ', num2str(sum(nan_data(:) == 1)), ...
            ' whole slices are nan\n'])
        if spte.debug == 1
            figure();
            imagesc(nan_data);
            keyboard;
        end
        
    end
    
    clear nan_data
    
    % Get stim init-end
    if ~isfield(iDat, 'lstEn')
        
        if repnum(i) == repnum(1)
            
            % single lstEn for all trials /reps per fly name
            % convert from sampling points to seconds
            
            if ~isfield(lStim, 'lstEn') && contains(fDat.DataType, 'opto')
                lStim.lstEn = optostim_init_end(lStim.trace, lStim)/lStim.fs;
            end
            
            lstEn = (lStim.lstEn)/lStim.fs;
            
        end
        
        iDat.lstEn = lstEn - lstEn(1, 1);
        
    else
        
        lstEn = iDat.lstEn;
        fprintf('Already added lstEn to iDat ')
        
    end
    
    % resample temporally green channel (from plane to plane)
    iDatcopy = iDat;
    
    if ~isfield(iDat, 'tResample') || iDat.tResample == 0
        
        % get new time stamps
        iniT = reshape(iDat.fstEn(:, 2), [iDat.FrameN iDat.StackN])/lStim.fs;
        iniT = iniT - lstEn(1, 1);
        
        % chop stim trace
        itIdx = round((endT(1) + lstEn(1, 1))*lStim.fs);
        lStim.trace = lStim.trace(itIdx:end);
        
        % resample PMT_fscore
        if isfield(iDat, 'PMT_fscore')
            iDat.PMT_fscore = interp1(iniT(1, :), iDat.PMT_fscore, endT);
        end
        
        % show initial time and final re-slicing
        fprintf(['time_i: ', num2str([max(iniT(:, 1)), ...
            min(iniT(:, end))]), ' time_e: ', num2str(endT([1 end])), ' '])
        fprintf('Time green ');
        
        Data = interp3DxTixTj(Data, [1 1], [1 1], iniT, endT);
        
        iDat.tResample = 1;
        iDat.Tres = endT;
        iDat.StackN = length(endT);
        
    else
        
        fprintf('Already TixTj resampled ');
        iniT = [];
        
    end
    
    Data = single(Data);
    
    % Report when whole volumes are nan
    nan_data = isnan(Data);
    nan_data = squeeze(sum(sum(sum(nan_data, 1), 2), 3)) ...
        == size(Data, 1)*size(Data, 2)*size(Data, 3);
    
    if sum(nan_data(:) == 1) ~= 0
        
        fprintf([' ', num2str(sum(nan_data(:) == 1)), ...
            ' whole volumes are nan\n'])
        if spte.debug == 1
            figure();
            plot(nan_data);
            keyboard;
        end
        
    end
    
    clear nan_data
  
    % save extra metadata for Opto
    if contains(fDat.DataType, 'opto') && ~contains(fDat.DataType, 'prv')
        
        if ~exist('cDat', 'var')
            cDat.minInit = nan;
            cDat.minEnd = nan;
            cDat.CorType = 'LTM';
            cDat.LTM = [];
            fprintf('save empty cDat (opto without LED denoising)\n');
            save(strrep(rep2run, '_rawdata', '_metadata'), 'cDat', '-append');
        end
        
        lStim.sPars.led_mini = cDat.minInit;
        lStim.sPars.led_mine = cDat.minEnd;
        lStim.sPars.led_delta = eval(['median(cDat.', cDat.CorType, ')']);
        
    end
    
    saveIm(Data, rep2run);
    clear Data
    
    % update iDat and lStim
    save(strrep(rep2run, '_rawdata', '_metadata'), 'iDat', 'lStim', '-append');
    
    % save wDat
    if spte.wDat_gen
        save_wDat(strrep(rep2run, '_rawdata', '_metadata'), '3DxT', ...
            spte.direction, spte.bkgate, spte.fshift, spte.blowcap);
    end
    
    % resample temporally ref channel (from plane to plane)
    if exist('Data_ref', 'var')
        
        if ~isfield(iDatcopy, 'tResample') || iDatcopy.tResample == 0
            fprintf('Time red ');
            Data_ref = interp3DxTixTj(Data_ref, [1 1], [1 1], iniT, endT);
        end
        
        Data_ref = single(Data_ref);
        saveIm(Data_ref, strrep(rep2run, 'raw', 'ref'));
        clear iDatcopy Data_ref
        
    end
    
    fprintf('\n')
    clear iDat lStim iniT Data inpres nantest rep2run
    
end

end

function saveIm(Data, filename)
% saveIm: function for saving 3DxT data
%
% Usage:
%   saveIm(Data, filename)
%
% Args:
%   Data: 3DxT matrix
%   filename: name to use for saving

save(filename, 'Data', '-v7.3');

end

function start_end_time = getInitEnd(f2run, repnum, newtimeres)
% getInitEnd: get start and end of trials across segments
%
% Usage:
%   start_end_time = getInitEnd(f2run, repnum, newtimeres)
%
% Args:
%   f2run: file name
%   repnum: reps to load and align
%   newtimeres: new temporal resolution in seconds
% 
% Notes
% Check for already run files: 
% 1) automatically generate time boundaries or 2) and retrieve Tres and iniT

load([f2run, '_', num2str(repnum(1)), '_metadata.mat'], 'iDat')
check_gate = 0;

if isfield(iDat, 'tResample')
    check_gate = iDat.tResample;
end

if ~check_gate

    for i = 1:numel(repnum)
        load([f2run, '_', num2str(repnum(i)), '_metadata.mat'], 'iDat', 'lStim')
        iniT_t = (iDat.fstEn(:, 2)' - lStim.lstEn(1, 1))/lStim.fs;
        iniT_t = reshape(iniT_t, [iDat.FrameN iDat.StackN]);
        iniT(i, :) = [max(iniT_t(:, 1)), min(iniT_t(:, end))];
    end

    start_end_time(1) = ceil(max(iniT(:, 1))/newtimeres)*newtimeres;
    start_end_time(2) = floor(min(iniT(:, 2))/newtimeres)*newtimeres;

else

    fprintf('File already resampled copying time\n')
    start_end_time = iDat.Tres([1 end]);
    iniT = start_end_time;

end

fprintf(['Target time ends are: ', num2str(start_end_time(1)), ' ', ...
    num2str(start_end_time(2)), ...
    ' (automatically generated), original time ends: ', ...
    num2str(max(iniT(:, 1))), ' ', num2str(min(iniT(:, 2))), '\n']);

clear iniT iniT_t iDat lStim

end

function [iIm, nan_idx, nan_idx_2_i] = ...
    spa_temp_resamp_int(iIm, nan2val, ...
    oldres, newres, nan_idx, nan_idx_2_i)
% spa_temp_resamp_int: internal function that performs the resampling, 
% it replaces all nans by negative values then it puts them back
%
% Usage:
%   [iIm, nan_idx, nan_idx_2_i] = ...
%       spa_temp_resamp_int(iIm, nan2val, ...
%       oldres, newres, nan_idx, nan_idx_2_i)
%
% Args:
%   iIm: 3DxT matrix
%   nan2val: value to replace nans
%   oldres: resolution of input data
%   newres: target resolution of data
%   nan_idx: nan pixels before resampling
%   nan_idx_2_i: nan pixels after resampling

% get original size
Im_siz_i = size(iIm);

% from XYZT --> RT
iIm = reshape(iIm, prod(Im_siz_i(1:end-1)), Im_siz_i(end));

% get min val
iIm_min = min(iIm(:));

% get nan val
if ~exist('nan_idx', 'var') || isempty(nan_idx)
    nan_idx = max(isnan(iIm), [], 2);
end

% replace nans --> nan2val
iIm(nan_idx, :) = nan2val;

% from RT --> XYZT
iIm = reshape(iIm, Im_siz_i);

% interpolate XY
iIm = interp3DxT(iIm, oldres, newres, 2);

% gt new size
Im_siz_o = size(iIm);

% from XYZT --> RT
iIm = reshape(iIm, prod(Im_siz_o(1:end-1)), Im_siz_o(end));

if ~exist('nan_idx_2_i', 'var') || isempty(nan_idx)
    nan_idx_2_i = max(iIm < iIm_min, [], 2);
end

% replace nan2val --> nans
iIm(nan_idx_2_i, :) = nan;

% from RT --> XYZT
iIm = reshape(iIm, Im_siz_o);

end

function [lStim, iDat, Y1, Y2] = ...
    volumeprunner(lStim, iDat, stack2del, Y1, Y2)
% volumeprunner: prunes the original Data
%   and adjust lStim and iDat metadata accordingly
%
% Usage:
%   [lStim, iDat] =  volumeprunner(lStim, iDat, data_obj)
%
% Args:
%   lStim: stimuli metadata variable
%   iDat: image metadata variable
%   stack2del: timepoints to delete
%   Y1: 3DxT data
%   Y2: 3DxT data (another channel)

if ~exist('Y1', 'var'); Y1 = []; end
if ~exist('Y2', 'var'); Y2 = []; end

% Prune Data
if ~isempty(stack2del)
    
    % load data    
    if iDat.StackN > 1
        
        if ~isempty(Y1)
            Y1(:, :, :, stack2del, :) = [];
        end
        
        if ~isempty(Y2)
            Y2(:, :, :, stack2del, :) = [];
        end
        
    else
        
        if ~isempty(Y1)
            Y1(:, :, :, stack2del) = [];
        end
        
        if ~isempty(Y2)
            Y2(:, :, :, stack2del) = [];
        end
        
    end

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
    lStim.trace = lStim.trace(iDat.fstEn(1, 1):iDat.fstEn(end, 2));
    lStim.lstEn = lStim.lstEn - iDat.fstEn(1, 1) + 1;
    
    % Zero frame init
    iDat.fstEn = iDat.fstEn - iDat.fstEn(1, 1) + 1;
    
    % Get frame times and prune them too
    iDat.StackN = iDat.StackN - numel(stack2del);
    
    % Update sstEn
    if iDat.StackN > 1
        preInit = min(reshape(iDat.fstEn(:, 1), [iDat.FrameN, iDat.StackN]), [], 1)';
        preEnd =  max(reshape(iDat.fstEn(:, 1), [iDat.FrameN, iDat.StackN]), [], 1)';
        iDat.sstEn = [preInit, preEnd];
        clear preInit preEnd
    end
    
end

end
