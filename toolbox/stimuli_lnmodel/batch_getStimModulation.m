function batch_getStimModulation(FolderName, FileName, iparams)
% batch_getStimModulation: Use Linear model to fit calcium traces
%
% Usage:
%   batch_getStimModulation(FolderName, FileName, iparams)
%
% Args:
%   FolderName: folders to use
%   FileName: files to use
%   iparams: parameters
%       (cDir: current directory)
%       (fo2reject: folders to reject)
%       (fi2reject: files to reject)
%       (fsuffix: suffix of files with roi varaiable)
%           (default, '_prosroi')
%       (metsuffix: suffix of files with wDat variable)
%           (default, '_prosmetadata')
%       (field2load: metadata variable 'wDat' or 'iDat')
%           (default, 'wDat')
%       (redo: redo even if sMod exis, redo raw, redo shuffle)
%           (default, [0 0 0])
%       (varname: output variable name
%           (default, 'sMod')
%       %%%%%%%%%%%% filter related %%%%%%%%%%%%
%       (filterlength: filter length in seconds)
%           (default, 10)
%       (customtrial_init: number of seconds before stimuli 
%           start to count trial start (abs number))
%           (default, [])
%       (customtrial_end: number of seconds after stimuli 
%           end to count trial end (abs number))
%           (default, [])
%       (stim2use: which stimuli to use)
%           (default, [], or all)
%       (trials2use: trials to use, by default it 
%           uses the min number of trials across stimuli)
%           (default, [], or min number of trials across stimuli)
%       (minsize: minimun size of chunks in seconds)
%           (default, 5)
%       (per2use: percentage of trials to use for training [0-1])
%           (default, 0.8)
%       (type2run: type of shuffle to perform (as, af, cs)
%           (as: Random Shuffle)
%           (af: Amplitude Adjusted Fourier Transform)
%           (cs: Random chunk Shuffle)
%           (default, default [0 0 1])
%       (btn: number of permutations)
%           (default, 10^4)
%       %%%%%%%%%%%% parpool & server related %%%%%%%%%%%%
%       (serId: server id)
%           (default, 'int')
%       (corenum: number of cores)
%           (default, 4)
%       (chunksiz: number of chunks for parpool)
%           (default, 80)
%       (febgate: gate to generate filter error bounds)
%           (default, 0)
%       (oDir: output directory to save summary results)
%           (default, pwd)
%
% Notes:
% for info about regression algorithm see:
%   runRidgeOnly (ridge regression using empirical Bayes)
% for info about surrogate data see https://en.wikipedia.org/wiki/Surrogate_data_testing
% for info about crossvalidation see https://en.wikipedia.org/wiki/Cross-validation_(statistics)
% 
% 01/18/19:
%   add option to define trial start and end (customtrial_init, and customtrial_end fields)
%   add option to define trials to use (trials2use)

% Default params
pSM = []; 
pSM.cDir = pwd;
pSM.fo2reject = {'.', '..', 'preprocessed', ...
    'BData', 'rawtiff', 'motcor', 'stitch', ...
    'dfrel_vid', 'smod', 'roicov'}; 
pSM.fi2reject = {'Zstack'};
pSM.fsuffix = '_prosroi'; 
pSM.metsuffix = '_prosmetadata';
pSM.field2load = 'wDat'; 
pSM.redo = [0 0 0];
pSM.varname = 'sMod';
pSM.filterlength = 10;
pSM.customtrial_init = [];
pSM.customtrial_end = [];
pSM.stim2use = [];
pSM.trials2use = [];
pSM.minsize = 5;
pSM.per2use = 0.8;
pSM.type2run = [0 0 1];
pSM.btn = 10^4;
pSM.serverid = 'int';
pSM.corenum = 4;
pSM.chunksiz = 80;
pSM.febgate = 0;

pSMplot = [];
pSMplot.hbins = -1:0.01:1;
pSMplot.hbinsp = 0:0.01:1.1;
pSMplot.prct2use = 30;
pSMplot.fdr = 0.01;
pSMplot.mccor_method = 'dep';
pSMplot.dir_depth = 1;

% update variables
if ~exist('FileName', 'var'); FileName = []; end
if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('iparams', 'var'); iparams = []; end
pSM = loparam_updater(pSM, iparams);

pSMplot.fsuffix = [pSM.fsuffix, '.mat'];
pSMplot.fmetsuffix = [pSM.metsuffix, '.mat'];

% start pararell pool if not ready yet
ppobj = setup_parpool(pSM.serverid, pSM.corenum);

% Selecting folders
f2run = dir;
f2run = str2match(FolderName, f2run);
f2run = str2rm(pSM.fo2reject, f2run);
f2run = {f2run.name};
fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']); 
    
    cd(f2run{i});
    
    oDir = [pwd, filesep, 'smod'];
    
    runperfolder(FileName, pSM);
    cd(pSM.cDir);
    
    batch_plot_sMod_results(FileName, oDir, pSMplot)
    
    fprintf('\n')
    
end

delete_parpool(ppobj);

fprintf('... Done\n')

end

function runperfolder(FileName, pSM)
% runperfolder: run all files per folder
%
% Usage:
%   runperfolder(fname)
%
% Args:
%   fname: file name pattern
%   pSM: parameter variable

% Files to load
f2run = rdir(['.', filesep, '*', pSM.fsuffix, '*.mat']);
f2run = str2match(FileName, f2run); 
f2run = {f2run.name};

fprintf('Estimating stimuli modulation for all ROIs\n')
fprintf(['Running n-files : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running file : ', strrep(f2run{i}, filesep, ' '), '\n']);
    getStimModulation(f2run{i}, pSM);
    
end

fprintf('****** Done ******\n')

end

function getStimModulation(filename, pSM)
% getStimModulation: run stimuli modulation for each file
%
% Usage:
%   runperfile(f2run, pSM)
%
% Args:
%   f2run: file name
%   pSM: parameter variable

global pSM

t0 = stic;

% 1) load required varaibles 'roi', and 'pSM.field2load'
%   compatible with old format

load(filename, 'roi');
load(strrep(filename, pSM.fsuffix, pSM.metsuffix), ...
    pSM.field2load)

if exist('iDat', 'var')
    
    load(strrep(filename, pSM.fsuffix, pSM.metsuffix), 'fDat')
    wDat.sTime = iDat.lstEn;
    wDat.fTime = iDat.Tres;
    wDat.datatype = fDat.DataType;
    
end

% 2) convert filterlenght from seconds to timepoints

dt = wDat.fTime(2) - wDat.fTime(1);

if isempty(pSM.filterlength)
    
    % Automatically calculate the max lag based on 
    %   the length of the stimuli divided by the sampling rate
    filterlength_tp = round(round(wDat.sTime(1, 2) ...
        - wDat.sTime(1, 1))/dt);
    
else
    
    filterlength_tp = round(pSM.filterlength/dt);
    
end

sMod.stim = getStimVect(wDat); 
[K, ~] = size(roi.filtered.dfof);

fprintf(['ROIs to process: ', num2str(K), '\n'])

% 3) Get number of stimuli presented

stim_ntype = 1;

if ~contains(wDat.datatype, 'opto') || contains(wDat.datatype, 'prv')
    
    stim_name = cellfun(@(x, y) [x '_AMP_' num2str(y)], ...
        wDat.sPars.name(1, :), chunk2cell(wDat.sPars.int(1, :), 1), ...
        'UniformOutput', false);
    [stim_name, ~, stim_u_idx] = unique(stim_name, 'stable');
    stim_all_idx = stim_u_idx(wDat.sPars.order(1, :));
    stim_all_idx = stim_all_idx(1:size(wDat.sTime, 1), 1);
    
    clear stim_u_idx
    
    stim_count = hist(stim_all_idx, 1:1:max(stim_all_idx(:)));
    stim_ntype = numel(unique(stim_name));
    
end

% 4) Define stim to use for model fitting

if isempty(pSM.stim2use)
    pSM.stim2use = 1:stim_ntype;
end

% 4.1) use predefined trial
if isfield(wDat, 'trialsToUse') && ~isempty(wDat.trialsToUse) ...
    || isfield(wDat, 'trials2use') && ~isempty(wDat.trials2use)

    try pSM.trials2use = wDat.trialsToUse; end
    try pSM.trials2use = wDat.trials2use; end
    
    % update number of trials
    if ~isempty(pSM.trials2use)
        if sum(pSM.trials2use > min(stim_count)) == 0
            stim_count(:) = numel(pSM.trials2use);
        else
            fprintf('Error trials requested are outside range\n')
            return
        end
    end
    
end

% 5) Build matrix of stimuli used for model with all lags

% collect stimuli to use
StimE_M = [];

if stim_ntype ~= 1
    jj = 1;
    for j = pSM.stim2use
        StimE(jj, :) = getStimVect(wDat, 1, pSM.trials2use, j);
        jj = jj + 1;
    end
else
    StimE(1, :) = sMod.stim;
end

% build matrix with different lags (to match filter length)
for i = 1:size(StimE, 1)
    [StimE_temp{i, 1}, ~] = ...
        buffer(StimE(i,:), filterlength_tp, filterlength_tp - 1);
end

% compile StimE_M: t x tau x type
StimE_M = cell2mat(StimE_temp);
sMod.stimM = StimE_M';
clear StimE_M StimE jj i j

display(['Stimuli to use: ', num2str(pSM.stim2use)])
display(['Variable name to use: ', num2str(pSM.varname)])

% 6) Get stimuli and trials to use for fitting and testing
% discretize stim in trials:
stim_trial_idx = getTrialVect(wDat, pSM.stim2use, ...
    pSM.customtrial_init, pSM.customtrial_end, pSM.trials2use);

% get the minimum number of trials of pSM.stim2use
%   to determine number of trials for training data
if exist('stim_count', 'var')
    
    min_ntrial = min(stim_count(pSM.stim2use));
    ntrial_train = floor(min_ntrial*pSM.per2use);
    
    % update sMod.tidx2use
    idx2change = unique(stim_trial_idx, 'stable');
    idx2change(idx2change == 0) = [];
    stim_trial_idx_b = changem(stim_trial_idx, ...
        1:numel(idx2change), idx2change);
    
    sMod.tidx2use = find(stim_trial_idx ~= 0 & ...
        stim_trial_idx_b <= min_ntrial);
    clear min_ntrial stim_count
    clear stim_trial_idx_b idx2change
    
else
    
    sMod.tidx2use = find(stim_trial_idx ~= 0);
    ntrial_train = floor((size(wDat.sTime, 1)/stim_ntype)*pSM.per2use);
    
end

% prune un-used stims
stim_trial_idx = stim_trial_idx(sMod.tidx2use);
sMod.stimM = sMod.stimM(sMod.tidx2use, :);

% build matrix with indexes of timepoints used
%   as training data (testing data is ~train_idx)
trial_comb = nchoosek(unique(stim_trial_idx), ntrial_train);
train_idx = zeros(size(trial_comb, 1), size(stim_trial_idx, 2), 'logical');

for iter_i = 1:size(trial_comb, 1)
    for t_i = 1:size(trial_comb, 2)
       train_idx(iter_i, stim_trial_idx == trial_comb(iter_i, t_i)) = 1;
    end
end

clear trial_comb ntrial_train iter_i t_i

% 7) Automatically calculate the minsize and 
%   number of splits base on the trial length
% transform minsize to timepoints
chunk_minsize = pSM.minsize/dt;
chunk_splitn = floor(numel(sMod.tidx2use)/chunk_minsize);

clear stim_trial_idx

% 8) Collect Ca trace and Stim
% prune unused stims
catrace_in = roi.filtered.dfof(:, sMod.tidx2use);

stim_in = sMod.stimM; 
stimSiz = size(stim_in, 2);
stim_bin = sMod.stim(1, sMod.tidx2use);

stocf(t0, 'Time consumed so far: ')

% 9) using LN model to fit the traces, 
%   find the pearson correlation to modeled
%   trace and bootstrap the pearson correlation

% Run Raw data
fprintf('Running Raw Data\n')
vList = whos('-file', filename);
sgate = sum(strcmp({vList.name}, pSM.varname));
if sgate && ~pSM.redo(1)
    fprintf([pSM.varname ' already exist in mat file\n']); 
    return;
end
tgate = exist([strrep(filename, ...
    [pSM.fsuffix, '.mat'], ''), '_temp.mat'], 'file');

if ~tgate || pSM.redo(2)
    
    % calculate filters and pearson correlations for raw data
    [pcor_raw, lFilter] = runridgeraw_chunks(...
        train_idx, catrace_in, stim_in, stimSiz, ...
        pSM.chunksiz, pSM.corenum);
    save([strrep(filename, ...
        [pSM.fsuffix, '.mat'], ''), ...
        '_temp.mat'], 'pcor_raw', 'lFilter', '-v7.3')
    
else
    
    load([strrep(filename, ...
        [pSM.fsuffix, '.mat'], ''), ...
        '_temp.mat'], 'pcor_raw', 'lFilter')
    
end

stocf(t0, 'Time consumed so far: ')

% name of temporary file
filename_temp = [strrep(filename, [pSM.fsuffix, '.mat'], ''), ...
    '_temp.mat'];

% Run (as) Random Shuffle
if pSM.type2run(1)
    
    iFilter = squeeze(mean(lFilter, 2));
    
    fprintf('Running Phase # 1: Random shuffle\n'); 
    lgate = 1;
    
    if tgate
        vList = whos('-file', filename_temp);
        lgate = ~sum(contains({vList.name}, 'pcor_as'));
    end
    
    [pcor_as, lFilter_as_mean, ...
        lFilter_as_med, lFilter_as_sd] = ...
        runridgeshuffle(train_idx(1, :), ...
        catrace_in, stim_in, stimSiz, ...
        iFilter, lgate, stim_bin, pSM.redo(3), ...
        filename_temp, pSM.febgate, ...
        pSM.chunksiz, pSM.corenum, pSM.btn);
    % edit runridgeshuffle

    save(filename_temp, 'pcor_as', 'lFilter_as_mean', ...
        'lFilter_as_med', 'lFilter_as_sd', '-append')
    
end
stocf(t0, 'Time consumed so far: ')

% Run (af) Amplitude Adjusted Fourier Transform
if pSM.type2run(2)
    
    iFilter = squeeze(mean(lFilter, 2));
    
    fprintf('Running Phase # 2: Random Phases\n');
    lgate = 1;
    
    if tgate
        vList = whos('-file', filename_temp);
        lgate = ~sum(contains({vList.name}, 'pcor_afs'));
    end
    
    [pcor_afs, lFilter_afs_mean, ...
        lFilter_afs_med, lFilter_afs_sd] = ...
        runridgeshuffle_af(train_idx(1, :), ...
        catrace_in, stim_in, stimSiz, iFilter, ...
        lgate, stim_bin, pSM.redo(3), filename_temp, ...
        pSM.febgate, pSM.chunksiz, pSM.corenum, pSM.btn);
    % edit runridgeshuffle_af

    save(filename_temp, 'pcor_afs', 'lFilter_afs_mean', ...
        'lFilter_afs_med', 'lFilter_afs_sd', '-append')
    
end
stocf(t0, 'Time consumed so far: ')

% Run (cs) Random chunk Shuffle
if pSM.type2run(3)
    
    iFilter = squeeze(mean(lFilter, 2));

    fprintf('Running Phase # 3: Random shuffle in chunks\n')
    lgate = 1;
    
    if tgate
        vList = whos('-file', filename_temp);
        lgate = ~sum(contains({vList.name}, 'pcor_cs'));
    end
    
    [pcor_cs, lFilter_cs_mean, ...
        lFilter_cs_med, lFilter_cs_sd] = ...
        runridgeshuffle_chunks(train_idx(1, :), ...
        catrace_in, stim_in, stimSiz, iFilter, ...
        lgate, stim_bin, pSM.redo(3), filename_temp, ...
        pSM.febgate, chunk_minsize, chunk_splitn, ...
        pSM.chunksiz, pSM.corenum, pSM.btn);
    % edit runridgeshuffle_chunks
    
    save(filename_temp, 'pcor_cs', ...
        'lFilter_cs_mean', 'lFilter_cs_med', ...
        'lFilter_cs_sd', '-append')
    
end
stocf(t0, 'Time consumed so far: ')

% delete temp file:

delete([strrep(filename, [pSM.fsuffix, '.mat'], ''), '_temp.mat'])

% Save fields

% always updates this fields
sMod.maxlag = filterlength_tp; 
sMod.btn = pSM.btn; 
sMod.stim2use = pSM.stim2use;
sMod.CC_raw = pcor_raw; 
sMod.lFilter = lFilter;

if pSM.type2run(1)
    sMod.CC_as = pcor_as;
    sMod.lFilter_as_mean = lFilter_as_mean; 
    sMod.lFilter_as_med = lFilter_as_med;
    sMod.lFilter_as_sd = lFilter_as_sd; 
end

if pSM.type2run(2)
    sMod.CC_afs = pcor_afs;
    sMod.lFilter_afs_mean = lFilter_afs_mean; 
    sMod.lFilter_afs_med = lFilter_afs_med;
    sMod.lFilter_afs_sd = lFilter_afs_sd;
end

if pSM.type2run(3)
    sMod.CC_cs = pcor_cs;
    sMod.lFilter_cs_mean = lFilter_cs_mean; 
    sMod.lFilter_cs_med = lFilter_cs_med;
    sMod.lFilter_cs_sd = lFilter_cs_sd;
end 

% sMod.snr = std(CaPred, [], 2).^2./std(CaRaw - CaPred, [], 2).^2;

% save variable with the custom name
eval([pSM.varname, ' = sMod']);
save(filename, pSM.varname, '-append')

end

function [pcor_raw, lFilter] = runridgeraw_chunks(...
    train_idx, catrace_in, stim_in, stimSiz, ...
    chunksiz, corenum)
% runridgeraw_chunks: Run Ridge regression to all 
%   the possible combinations of train and test data
%
% Usage:
%   [pcor_raw, lFilter] = runridgeraw_chunks(...
%       train_idx, catrace_in, stim_in, stimSiz, K, ...
%       chunksiz, corenum)
%
% Args:
%   train_idx: indexes of train timepoints [n, T]
%       n: different train arragaments, T: time.
%   catrace_in: time traces [n, T]
%   stim_in: stimuli to use for prediction (at different lags) [T, m]
%       m: number of weights.
%   stimSiz: stimuli size (in case stimM is composed of many stimuli types)
%   chunksiz: number of chunks for parpool
%   corenum: number of cores
%
% Outputs:
%   pcor_raw: pearson correlation predicted vs raw for test indeces
%   lFilter: filter per train arragement

siz1 = size(catrace_in, 1);

[~, ~, chunk_idx] = ppool_makechunks(...
    chunksiz, corenum, siz1);

for i = 1:numel(chunk_idx)
    
    batch2run = chunk_idx{i};
    
    parfor ii = 1:numel(batch2run)
        
        roi_i_l{ii, 1} = (batch2run(ii):min(batch2run(ii) ...
            + chunksiz - 1, siz1));
        k = 1;
        t_lFilter = [];
        t_pcor_raw = [];
        
        for iii = roi_i_l{ii, 1}
            
            if i == 1 && ii == 1 && k == 1
                t0_ = stic;
            end
            
            [t_pcor_raw(k, :), t_lFilter(:, :, k)] = ...
                getsMod_ridgeperiter_raw(train_idx, ...
                catrace_in(iii, :), stim_in, stimSiz);

            if i == 1 && ii == 1 && k == 1
                fprintf(['time per roi ', ...
                    num2str(stoc(t0_)), ' seconds\n']); 
            end
            
            k = k + 1;
            
        end
        
        tb_pcor_raw{ii, 1} = t_pcor_raw; 
        tb_lFilter{ii, 1} = t_lFilter;
        
    end
    
    pcor_raw(cat(2, roi_i_l{:}), :) = ...
        cat(1, tb_pcor_raw{:});
    lFilter(:, :, cat(2, roi_i_l{:})) = ...
        cat(3, tb_lFilter{:});
    
    clear tb_pcor_raw tb_lFilter roi_i_l
    
    if mod(i, 1) == 0
        fprintf('%2.1f%% of chunks completed \n', ...
            i*100/numel(chunk_idx));
    end
    
end

end

function [pcor_cs, lFilter_cs_mean, ...
    lFilter_cs_med, lFilter_cs_sd] = ...
    runridgeshuffle_chunks(train_idx, ...
    catrace_in, stim_in, stimSiz, ...
    iFilter, lgate, stim_bin, redo, ...
    filename, filterebound, chunk_minsize, ...
    chunk_splitn, chunksiz, corenum, btn)
% runridgeshuffle_chunks: Run Ridge regression to all 
%   the possible combinations of train and test data
%
% Usage:
%   [pcor_cs, lFilter_cs_mean, ...
%       lFilter_cs_med, lFilter_cs_sd] = ...
%       runridgeshuffle_chunks(train_idx, ...
%       catrace_in, stim_in, stimSiz, ...
%       iFilter, lgate, stim_bin, redo, ...
%       filename, filterebound, ...
%       chunk_minsize, chunk_splitn)
%
% Args:
%   train_idx: indexes of train timepoints [n, T]
%       n: different train arragaments, T: time.
%   catrace_in: time traces [n, T]
%   stim_in: stimuli to use for prediction (at different lags) [T, m]
%       m: number of weights.
%   stimSiz: stimuli size (in case stimM is composed of many stimuli types)
%   iFilter: input filter to use for prediction (LN_pcor) [m, 1],
%       m: number of weights
%   lgate: flag to generate permuted data
%   stim_bin: binary vector of stimuli
%   redo: redo flag
%   filename: name of temporary file to generate
%   filterebound: filter error bounds 
%       (get filter for shuffle data, only collects mean, median and sd)
%   chunk_minsize: minimun size of time chunk
%   chunk_splitn: number of times to split catrace_in in the time domain
%   chunksiz: number of chunks for parpool
%   corenum: number of cores
%   btn: number of permutations
%
% Outputs:
%   pcor_cs: pearson correlation predicted vs raw for test indecex
%   lFilter_cs_mean: mean estimated filter
%   lFilter_cs_med: median estimated filter
%   lFilter_cs_sd: std of estimated filters
%
% See also getsMod_ridgeperiter_ashuffle

siz1 = size(catrace_in, 1);
dataObj = matfile(filename, 'Writable', true);

if lgate || redo
    % initialize parameters
    dataObj.pcor_cs = [];
    dataObj.lFilter_cs_mean = [];
    dataObj.lFilter_cs_med = [];
    dataObj.lFilter_cs_sd = [];

    % Get shuffle time
    [~, rperm_chunkIdx] = randchunkper(...
        catrace_in(1, :), chunk_splitn, ...
        round(chunk_minsize), btn, stim_bin);
    dataObj.rperm_chunkIdx = rperm_chunkIdx;

    % make chunks to run
    [~, ~, chunk_idx] = ppool_makechunks(...
        chunksiz, corenum, siz1);
else
    % update initial index
    vect_init = size(dataObj.pcor_cs, 1) + 1;

    % Get shuffle time
    rperm_chunkIdx = dataObj.rperm_chunkIdx;

    % make chunks to run
    [~, ~, chunk_idx] = ppool_makechunks(...
        chunksiz, corenum, siz1, vect_init);       
end

% new implementation of shuffle with
%   parpool saving file in a regular basis
train_idx_l = train_idx(1, :);

for i = 1:numel(chunk_idx)
    
    if i == 1
        t0 = stic;
    end
    
    batch2run = chunk_idx{i};
    
    parfor ii = 1:numel(batch2run)
        
        roi_i_l{ii, 1} = (batch2run(ii):min(batch2run(ii)...
            + chunksiz - 1, siz1));
        k = 1;
        t_pcor_cs_i = [];
        t_fmean_cs_i = [];
        t_fmed_cs_i = [];
        t_fsd_cs_i = [];
        
        for iii = roi_i_l{ii, 1}
            
            if i == 1 && ii == 1 && k == 1
                t0_ = stic;
            end

            [t_pcor_cs_i(k, :), t_fmean_cs_i(:, k), ...
                t_fmed_cs_i(:, k), t_fsd_cs_i(:, k)] = ...
                getsMod_ridgeperiter_ashuffle(...
                train_idx_l, catrace_in(iii, :), stim_in, ...
                stimSiz, rperm_chunkIdx, iFilter(:, iii), ...
                filterebound)

            if i == 1 && ii == 1 && k == 1
                fprintf(['time per roi', ...
                    num2str(stoc(t0_)), ' seconds\n']); 
                fprintf(['Estimated time ', ...
                    num2str(numel(chunk_idx)*chunksiz*stoc(t0_)/3600), ...
                    ' hours\n']); 
            end
        
            k = k + 1;
       end
        
        tb_pcor_cs_i{ii, 1} = t_pcor_cs_i;
        tb_fmean_cs_i{ii, 1} = t_fmean_cs_i; 
        tb_fmed_cs_i{ii, 1} = t_fmed_cs_i;
        tb_fsd_cs_i{ii, 1} = t_fsd_cs_i;
        
    end
    
    dataObj.pcor_cs(cat(2, roi_i_l{:}), 1:btn) = ...
        cat(1, tb_pcor_cs_i{:});
    dataObj.lFilter_cs_mean(1:stimSiz, cat(2, roi_i_l{:})) = ...
        cat(2, tb_fmean_cs_i{:});
    dataObj.lFilter_cs_med(1:stimSiz, cat(2, roi_i_l{:})) = ...
        cat(2, tb_fmed_cs_i{:}); 
    dataObj.lFilter_cs_sd(1:stimSiz, cat(2, roi_i_l{:})) = ...
        cat(2, tb_fsd_cs_i{:});
    
    clear tb_pcor_cs_i roi_i_l tb_fmean_cs_i ...
        tb_fmed_cs_i tb_fsd_cs_i
    
    if i == 1
        fprintf(['Estimated time ', ...
            num2str(numel(chunk_idx)*stoc(t0)/3600), ' hours\n']);
    end
    
    if mod(i, 1) == 0
        fprintf('%2.1f%% of chunks completed \n', ...
            i*100/numel(chunk_idx));
    end
    
end

pcor_cs = dataObj.pcor_cs;
lFilter_cs_mean = dataObj.lFilter_cs_mean;
lFilter_cs_med = dataObj.lFilter_cs_med;
lFilter_cs_sd = dataObj.lFilter_cs_sd;

end

function [pcor_as, lFilter_as_mean, ...
    lFilter_as_med, lFilter_as_sd] = ...
    runridgeshuffle(train_idx, catrace_in, ...
    stim_in, stimSiz, iFilter, lgate, ...
    stim_bin, redo, filename, filterebound, ...
    chunksiz, corenum, btn)
% runridgeshuffle: Run Ridge regression to all 
%   the possible combinations of train and test data
%
% Usage:
%   [pcor_as, lFilter_as_mean, ...
%       lFilter_as_med, lFilter_as_sd] = ...
%       runridgeshuffle(train_idx, catrace_in, ...
%       stim_in, stimSiz, iFilter, lgate, ...
%       stim_bin, redo, filename, filterebound, ...
%       chunksiz, corenum, btn)
%
% Args:
%   train_idx: indexes of train timepoints [n, T]
%       n: different train arragaments, T: time.
%   catrace_in: time traces [n, T]
%   stim_in: stimuli to use for prediction (at different lags) [T, m]
%       m: number of weights.
%   stimSiz: stimuli size (in case stimM is composed of many stimuli types)
%   iFilter: input filter to use for prediction (LN_pcor) [m, 1],
%       m: number of weights
%   lgate: flag to generate permuted data
%   stim_bin: binary vector of stimuli
%   redo: redo flag
%   filename: name of temporary file to generate
%   filterebound: filter error bounds 
%       (get filter for shuffle data, only collects mean, median and sd)
%   chunk_minsize: minimun size of time chunk
%   chunk_splitn: number of times to split catrace_in in the time domain
%   chunksiz: number of chunks for parpool
%   corenum: number of cores
%   btn: number of permutations
%
% Outputs:
%   pcor_as: pearson correlation predicted vs raw for test indecex
%   lFilter_as_mean: mean estimated filter
%   lFilter_as_med: median estimated filter
%   lFilter_as_sd: std of estimated filters
%
% See also getsMod_ridgeperiter_ashuffle

siz1 = size(catrace_in, 1);
T = size(catrace_in, 2);

dataObj = matfile(filename, ...
    'Writable', true);

if lgate || redo
    % initialize parameters
    dataObj.pcor_as = [];
    dataObj.lFilter_as_mean = [];
    dataObj.lFilter_as_med= [];
    dataObj.lFilter_as_sd = [];

    % Get shuffle time
    rperm_chunkIdx = randshuffleper(T, btn, stim_bin);
    dataObj.rperm_chunkIdx = rperm_chunkIdx;

    % make chunks to run
    [~, ~, chunk_idx] = ...
        ppool_makechunks(chunksiz, corenum, siz1);
    
else
    
    % update initial index
    vect_init = size(dataObj.pcor_as, 1) + 1;

    % Get shuffle time
    rperm_chunkIdx = dataObj.rperm_chunkIdx;

    % make chunks to run
    [~, ~, chunk_idx] = ...
        ppool_makechunks(chunksiz, corenum, siz1, vect_init);
    
end

% new implementation of shuffle with parpool
%   saving file in a regular basis
train_idx_l = train_idx(1, :);

for i = 1:numel(chunk_idx)
    
    if i == 1
        t0 = stic;
    end
    batch2run = chunk_idx{i};
    
    parfor ii = 1:numel(batch2run)
        
        roi_i_l{ii, 1} = (batch2run(ii):min(batch2run(ii)...
            + chunksiz - 1, siz1));
        k = 1;
        t_pcor_as_i = [];
        t_fmean_as_i = [];
        t_fmed_as_i = [];
        t_fsd_as_i = [];
        
        for iii = roi_i_l{ii, 1}
            
            if i == 1 && ii == 1 && k == 1
                t0_ = stic;
            end

            [t_pcor_as_i(k, :), t_fmean_as_i(:, k), ...
                t_fmed_as_i(:, k), t_fsd_as_i(:, k)] = ...
                getsMod_ridgeperiter_ashuffle(train_idx_l, ...
                catrace_in(iii, :), stim_in, stimSiz, ...
                rperm_chunkIdx, iFilter(:, iii), filterebound);
            
            if i == 1 && ii == 1 && k == 1
                fprintf(['time per roi', ...
                    num2str(stoc(t0_)), ' seconds\n']); 
                fprintf(['Estimated time ', ...
                    num2str(numel(chunk_idx)*chunksiz*stoc(t0_)/3600), ...
                    ' hours\n']); 
            end
            
            k = k + 1;
            
        end
        
        tb_pcor_as_i{ii, 1} = t_pcor_as_i;
        tb_fmean_as_i{ii, 1} = t_fmean_as_i; 
        tb_fmed_as_i{ii, 1} = t_fmed_as_i;
        tb_fsd_as_i{ii, 1} = t_fsd_as_i;
        
    end
    
    dataObj.pcor_as(cat(2, roi_i_l{:}), 1:btn) = ...
        cat(1, tb_pcor_as_i{:});
    dataObj.lFilter_as_mean(1:stimSiz, cat(2, roi_i_l{:})) = ...
        cat(2, tb_fmean_as_i{:});
    dataObj.lFilter_as_med(1:stimSiz, cat(2, roi_i_l{:})) = ...
        cat(2, tb_fmed_as_i{:}); 
    dataObj.lFilter_as_sd(1:stimSiz, cat(2, roi_i_l{:})) = ...
        cat(2, tb_fsd_as_i{:});
    
    clear tb_pcor_as_i roi_i_l ...
        tb_fmean_as_i tb_fmed_as_i tb_fsd_as_i
    
    if i == 1
        fprintf(['Estimated time ', ...
            num2str(numel(chunk_idx)*stoc(t0)/3600), ...
            ' hours\n']);
    end
    
    if mod(i, 1) == 0
        fprintf('%2.1f%% of chunks completed \n', ...
            i*100/numel(chunk_idx));
    end
    
end

pcor_as = dataObj.pcor_as;
lFilter_as_mean = dataObj.lFilter_as_mean;
lFilter_as_med = dataObj.lFilter_as_med;
lFilter_as_sd = dataObj.lFilter_as_sd;

end

function [pcor_afs, lFilter_afs_mean, ...
    lFilter_afs_med, lFilter_afs_sd] = ...
    runridgeshuffle_af(train_idx, ...
    catrace_in, stim_in, stimSiz, iFilter, ...
    lgate, stim_bin, redo, filename, filterebound, ...
    chunksiz, corenum, btn)
% runridgeshuffle: Run Ridge regression to all 
%   the possible combinations of train and test data
%
% Usage:
%   [pcor_afs, lFilter_afs_mean, ...
%   	lFilter_afs_med, lFilter_afs_sd] = ...
%       runridgeshuffle_af(train_idx, ...
%       catrace_in, stim_in, stimSiz, iFilter, ...
%       lgate, stim_bin, redo, filename, filterebound, ...
%       chunksiz, corenum, btn)
%
% Args:
%   train_idx: indexes of train timepoints [n, T]
%       n: different train arragaments, T: time.
%   catrace_in: time traces [n, T]
%   stim_in: stimuli to use for prediction (at different lags) [T, m]
%       m: number of weights.
%   stimSiz: stimuli size (in case stimM is composed of many stimuli types)
%   iFilter: input filter to use for prediction (LN_pcor) [m, 1],
%       m: number of weights
%   lgate: flag to generate permuted data
%   stim_bin: binary vector of stimuli
%   redo: redo flag
%   filename: name of temporary file to generate
%   filterebound: filter error bounds 
%       (get filter for shuffle data, only collects mean, median and sd)
%   chunksiz: number of chunks for parpool
%   corenum: number of cores
%   btn: number of permutations
%
% Outputs:
%   pcor_afs: pearson correlation predicted vs raw for test indecex
%   lFilter_afs_mean: mean estimated filter
%   lFilter_afs_med: median estimated filter
%   lFilter_afs_sd: std of estimated filters
%
% See also getsMod_ridgeperiter_afshuffle

K = size(catrace_in, 1);
dataObj = matfile(filename, 'Writable', true);

if lgate || redo
    
    % initialize parameters
    dataObj.pcor_afs = [];
    dataObj.lFilter_afs_mean = [];
    dataObj.lFilter_afs_med = [];
    dataObj.lFilter_afs_sd = [];

    % make chunks to run
    [~, ~, chunk_idx] = ...
        ppool_makechunks(chunksiz, ...
        corenum, K);
    
else
    
    % update initial index
    vect_init = size(dataObj.pcor_afs, 1) + 1;

    % make chunks to run
    [~, ~, chunk_idx] = ...
        ppool_makechunks(chunksiz, ...
        corenum, K, vect_init);  
    
end

% new implementation of shuffle with parpool
%   saving file in a regular bases
train_idx_l = train_idx(1, :);

for i = 1:numel(chunk_idx)
    
    if i == 1
        t0 = stic;
    end
    batch2run = chunk_idx{i};
    
    parfor ii = 1:numel(batch2run)
        
        roi_i_l{ii, 1} = (batch2run(ii):min(batch2run(ii) ...
            + chunksiz - 1, K));
        k = 1;
        t_pcor_afs_i = [];
        t_fmean_afs_i = [];
        t_fmed_afs_i = [];
        t_fsd_afs_i = [];
        
        for iii = roi_i_l{ii, 1}
            
            if i == 1 && ii == 1 && k == 1
                t0_ = stic;
            end
            
            [t_pcor_afs_i(k, :), t_fmean_afs_i(:, k), ...
                t_fmed_afs_i(:, k), t_fsd_afs_i(:, k)] = ...
                getsMod_ridgeperiter_afshuffle(train_idx_l, ...
                catrace_in(iii, :), stim_in, stimSiz, ...
                iFilter(:, iii), btn, stim_bin, filterebound);
            
            if i == 1 && ii == 1 && k == 1
                fprintf(['time per roi', ...
                    num2str(stoc(t0_)), ' seconds\n']); 
                fprintf(['Estimated time ', ...
                    num2str(numel(chunk_idx)*chunksiz*stoc(t0_)/3600), ...
                    ' hours\n']); 
            end
            
            k = k + 1;
            
        end
        
        tb_pcor_afs_i{ii, 1} = t_pcor_afs_i;
        tb_fmean_afs_i{ii, 1} = t_fmean_afs_i; 
        tb_fmed_afs_i{ii, 1} = t_fmed_afs_i;
        tb_fsd_afs_i{ii, 1} = t_fsd_afs_i;
        
    end
    
    dataObj.pcor_afs(cat(2, roi_i_l{:}), 1:btn) = ...
        cat(1, tb_pcor_afs_i{:});
    dataObj.lFilter_afs_mean(1:stimSiz, cat(2, roi_i_l{:})) = ...
        cat(2, tb_fmean_afs_i{:});
    dataObj.lFilter_afs_med(1:stimSiz, cat(2, roi_i_l{:})) = ...
        cat(2, tb_fmed_afs_i{:}); 
    dataObj.lFilter_afs_sd(1:stimSiz, cat(2, roi_i_l{:})) = ...
        cat(2, tb_fsd_afs_i{:});
    
    clear tb_pcor_afs_i roi_i_l tb_fmean_afs_i ...
        tb_fmed_afs_i tb_fsd_afs_i
    
    if i == 1
        fprintf(['Estimated time ', ...
            num2str(numel(chunk_idx)*stoc(t0)/3600), ...
            ' hours\n']);
    end
    
    if mod(i, 1) == 0
        fprintf('%2.1f%% of chunks completed \n', ...
            i*100/numel(chunk_idx));
    end
    
end

pcor_afs = dataObj.pcor_afs;
lFilter_afs_mean = dataObj.lFilter_afs_mean;
lFilter_afs_med = dataObj.lFilter_afs_med;
lFilter_afs_sd = dataObj.lFilter_afs_sd;

end
