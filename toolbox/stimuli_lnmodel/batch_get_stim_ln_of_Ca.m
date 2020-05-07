function batch_get_stim_ln_of_Ca(FolderName, FileName, iparams)
% batch_get_stim_ln_of_Ca: Use Linear model to fit calcium traces
%   using stimuli as regressors
%
% Usage:
%   batch_get_stim_ln_of_Ca(FolderName, FileName, iparams)
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
%       (redo: redo even if sti_ln_ca exist, redo raw, redo shuffle)
%           (default, [0 0 0])
%       (varname: output variable name
%           (default, 'sti_ln_ca')
%       %%%%%%%%%%%% filter related %%%%%%%%%%%%
%       (filterlength_stim: filter length in seconds (internally calculated))
%           (default, [])
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
%       (btn: number of permutations)
%           (default, 10^4)
%       (btn: number of permutations)
%           (default, 10^4)
%       (pst_time: peri stimulus time to use to extract trial structure from whole trace (seconds))
%           (default, [-10 10])
%       %%%%%%%%%%%% significance related %%%%%%%%%%%%
%       (fdr: false discovery rate)
%            (default, 0.01)
%       (pval_cor_method: multiple comparison correction to use: dep, pdep, bh)
%            (default, 'dep')
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
% To do:
%   confirm it works with opto
%   make sure it is robust to randomized stimuli (specially the stim2use variable)

% Default params
pSM = []; 
pSM.cDir = pwd;
pSM.fo2reject = {'.', '..', 'preprocessed', ...
    'BData', 'rawtiff', 'motcor', 'stitch', ...
    'dfrel_vid', 'sti_ln_ca', 'sti_mot_ln_ca', ...
    'sti_ln_mot', 'roicov'}; 
pSM.fi2reject = {'Zstack'};
pSM.fsuffix = '_prosroi'; 
pSM.metsuffix = '_prosmetadata';
pSM.field2load = 'wDat'; 
pSM.redo = [0 0 0];
pSM.varname = 'sti_ln_ca';
pSM.filterlength_stim = [];
pSM.customtrial_init = [];
pSM.customtrial_end = [];
pSM.stim2use = [];
pSM.trials2use = [];
pSM.minsize = 5;
pSM.per2use = 0.8;
pSM.btn = 10^4;
pSM.pst_time = [-10 10];
pSM.fdr = 0.01;
pSM.pval_cor_method = 'dep';
pSM.serverid = 'int';
pSM.corenum = 4;
pSM.chunksiz = 80;
pSM.febgate = 0;

pSM_plot = [];
pSM_plot.hbins = -1:0.1:1;
pSM_plot.hbinsp = 0:0.01:1.1;
pSM_plot.prct2use = 30;
pSM_plot.fdr = 0.01;
pSM_plot.mccor_method = 'dep';
pSM_plot.dir_depth = 1;
pSM_plot.var2load = 'sti_ln_ca';
pSM_plot.coeffname = 'pearson correlation';
pSM_plot.coeffrange_im = {[0 .12], [0 .4]};
pSM_plot.coeffths = 0.1;

% update variables
if ~exist('FileName', 'var'); FileName = []; end
if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('iparams', 'var'); iparams = []; end
pSM = loparam_updater(pSM, iparams);

pSM_plot.fsuffix = [pSM.fsuffix, '.mat'];
pSM_plot.fmetsuffix = [pSM.metsuffix, '.mat'];

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
    
    oDir = [pwd, filesep, 'sti_ln_ca'];
    
    runperfolder(FileName, pSM);
    cd(pSM.cDir);
    
    batch_plot_stim_ln_of_Ca_results(f2run{i}, FileName, oDir, pSM_plot)
    
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
%   pSM: internal variable

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
%   pSM: internal variable

t0 = stic;

% 1) load required variables 'roi', and 'pSM.field2load'
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

% 2) get dt in seconds

dt = wDat.fTime(2) - wDat.fTime(1);

% get number or ROIs
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
        
    stim_count = hist(stim_all_idx, 1:1:max(stim_all_idx(:)));
    stim_ntype = numel(unique(stim_name));
    
    % get stimulus post baseline time to define filter length
    if ~isempty(pSM.stim2use)
        stim_u_idx = pSM.stim2use;
    end
    
    for i = 1:max(stim_u_idx)
        stim_post_lag(i, 1) = min(wDat.sPars.basPost(stim_u_idx == i));
    end
    
    % get min stimuli witdh
    stim_width = wDat.sTime(:, 2) - wDat.sTime(:, 1);
    stim_width = stim_width(ismember(stim_all_idx, stim_u_idx));
    stim_width = ceil(min(stim_width)*1.1);
    
    clear stim_u_idx stim_name i
    
end

% 4) use predefined trial
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

% 5) generate filter length for stim (by default pre stimuli is 0)
if isempty(pSM.filterlength_stim)
    
    % Automatically calculate the max lag based on 
    %   the length of the stimuli divided by the sampling rate
    filterlength_tp(:, 2) = round(round(stim_post_lag)/dt);
    filterlength_tp(:, 1) = 0;
    
else
    
    filterlength_tp = round(pSM.filterlength_stim/dt);
    filterlength_tp = repmat(filterlength_tp, [stim_ntype, 1]);
    
end

% 6) Define stim to use for model fitting
if isempty(pSM.stim2use)
    pSM.stim2use = 1:stim_ntype;
end

% 7) Build matrix of stimuli used for model with all lags
%   get a binary vector with stimuli presentation
sti_ln_ca.stim = getStimVect(wDat); 

% generate binary vectors for each stimulus to use pSM.stim2use
if stim_ntype ~= 1
    jj = 1;
    for j = pSM.stim2use
        stim_m(jj, :) = getStimVect(wDat, 1, pSM.trials2use, j);
        jj = jj + 1;
    end
else
    stim_m(1, :) = sti_ln_ca.stim;
end

% build stimuli matrix with different lags (to match filter length)
[sti_ln_ca.stimM, ~] = build_matrix_with_lags(stim_m, filterlength_tp, 0);
% ca.stimM: t x tau x type
clear  stim_m jj i j

display(['Stimuli to use: ', num2str(pSM.stim2use)])
display(['Variable name to use: ', num2str(pSM.varname)])

% 8) discretize timepoints into trials:
stim_trial_idx = getTrialVect(wDat, pSM.stim2use, ...
    pSM.customtrial_init, pSM.customtrial_end, pSM.trials2use);

% get the minimum number of trials of pSLM.stim2use
%   to determine number of trials for training data
if exist('stim_count', 'var')
    
    min_ntrial = min(stim_count(pSM.stim2use));
    ntrial_train = floor(min_ntrial*pSM.per2use);
    
    % update sti_ln_ca.tidx2use
    idx2change = unique(stim_trial_idx, 'stable');
    idx2change(idx2change == 0) = [];
    stim_trial_idx_b = changem(stim_trial_idx, ...
        1:numel(idx2change), idx2change);
    
    sti_ln_ca.tidx2use = find(stim_trial_idx ~= 0 & ...
        stim_trial_idx_b <= min_ntrial);
    clear min_ntrial stim_count
    clear stim_trial_idx_b idx2change
    
else
    
    sti_ln_ca.tidx2use = find(stim_trial_idx ~= 0);
    ntrial_train = floor((size(wDat.sTime, 1)/stim_ntype)*pSM.per2use);
    
end

% prune un-used stims
stim_trial_idx = stim_trial_idx(sti_ln_ca.tidx2use);
sti_ln_ca.stimM = sti_ln_ca.stimM(sti_ln_ca.tidx2use, :);

% 9) define trials for training and testing
%   build matrix with indexes of timepoints used
%   as training data (testing data is ~train_idx)
trialidx_u = unique(stim_trial_idx);
test_comb = chunk2cell(randperm (numel(trialidx_u)), numel(trialidx_u) - ntrial_train)';
trial_comb = cf(@(x) setdiff(trialidx_u, x), test_comb);
train_idx = zeros(size(trial_comb, 1), size(stim_trial_idx, 2), 'logical');

for iter_i = 1:size(trial_comb)
    for t_i = 1:length(trial_comb{iter_i})
       train_idx(iter_i, stim_trial_idx == trial_comb{iter_i}(t_i)) = 1;
    end
end

clear trial_comb test_comb ...
    ntrial_train iter_i t_i trialidx_u

% 10) Automatically calculate the minsize and 
%   number of splits base on the trial length
% transform minsize to timepoints
chunk_minsize = pSM.minsize/dt;
chunk_splitn = floor(numel(sti_ln_ca.tidx2use)/chunk_minsize);

clear stim_trial_idx

% 11) Collect Ca trace and Stim
% prune unused stims
Ca_in = roi.filtered.dfof(:, sti_ln_ca.tidx2use);

stocf(t0, 'Time consumed so far: ')

% 12) using LN model to fit Ca traces

% Run Raw data
fprintf('Running Raw Data\n')
vList = whos('-file', filename);
sgate = sum(strcmp({vList.name}, pSM.varname));
if sgate && ~pSM.redo(1)
    fprintf([pSM.varname ' already exist in mat file\n']); 
    return;
end

% name of temporary file
filename_temp = [strrep(filename, [pSM.fsuffix, '.mat'], ''), ...
    '_temp.mat'];
tgate = exist(filename_temp, 'file');

if ~tgate || pSM.redo(2)
    
    % calculate filters and predicted traces for Ca
    [~, lFilter, Ca_pred] = ridgeregres_per_row_crossval(...
        train_idx, Ca_in, sti_ln_ca.stimM, pSM.chunksiz, pSM.corenum, 'bayes');
    pcor_raw = diag(corr(zscorebigmem(Ca_in)', zscorebigmem(Ca_pred)'));
    save(filename_temp, 'pcor_raw', 'lFilter', 'Ca_pred', '-v7.3')
    
else
    
    load(filename_temp, 'pcor_raw', 'lFilter', 'Ca_pred')
    
end

stocf(t0, 'Time consumed so far: ')

% 13) testing model significance
% Run (cs) Random chunk Shuffle

fprintf('Running Phase # 3: Random shuffle in chunks\n')
lgate = 1;

if tgate
    vList = whos('-file', filename_temp);
    lgate = ~sum(contains({vList.name}, 'pcor_cs'));
end

add_stim_lag = [-ceil(4/dt) ceil(stim_width/dt)];
shuffle2use = 2;
fprintf(['adding stimuli lag of: ', num2str(add_stim_lag), ' (timestamps) \n'])

[~, pcor_cs] = get_explained_variance_shuffle(...
    Ca_in, Ca_pred, lgate, sti_ln_ca.stim(1, sti_ln_ca.tidx2use), ...
    pSM.redo(3), filename_temp, ...
    chunk_minsize, chunk_splitn, ...
    pSM.chunksiz, pSM.corenum, pSM.btn, ...
    shuffle2use, add_stim_lag, []);

save(filename_temp, 'pcor_cs', '-append')
    
stocf(t0, 'Time consumed so far: ')

% get p-values
sti_ln_ca.pval = calculate_pval_2(pcor_raw, ...
    pcor_cs, pSM.fdr, pSM.pval_cor_method);

% 14) calculating ceiling of explained variance using mean-per-trial
fprintf('Split trace into trials\n')
wDat_edit = wDat;
wDat_edit.fTime = wDat.fTime(sti_ln_ca.tidx2use);
[roi_per_trial, ~, ~, ~, ~, ~, ~, ~, ~, ~] = trace2trials(wDat_edit, ...
    zscorebigmem(Ca_in), pSM.pst_time, stim_all_idx, 1, []);
stocf(t0, 'Time consumed so far: ')

sti_ln_ca.eVar_mean = cell2mat(cf(@(x, y) corr(horz(x')', horz(y')').^2, ...
    roi_per_trial, cf(@(x) repmat(mean(x, 1), [size(x, 1), 1]), roi_per_trial)));

% 15) calculating ceiling of explained variance using median-per-trial
sti_ln_ca.eVar_med = cell2mat(cf(@(x, y) corr(horz(x')', horz(y')').^2, ...
    roi_per_trial, cf(@(x) repmat(median(x, 1), [size(x, 1), 1]), roi_per_trial)));

% delete temp file:
delete(filename_temp)

% 16) Save fields
sti_ln_ca.maxlag = filterlength_tp;
sti_ln_ca.btn = pSM.btn;
sti_ln_ca.stim2use = pSM.stim2use;
sti_ln_ca.CC_raw = pcor_raw;
sti_ln_ca.lFilter = lFilter;
sti_ln_ca.CC_cs = pcor_cs;
sti_ln_ca.train_idx = train_idx;

% save variable with the custom name
eval([pSM.varname, ' = sti_ln_ca']);
save(filename, pSM.varname, '-append')

end
