function batch_get_stim_ln_of_motor(FolderName, FileName, iparams)
% batch_get_stim_ln_of_motor: Use Linear model to fit motor
%   variables using stimuli delivered as regressors
%
% Usage:
%   batch_get_stim_ln_of_motor(FolderName, FileName, iparams)
%
% Args:
%   FolderName: folders to use
%   FileName: files to use
%   iparams: parameters
%       (cDir: current directory)
%       (fo2reject: folders to reject)
%       (fi2reject: files to reject)
%       (metsuffix: suffix of files with wDat variable)
%           (default, '_prosmetadata')
%       (redo: redo even if sti_ln_mot exist, redo raw, redo shuffle)
%           (default, [0 0 0])
%       (varname: output variable name
%           (default, 'sti_ln_mot')
%       %%%%%%%%%%%% filter related %%%%%%%%%%%%
%       (filterlength_stim: filter length in seconds)
%           (default, 10)
%       (customtrial_init: number of seconds before stimuli 
%           start to count trial start (abs number))
%           (default, [])
%       (customtrial_end: number of seconds after stimuli 
%           end to count trial end (abs number))
%           (default, [])
%       (stim2use: which stimuli to use)
%           (default, [], or all)
%       (mot2use: motor variables to use)
%           (default, {'speed', 'fV', 'lV', 'yaw', 'pitch', 'roll'})
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
%       (pst_time: peri stimulus time to use to extract trial structure from whole trace (seconds))
%           (default, [-10 10])
%       %%%%%%%%%%%% significance related %%%%%%%%%%%%
%       (fdr: false discovery rate)
%            (default, 0.01)
%       (pval_cor_method: multiple comparison correction to use: dep, pdep, bh)
%            (default, 'dep')
%       %%%%%%%%%%%% motor variable related %%%%%%%%%%%%
%       (ball_radious: depth of directory search)
%           (default, 4.5)
%       (sig: std of gaussian kernel to smooth motor variable)
%           (default, 3)
%       (siz: size of kernel  to smooth motor variable)
%           (default, 10)
%       %%%%%%%%%%%% save summary plots %%%%%%%%%%%%
%       (oDir: output directory to save summary results)
%           (default, [pwd, filesep, 'sti_ln_mot'])
%
% Notes:
% for info about regression algorithm see:
%   runRidgeOnly (ridge regression using empirical Bayes)
% for info about surrogate data see https://en.wikipedia.org/wiki/Surrogate_data_testing
% for info about crossvalidation see https://en.wikipedia.org/wiki/Cross-validation_(statistics)

% Default params
pSMMot = []; 
pSMMot.cDir = pwd;
pSMMot.fo2reject = {'.', '..', 'preprocessed', ...
    'BData', 'rawtiff', 'motcor', 'stitch', ...
    'dfrel_vid', 'sti_ln_ca', 'sti_mot_ln_ca', ...
    'sti_ln_mot', 'roicov'};
pSMMot.fi2reject = {'Zstack'};
pSMMot.metsuffix = '_metadata';
pSMMot.redo = [0 0 0];
pSMMot.varname = 'sti_ln_mot';
pSMMot.filterlength_stim = [];
pSMMot.customtrial_init = [];
pSMMot.customtrial_end = [];
pSMMot.stim2use = [];
pSMMot.mot2use = {'speed', 'fV', 'lV', 'yaw', 'pitch', 'roll'};
pSMMot.trials2use = [];
pSMMot.minsize = 5;
pSMMot.per2use = 0.8;
pSMMot.type2run = [0 0 1];
pSMMot.btn = 10^4;
pSMMot.pst_time = [-10 10];
pSMMot.fdr = 0.01;
pSMMot.pval_cor_method = 'dep';
pSMMot.serverid = 'int';
pSMMot.corenum = 4;
pSMMot.chunksiz = 80;
pSMMot.febgate = 0;
pSMMot.ball_radious = 4.5;
pSMMot.sig = 3;
pSMMot.siz = 10;

pSMM_plot = [];
pSMM_plot.hbins = -1:0.01:1;
pSMM_plot.hbinsp = 0:0.01:1.1;
pSMM_plot.prct2use = 30;
pSMM_plot.fdr = 0.01;
pSMM_plot.mccor_method = 'dep';
pSMM_plot.dir_depth = 1;

% update variables
if ~exist('FileName', 'var'); FileName = []; end
if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('iparams', 'var'); iparams = []; end
pSMMot = loparam_updater(pSMMot, iparams);

pSMM_plot.fmetsuffix = [pSMMot.metsuffix, '.mat'];

% start pararell pool if not ready yet
ppobj = setup_parpool(pSMMot.serverid, pSMMot.corenum);

% Selecting folders
f2run = dir;
f2run = str2match(FolderName, f2run);
f2run = str2rm(pSMMot.fo2reject, f2run);
f2run = {f2run.name};
fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']); 
    
    cd(f2run{i});
    
    oDir = [pwd, filesep, 'sti_ln_mot'];
    
    runperfolder(FileName, pSMMot);
    
    cd(pSMMot.cDir);
    
    batch_plot_stim_ln_of_motor_results(FolderName, FileName, oDir, pSMM_plot)
    
    fprintf('\n')
    
end

delete_parpool(ppobj);

fprintf('... Done\n')

end

function runperfolder(FileName, pSMMot)
% runperfolder: run all files per folder
%
% Usage:
%   runperfolder(fname)
%
% Args:
%   fname: file name pattern
%   pSMMot: internal variable

% Files to load
f2run = rdir(['.', filesep, '*', pSMMot.metsuffix, '*.mat']);
f2run = str2match(FileName, f2run); 
f2run = {f2run.name};

fprintf('Estimating stimuli modulation for all ROIs\n')
fprintf(['Running n-files : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running file : ', strrep(f2run{i}, filesep, ' '), '\n']);
    gen_ln_model(f2run{i}, pSMMot);
        
end

fprintf('****** Done ******\n')

end

function gen_ln_model(filename, pSMMot)
% gen_ln_model: run stimuli modulation for each file
%
% Usage:
%   gen_ln_model(f2run, pSMMot)
%
% Args:
%   f2run: file name
%   pSMMot: internal variable

t0 = stic;

% 1) load required varaibles 'roi', and 'pSLM.field2load'
%   compatible with old format
load(filename, 'wDat')

% 2) get dt in seconds
dt = wDat.fTime(2) - wDat.fTime(1);

% 3) load motor variables to fit (fictrac and SVD of video)
motor_fictrac = load_edit_fictrac_motor_var(wDat, pSMMot.ball_radious, ...
    pSMMot.mot2use, pSMMot.sig, pSMMot.siz, 1);

[motor_SVD, ~, ~, ~] = ...
    load_edit_video_SVD(wDat, strrep(filename, '_metadata.mat', '_proc.mat'), ...
    pSMMot.sig, pSMMot.siz, 1);

motor_all = [motor_fictrac; motor_SVD];
clear motor_fictrac motor_SVD

[K, ~] = size(motor_all);

fprintf(['motor var to process: ', num2str(K), '\n'])

% 4) Get number of stimuli presented
stim_ntype = 1;

if ~contains(wDat.datatype, 'opto') || contains(wDat.datatype, 'prv')
    
    % get stimuli name
    stim_name = cellfun(@(x, y) [x '_AMP_' num2str(y)], ...
        wDat.sPars.name(1, :), chunk2cell(wDat.sPars.int(1, :), 1), ...
        'UniformOutput', false);
    [stim_name, ~, stim_u_idx] = unique(stim_name, 'stable');
     
    % get number of reps per stim
    stim_all_idx = stim_u_idx(wDat.sPars.order(1, :));
    stim_all_idx = stim_all_idx(1:size(wDat.sTime, 1), 1);
    stim_count = hist(stim_all_idx, 1:1:max(stim_all_idx(:)));
    stim_ntype = numel(unique(stim_name));
    
    % get stimulus post baseline time to define filter length
    for i = 1:max(stim_u_idx)
        stim_post_lag(i, 1) = min(wDat.sPars.basPost(stim_u_idx == i));
    end
    clear stim_u_idx i
       
end

% 5) use predefined trials
if isfield(wDat, 'trialsToUse') && ~isempty(wDat.trialsToUse) ...
    || isfield(wDat, 'trials2use') && ~isempty(wDat.trials2use)

    try pSMMot.trials2use = wDat.trialsToUse; end
    try pSMMot.trials2use = wDat.trials2use; end
    
    % update number of trials
    if ~isempty(pSMMot.trials2use)
        if sum(pSMMot.trials2use > min(stim_count)) == 0
            stim_count(:) = numel(pSMMot.trials2use);
        else
            fprintf('Error trials requested are outside range\n')
            return
        end
    end
    
end

% 6) generate filter length for stim (by default pre stimuli is 0)
if isempty(pSMMot.filterlength_stim)
    
    % Automatically calculate the max lag based on 
    %   the length of the stimuli divided by the sampling rate
    filterlength_tp(:, 2) = round(round(stim_post_lag)/dt);
    filterlength_tp(:, 1) = 0;
    
else
    
    filterlength_tp = round(pSMMot.filterlength_stim/dt);
    filterlength_tp = repmat(filterlength_tp, [1, stim_ntype]);
    
end

% 7) Define stim to use for model fitting
if isempty(pSMMot.stim2use)
    pSMMot.stim2use = 1:stim_ntype;
end

% 8) Build matrix of stimuli used for model with all lags
%   get a binary vector with stimuli presentation
sti_ln_mot.stim = getStimVect(wDat); 

% generate binary vectors for each stimulus to use pSLM.stim2use
if stim_ntype ~= 1
    jj = 1;
    for j = pSMMot.stim2use
        stim_m(jj, :) = getStimVect(wDat, 1, pSMMot.trials2use, j);
        jj = jj + 1;
    end
else
    stim_m(1, :) = sti_ln_mot.stim;
end

% build stimuli matrix with different lags (to match filter length)
[sti_ln_mot.stimM, ~] = build_matrix_with_lags(stim_m, filterlength_tp, 0);
% compile sti_ln_mot.stimM: t x tau x type
clear stim_m jj i j

display(['Stimuli to use: ', num2str(pSMMot.stim2use)])
display(['Variable name to use: ', num2str(pSMMot.varname)])

% 9) discretize timepoints into trials:
stim_trial_idx = getTrialVect(wDat, pSMMot.stim2use, ...
    pSMMot.customtrial_init, pSMMot.customtrial_end, pSMMot.trials2use);

% get the minimum number of trials of pSM.stim2use
%   to determine number of trials for training data
if exist('stim_count', 'var')
    
    min_ntrial = min(stim_count(pSMMot.stim2use));
    ntrial_train = floor(min_ntrial*pSMMot.per2use);
    
    % update sti_ln_mot.tidx2use
    idx2change = unique(stim_trial_idx, 'stable');
    idx2change(idx2change == 0) = [];
    stim_trial_idx_b = changem(stim_trial_idx, ...
        1:numel(idx2change), idx2change);
    
    sti_ln_mot.tidx2use = find(stim_trial_idx ~= 0 & ...
        stim_trial_idx_b <= min_ntrial);
    clear min_ntrial stim_count
    clear stim_trial_idx_b idx2change
    
else
    
    sti_ln_mot.tidx2use = find(stim_trial_idx ~= 0);
    ntrial_train = floor((size(wDat.sTime, 1)/stim_ntype)*pSMMot.per2use);
    
end

% prune un-used stims
stim_trial_idx = stim_trial_idx(sti_ln_mot.tidx2use);
sti_ln_mot.stimM = sti_ln_mot.stimM(sti_ln_mot.tidx2use, :);

% 10) define trials for training and testing
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
    iter_i t_i trialidx_u

% 11) Automatically calculate the minsize and 
%   number of splits base on the trial length
%   transform minsize to timepoints
chunk_minsize = pSMMot.minsize/dt;
chunk_splitn = floor(numel(sti_ln_mot.tidx2use)/chunk_minsize);

clear stim_trial_idx

% 12) Collect variables to use
motvar_in = motor_all(:, sti_ln_mot.tidx2use);
stocf(t0, 'Time consumed so far: ')

% 13) using LN model to fit motor variables

% Run Raw data
fprintf('Running Raw Data\n')
vList = whos('-file', filename);
sgate = sum(strcmp({vList.name}, pSMMot.varname));
if sgate && ~pSMMot.redo(1)
    fprintf([pSMMot.varname ' already exist in mat file\n']); 
    return;
end

% name of temporary file
filename_temp = [strrep(filename, [pSMMot.metsuffix, '.mat'], ''), ...
    '_temp.mat'];
tgate = exist(filename_temp, 'file');

if ~tgate || pSMMot.redo(2)
    
    % calculate filters and predicted traces for motor
    [~, lFilter, motvar_pred] = ridgeregres_per_row_crossval(...
        train_idx, motvar_in, sti_ln_mot.stimM, pSMMot.chunksiz, pSMMot.corenum, 'bayes');
    
    eVar = diag(corr(zscorebigmem(motvar_in)', zscorebigmem(motvar_pred)'));
    
    save(filename_temp, 'eVar', 'lFilter', 'motvar_pred', '-v7.3')
    
else
    
    load(filename_temp, 'eVar', 'lFilter', 'motvar_pred')
    
end

stocf(t0, 'Time consumed so far: ')

% 14) testing model significance
% Run (cs) Random chunk Shuffle
fprintf('Running Random shuffle in chunks\n')
lgate = 1;

if tgate
    vList = whos('-file', filename_temp);
    lgate = ~sum(contains({vList.name}, 'eVar_cs'));
end

stim_width = wDat.sTime(:, 2) - wDat.sTime(:, 1);
stim_width = max(stim_width*1.2);
add_stim_lag = [-ceil(6/dt) ceil(stim_width/dt)];
shuffle2use = 2;
fprintf(['adding stimuli lag of: ', ...
    num2str(add_stim_lag), ' (timestamps) \n'])

[~, eVar_cs] = get_explained_variance_shuffle(...
    motvar_in, motvar_pred, lgate, ...
    sti_ln_mot.stim(1, sti_ln_mot.tidx2use), ...
    pSMMot.redo(3), filename_temp, ...
    chunk_minsize, chunk_splitn, ...
    pSMMot.chunksiz, pSMMot.corenum, ...
    pSMMot.btn, shuffle2use, add_stim_lag, []);

% plot traces
% i = 1;
% plot(zscorebigmem(motvar_in(i, :)), 'r')
% hold on
% plot(motvar_pred(i, :), 'b')

save(filename_temp, 'eVar_cs', '-append')
    
stocf(t0, 'Time consumed so far: ')

% get p-values
sti_ln_mot.pval = calculate_pval_2(eVar, ...
    eVar_cs, pSMMot.fdr, pSMMot.pval_cor_method);

% split motor variable traces into trials
fprintf('Split trace into trials\n')
wDat_edit = wDat;
wDat_edit.fTime = wDat.fTime(sti_ln_mot.tidx2use);
[sti_ln_mot.motvar_per_trial, ~, ~, ~, ~, ~, ~, ~, ...
    sti_ln_mot.oStimidx, sti_ln_mot.stimvect_] = ...
    trace2trials(wDat_edit, zscorebigmem(motvar_in), ...
    pSMMot.pst_time, stim_all_idx, 1, []);
stocf(t0, 'Time consumed so far: ')
sti_ln_mot.stim_name = stim_name;

% generate baseline norm trace
matrix2norm = sti_ln_mot.motvar_per_trial;
for i = 1:numel(sti_ln_mot.stim_name)
    time_bs = [find(sti_ln_mot.oStimidx == i, 1, 'first'), ...
        find(sti_ln_mot.stimvect_ == i, 1,'first')];
    time_bs = time_bs(1):(time_bs(2)-1);
    time2edit = sti_ln_mot.oStimidx == i;
    out(:, i) = cf(@(x) (x(:, time2edit) - nanmean(x(:, time_bs), 2))./std(x(:, time_bs), [], 2), ...
            matrix2norm);
end

for i = 1:size(out, 1)
    sti_ln_mot.motvar_per_trial_norm{i, 1} = cell2mat(out(i, :));
end
clear matrix2norm

% compute correlation of the mean across half trials
trial_n = size(sti_ln_mot.motvar_per_trial{1}, 1);
group_perm_1 = nchoosek(1:trial_n, round(trial_n/2));
group_perm_1 = mat2cell(group_perm_1, ones(size(group_perm_1, 1), 1), size(group_perm_1, 2));
group_perm_2 = cf(@(x) setdiff(1:trial_n, x), group_perm_1);

% compute correlation across means
for i = 1:numel(group_perm_1)
    
    mean_g_1 = cf(@(x) nanmean(x(group_perm_1{i}, :), 1), ...
        sti_ln_mot.motvar_per_trial);
    mean_g_2 = cf(@(x) nanmean(x(group_perm_2{i}, :), 1), ...
        sti_ln_mot.motvar_per_trial);
    sti_ln_mot.mean_to_mean_corr(:, i) = cell2mat(cf(@(x, y) corr(x', y'), ...
        mean_g_1, mean_g_2));

end
clear mean_g_1 mean_g_2

% compute correlation across means
for i = 1:numel(group_perm_1)
    
    mean_g_1 = cf(@(x) nanmean(x(group_perm_1{i}, :), 1), ...
        sti_ln_mot.motvar_per_trial_norm);
    mean_g_2 = cf(@(x) nanmean(x(group_perm_2{i}, :), 1), ...
        sti_ln_mot.motvar_per_trial_norm);
    sti_ln_mot.mean_to_mean_norm_corr(:, i) = cell2mat(cf(@(x, y) corr(x', y'), ...
        mean_g_1, mean_g_2));

end
clear mean_g_1 mean_g_2

% delete temp file:
delete(filename_temp)

% 18) Save fields

% always updates this fields
sti_ln_mot.maxlag = filterlength_tp; 
sti_ln_mot.btn = pSMMot.btn; 
sti_ln_mot.stim2use = pSMMot.stim2use;
sti_ln_mot.eVar = eVar; 
sti_ln_mot.lFilter = lFilter;
sti_ln_mot.eVar_cs = eVar_cs;
sti_ln_mot.Y = motvar_in;
sti_ln_mot.train_idx = train_idx;

% save variable with the custom name
eval([pSMMot.varname, ' = sti_ln_mot']);
save(filename, pSMMot.varname, '-append')

% to generate predicted using sti_ln_mot
% Y_pred = get_predicted_signal(sti_ln_mot.train_idx, ...
%   sti_ln_mot.lFilter, sti_ln_mot.stimM);

end
