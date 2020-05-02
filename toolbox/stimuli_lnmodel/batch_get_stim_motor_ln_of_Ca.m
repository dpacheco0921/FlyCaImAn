function batch_get_stim_motor_ln_of_Ca(FolderName, FileName, iparams)
% batch_get_stim_motor_ln_of_Ca: make a linear model combining stimuli plus
% locomotor variables to explain ROI signal.
%
% Usage:
%   batch_get_stim_motor_ln_of_Ca(FolderName, FileName, iparams)
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
%       (redo: redo even if sti_mot_ln_ca exist, redo raw, redo shuffle)
%           (default, [0 0 0])
%       (varname: output variable name
%           (default, 'sti_mot_ln_ca')
%       %%%%%%%%%%%% filter related %%%%%%%%%%%%
%       (filterlength_stim: filter length in seconds before and after stimulus start)
%           (default, [])
%       (filterlength_motor: filter length in seconds before and after motor event)
%           (default, [-2 2])
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
%       %%%%%%%%%%%% motor variable related %%%%%%%%%%%%
%       (ball_radious: depth of directory search)
%           (default, 0)
%       (sig: std of gaussian kernel to smooth motor variable)
%       (siz: size of kernel  to smooth motor variable)
%       %%%%%%%%%%%% save summary plots %%%%%%%%%%%%
%       (oDir: output directory to save summary results)
%           (default, [pwd, filesep, 'sti_mot_ln_ca'])
%
% Notes:
% for info about regression algorithms see:
%   runRidgeOnly (ridge regression using empirical Bayes)
%   ridgeMML (ridge regression, regularization: marginal maximun likelihood)
% for info about surrogate data see https://en.wikipedia.org/wiki/Surrogate_data_testing
% for info about crossvalidation see https://en.wikipedia.org/wiki/Cross-validation_(statistics)

% Default params
pSLM = []; 
pSLM.cDir = pwd;
pSLM.fo2reject = {'.', '..', 'preprocessed', ...
    'BData', 'rawtiff', 'motcor', 'stitch', ...
    'dfrel_vid', 'sti_ln_ca', 'sti_mot_ln_ca', ...
    'sti_ln_mot', 'roicov'}; 
pSLM.fi2reject = {'Zstack'};
pSLM.fsuffix = '_prosroi'; 
pSLM.metsuffix = '_metadata';
pSLM.redo = [0 0 0];
pSLM.varname = 'sti_mot_ln_ca';
pSLM.filterlength_stim = [];
pSLM.filterlength_motor = [-2 2];
pSLM.customtrial_init = [];
pSLM.customtrial_end = [];
pSLM.stim2use = [];
pSLM.mot2use = {'speed', 'fV', 'lV', 'yaw', 'pitch', 'roll'};
pSLM.trials2use = [];
pSLM.minsize = 5;
pSLM.per2use = 0.8;
pSLM.type2run = [0 0 1];
pSLM.btn = 10^4;
pSLM.pst_time = [-10 10];
pSLM.fdr = 0.01;
pSLM.pval_cor_method = 'dep';
pSLM.serverid = 'int';
pSLM.corenum = 4;
pSLM.chunksiz = 80;
pSLM.febgate = 0;
pSLM.ball_radious = 9;
pSLM.sig = 3;
pSLM.siz = 10;

% update variables
if ~exist('FileName', 'var'); FileName = []; end
if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('iparams', 'var'); iparams = []; end
pSLM = loparam_updater(pSLM, iparams);

% start pararell pool if not ready yet
ppobj = setup_parpool(pSLM.serverid, pSLM.corenum);

% Selecting folders
f2run = dir;
f2run = str2match(FolderName, f2run);
f2run = str2rm(pSLM.fo2reject, f2run);
f2run = {f2run.name};
fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']); 
    
    cd(f2run{i});
    
    pSLM.oDir = [pwd, filesep, 'sti_mot_ln_ca'];
    if ~exist(pSLM.oDir, 'dir')
        mkdir(pSLM.oDir)
    end
    
    runperfolder(FileName, pSLM);
    cd(pSLM.cDir);
        
    fprintf('\n')
    
end

delete_parpool(ppobj);

fprintf('... Done\n')

end

function runperfolder(filename, pSLM)
% runperfolder: run all files per folder
%
% Usage:
%   runperfolder(filename)
%
% Args:
%   filename: file name pattern
%   pSLM: internal variable

% Files to load
f2run = rdir(['.', filesep, '*', pSLM.fsuffix, '*.mat']);
f2run = str2match(filename, f2run); 
f2run = {f2run.name};

fprintf('Estimating stimuli modulation for all ROIs\n')
fprintf(['Running n-files : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running file : ', strrep(f2run{i}, filesep, ' '), '\n']);
    gen_ln_model(f2run{i}, pSLM);
    
end

fprintf('****** Done ******\n')

end

function gen_ln_model(filename, pSLM)
% gen_ln_model: run stimuli modulation for each file
%
% Usage:
%   gen_ln_model(filename, pSLM)
%
% Args:
%   filename: file name
%   pSLM: internal variable

t0 = stic;

% 1) load required varaibles 'roi', and 'wDat'
%   compatible with old format
load(filename, 'roi');
load(strrep(filename, pSLM.fsuffix, pSLM.metsuffix), 'wDat')

% 2) get dt in seconds
dt = wDat.fTime(2) - wDat.fTime(1);

% get number or ROIs
[K, ~] = size(roi.filtered.dfof);

fprintf(['ROIs to process: ', num2str(K), '\n'])

% 3) get number of stimuli presented
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
    clear stim_u_idx stim_name i
       
end

% 4) use predefined trials
if isfield(wDat, 'trialsToUse') && ~isempty(wDat.trialsToUse) ...
    || isfield(wDat, 'trials2use') && ~isempty(wDat.trials2use)

    try pSLM.trials2use = wDat.trialsToUse; end
    try pSLM.trials2use = wDat.trials2use; end
    
    % update number of trials
    if ~isempty(pSLM.trials2use)
        if sum(pSLM.trials2use > min(stim_count)) == 0
            stim_count(:) = numel(pSLM.trials2use);
        else
            fprintf('Error trials requested are outside range\n')
            return
        end
    end
    
end

% 5) generate filter length for stim (by deafult pre stimuli is 0)
if isempty(pSLM.filterlength_stim)
    
    % Automatically calculate the max lag based on 
    %   the length of the stimuli divided by the sampling rate
    filterlength_tp(:, 2) = round(round(stim_post_lag)/dt);
    filterlength_tp(:, 1) = 0;
    
else
    
    filterlength_tp = round(pSLM.filterlength_stim/dt);
    filterlength_tp = repmat(filterlength_tp, [1, stim_ntype]);
    
end

% 6) Define stim to use for model fitting
if isempty(pSLM.stim2use)
    pSLM.stim2use = 1:stim_ntype;
end

% 7) Build matrix of stimuli used for model with all lags
%   get a binary vector with stimuli presentation
sti_mot_ln_ca.stim = getStimVect(wDat); 

% generate binary vectors for each stimulus to use pSLM.stim2use
if stim_ntype ~= 1
    jj = 1;
    for j = pSLM.stim2use
        stim_m(jj, :) = getStimVect(wDat, 1, pSLM.trials2use, j);
        jj = jj + 1;
    end
else
    stim_m(1, :) = sti_mot_ln_ca.stim;
end

% build stimuli matrix with different lags (to match filter length)
[sti_mot_ln_ca.stimM, idx_stimM] = build_matrix_with_lags(stim_m, filterlength_tp, 0);
% compile sti_mot_ln_ca.stimM: t x tau x type
sti_mot_ln_ca.varnames = chunk2cell(pSLM.stim2use, 1);
sti_mot_ln_ca.varnames = cf(@(x) ['stim_', num2str(x)], sti_mot_ln_ca.varnames)';
stim_type = ones(length(idx_stimM), 1);
clear stim_m jj

% 8) load motor variables directly from fictrac (fV, lV, V, yaw-vel)
motor_m = load_edit_fictrac_motor_var(wDat, pSLM.ball_radious, ...
    pSLM.mot2use, pSLM.sig, pSLM.siz, 1);
filterlength_tp_ = repmat(pSLM.filterlength_motor/dt, [size(motor_m, 1), 1]);
[sti_mot_ln_ca.stimM(:, end+1:end+sum(sum(abs(filterlength_tp_), 2))), idx_motorM] = ...
    build_matrix_with_lags(motor_m, filterlength_tp_, nan);
sti_mot_ln_ca.idx_stimM = [idx_stimM; idx_motorM + max(idx_stimM(:))];
sti_mot_ln_ca.stim_type = [stim_type; ones(length(idx_motorM), 1)*2];
sti_mot_ln_ca.varnames = [sti_mot_ln_ca.varnames; pSLM.mot2use'];
clear motorM idx_motorM idx_stimM stim_type

% 9) load additional motor variables SVD of cropped video
motor_SVD = load_edit_video_SVD(wDat, strrep(filename, pSLM.fsuffix, '_proc'), pSLM.sig, pSLM.siz, 1);
filterlength_tp_ = repmat([0 1], [size(motor_SVD, 1), 1]);
[sti_mot_ln_ca.stimM(:, end+1:end+sum(sum(abs(filterlength_tp_), 2))), idx_motorSVD] = ...
    build_matrix_with_lags(motor_SVD, filterlength_tp_, nan);
sti_mot_ln_ca.idx_stimM = [sti_mot_ln_ca.idx_stimM; idx_motorSVD + max(sti_mot_ln_ca.idx_stimM(:))];
sti_mot_ln_ca.stim_type = [sti_mot_ln_ca.stim_type; ones(length(idx_motorSVD), 1)*3];

motorSVD_varnames = chunk2cell(1:size(motor_SVD, 1), 1);
motorSVD_varnames = cf(@(x) ['SVD_', num2str(x)], motorSVD_varnames)';
sti_mot_ln_ca.varnames = [sti_mot_ln_ca.varnames; motorSVD_varnames];
clear motorM idx_motorSVD motorSVD_varnames

display(['Stimuli to use: ', num2str(pSLM.stim2use)])
display(['Motor variables to use: ', pSLM.mot2use])
display(['Variable name to use: ', num2str(pSLM.varname)])

% 10) discretize timepoints into trials:
stim_trial_idx = getTrialVect(wDat, pSLM.stim2use, ...
    pSLM.customtrial_init, pSLM.customtrial_end, pSLM.trials2use);

% get the minimum number of trials of pSLM.stim2use
%   to determine number of trials for training data
if exist('stim_count', 'var')
    
    min_ntrial = min(stim_count(pSLM.stim2use));
    ntrial_train = floor(min_ntrial*pSLM.per2use);
    
    % update sti_mot_ln_ca.tidx2use
    idx2change = unique(stim_trial_idx, 'stable');
    idx2change(idx2change == 0) = [];
    stim_trial_idx_b = changem(stim_trial_idx, ...
        1:numel(idx2change), idx2change);
    
    sti_mot_ln_ca.tidx2use = find(stim_trial_idx ~= 0 & ...
        stim_trial_idx_b <= min_ntrial);
    clear min_ntrial stim_count
    clear stim_trial_idx_b idx2change
    
else
    
    sti_mot_ln_ca.tidx2use = find(stim_trial_idx ~= 0);
    ntrial_train = floor((size(wDat.sTime, 1)/stim_ntype)*pSLM.per2use);
    
end

% remove nan timepoints
sti_mot_ln_ca.tidx2use = setdiff(sti_mot_ln_ca.tidx2use, ...
    find(isnan(sum(sti_mot_ln_ca.stimM, 2))));

% prune un-used stims
stim_trial_idx = stim_trial_idx(sti_mot_ln_ca.tidx2use);
sti_mot_ln_ca.stimM = sti_mot_ln_ca.stimM(sti_mot_ln_ca.tidx2use, :);

% 11) Automatically calculate the minsize and 
%   number of splits base on the trial length
%   transform minsize to timepoints
chunk_minsize = pSLM.minsize/dt;
chunk_splitn = floor(numel(sti_mot_ln_ca.tidx2use)/chunk_minsize);

clear stim_trial_idx

% 12) test if regressors are orthogonal
%   the resulting plot ranges from 0 to 1 for each regressor, with 1 being
%   fully orthogonal to all preceeding regressors in the matrix and 0 being
%   fully redundant.
rejIdx = false(1, size(sti_mot_ln_ca.stimM, 2));
[~, fullQRR] = qr(bsxfun(@rdivide, sti_mot_ln_ca.stimM, sqrt(sum(sti_mot_ln_ca.stimM.^2, 1))), 0);

% plot orthogonality of regressors
filename_ = strrep(strrep(filename, [pSLM.fsuffix, '.mat'], ''), ['.', filesep], '');
plot_orthogonal_test(fullQRR, filename_, pSLM.oDir)

% check if design matrix is full rank
if sum(abs(diag(fullQRR)) > max(size(sti_mot_ln_ca.stimM)) * ...
        eps(fullQRR(1))) < size(sti_mot_ln_ca.stimM, 2)
    temp = ~(abs(diag(fullQRR)) > max(size(sti_mot_ln_ca.stimM)) * eps(fullQRR(1)));
    fprintf(['Design matrix is rank-defficient. ', ...
        'Removing %d/%d additional regressors.\n'], sum(temp), sum(~rejIdx));
    % reject regressors that cause rank-defficient matrix
    rejIdx(~rejIdx) = temp;
    keyboard
end

% 13) Collect Ca trace and Stim
% prune unused stims
catrace_in = roi.filtered.dfof(:, sti_mot_ln_ca.tidx2use);

stim_in = sti_mot_ln_ca.stimM; 
stimSiz = size(stim_in, 2);
stim_bin = sti_mot_ln_ca.stim(1, sti_mot_ln_ca.tidx2use);

stocf(t0, 'Time consumed so far: ')

% 14) get data splitted in trials and get mean and median
fprintf('Split trace into trials\n')
wDat_edit = wDat;
wDat_edit.fTime = wDat.fTime(sti_mot_ln_ca.tidx2use);

[roi_per_trial, rel_time, ~, ~, ~, ~, ~, ~, trialvect_, stimvect_] = trace2trials(wDat_edit, ...
    zscorebigmem(catrace_in), pSLM.pst_time, stim_all_idx, 1, []);
[motor_per_trial, ~, ~, ~, ~, ~, ~, ~, ~, ~] = trace2trials(wDat_edit, ...
    stim_in', pSLM.pst_time, stim_all_idx, 1, []);
stocf(t0, 'Time consumed so far: ')

% zscore signal (for selected trials) per ROI
fprintf('Get variance (total, stimuli-related, motor-related)\n')
roi_per_trial = cf(@(x, y) reshape(zscorebigmem(x(:)'), y), roi_per_trial, ...
    cf(@(x) size(x), roi_per_trial));

% get total variance
var_total = cell2mat(cf(@(x) var(x(:), 1), roi_per_trial));

% get stim variance
var_stim_mean = cell2mat(cf(@(x) var(x(:), 1), cf(@(x) mean(x, 1), roi_per_trial)));
%evar_stim_mean = cell2mat(cf(@(x, y) corr(horz(x')', horz(y')').^2, ...
%    roi_per_trial, cf(@(x) repmat(mean(x, 1), [size(x, 1), 1]), roi_per_trial)));

% get residual (total-stim)
var_res_mean = cell2mat(cf(@(x) var(x(:), 1), cf(@(x) x - mean(x, 1), roi_per_trial)));

% fit motor to residual
roi_per_trial_res_mean = cf(@(x) horz(x'), cf(@(x) x - mean(x, 1), roi_per_trial));
roi_per_trial_res_mean = cell2mat(roi_per_trial_res_mean);

motor_per_trial_flat = cf(@(x) horz(x'), motor_per_trial);
motor_per_trial_flat = cell2mat(motor_per_trial_flat);

% re-define trials for training and testing
stim_trial_idx_2 = ones(size(roi_per_trial{1}));
for i = 1:size(stim_trial_idx_2, 1)
    stim_trial_idx_2(i, :) = stim_trial_idx_2(i, :)*i;
end
stim_trial_idx_2 = horz(stim_trial_idx_2');

trialidx_u = unique(stim_trial_idx_2);
test_comb = chunk2cell(randperm (numel(trialidx_u)), numel(trialidx_u) - ntrial_train)';
trial_comb = cf(@(x) setdiff(trialidx_u, x), test_comb);
train_idx = zeros(size(trial_comb, 1), size(stim_trial_idx_2, 2), 'logical');

for iter_i = 1:size(trial_comb)
    for t_i = 1:length(trial_comb{iter_i})
       train_idx(iter_i, stim_trial_idx_2 == trial_comb{iter_i}(t_i)) = 1;
    end
end

clear trial_comb test_comb ...
    ntrial_train iter_i t_i trialidx_u

% Run Raw data
fprintf('Running Raw Data\n')
vList = whos('-file', filename);
sgate = sum(strcmp({vList.name}, pSLM.varname));
if sgate && ~pSLM.redo(1)
    fprintf([pSLM.varname ' already exist in mat file\n']); 
    return;
end

% name of temporary file
filename_temp = [strrep(filename, [pSLM.fsuffix, '.mat'], ''), '_temp.mat'];
tgate = exist(filename_temp, 'file');

if ~tgate || pSLM.redo(2)
    
    % 15) using LN model to fit residual
    [~, LN_filter_motor, motor_prediction, motor_prediction_single, ~] = ...
        ridgeregres_per_row_crossval(...
        train_idx, roi_per_trial_res_mean, ...
        motor_per_trial_flat(sti_mot_ln_ca.stim_type ~= 1, :)', ...
        pSLM.chunksiz, pSLM.corenum, 'bayes', sti_mot_ln_ca.stim_type(sti_mot_ln_ca.stim_type ~= 1));
    
    save(filename_temp, 'LN_filter_motor', ...
        'motor_prediction', 'motor_prediction_single', '-v7.3')
    
else
    
    load(filename_temp, 'LN_filter_motor', ...
        'motor_prediction', 'motor_prediction_single')
    
end

% motor prediction per trial
% [motor_prediction_per_trial, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
%     trace2trials(wDat_edit, ...
%     motor_prediction, [-10 10], stim_all_idx, 1, []);

% get motor variance
var_motor = var(motor_prediction', 1, 1)';
evar_motor = corr(roi_per_trial_res_mean', motor_prediction');
evar_motor = diag(evar_motor).^2;

for i = 1:size(motor_prediction_single, 1)
    var_motor_single(:, i) = var(squeeze(motor_prediction_single(i, :, :)), 1, 1)';
    temp_ = corr(roi_per_trial_res_mean', squeeze(motor_prediction_single(i, :, :)));
    evar_motor_single(:, i) = diag(temp_).^2;
end

% plot filters
plot_filters_full_model(LN_filter_motor, ...
    sti_mot_ln_ca, evar_motor, var_stim_mean, var_total, filename_, pSLM.oDir)

% plot variance and explained variance in full vs 
%   single variable type model
plot_full_vs_single_vartype_model(...
    var_motor_single, var_motor, evar_motor_single, ...
    evar_motor, filename_, pSLM.oDir)

% plot stim var vs motor var
plot_evar_stim_vs_motor(var_stim_mean, evar_motor, var_total, ...
    motor_prediction, rel_time, stimvect_, roi_per_trial, ...
    roi_per_trial_res_mean, filename_, pSLM.oDir)

stocf(t0, 'Time consumed so far: ')

% 16) testing model significance
% Run (cs) Random chunk Shuffle
lgate = 1;
    
if tgate
    vList = whos('-file', filename_temp);
    lgate = ~sum(contains({vList.name}, 'eVar_cs'));
end
    
stim_bin_2 = repmat(stimvect_, [size(roi_per_trial{1}, 1) 1]);
stim_bin_2 = horz(stim_bin_2');
chunk_splitn_2 = floor(numel(stim_bin_2)/chunk_minsize);
sti_mot_ln_ca.stim_2 = stim_bin_2;
sti_mot_ln_ca.stim_trial = stimvect_;
sti_mot_ln_ca.rel_time_trial = rel_time;

add_stim_lag = [-10 10];
eVar_cs = get_explained_variance_shuffle(...
    roi_per_trial_res_mean, motor_prediction, lgate, ...
    stim_bin_2, pSLM.redo(3), filename_temp, ...
    chunk_minsize, chunk_splitn_2, ...
    pSLM.chunksiz, pSLM.corenum, pSLM.btn, 1, add_stim_lag);

stocf(t0, 'Time consumed so far: ')

% plot significance of motor explained variance
plot_sig_motor_mod(evar_motor, eVar_cs, filename_, pSLM.oDir)

% get p-values
sti_mot_ln_ca.pval = calculate_pval_2(evar_motor, ...
    eVar_cs, pSLM.fdr, pSLM.pval_cor_method);

% delete temp file:
%delete(filename_temp)

% 17) Save fields
sti_mot_ln_ca.filter_lags = [filterlength_tp; filterlength_tp_]; 
sti_mot_ln_ca.btn = pSLM.btn; 
sti_mot_ln_ca.stim2use = pSLM.stim2use;
sti_mot_ln_ca.var_total = var_total;
sti_mot_ln_ca.var_stim = var_stim_mean;
sti_mot_ln_ca.var_res = var_res_mean;
sti_mot_ln_ca.var_motor = var_motor;
sti_mot_ln_ca.evar_motor = evar_motor;
sti_mot_ln_ca.evar_motor_single = var_motor_single;
sti_mot_ln_ca.evar_motor_single = evar_motor_single;
sti_mot_ln_ca.evar_motor_cs = eVar_cs;
sti_mot_ln_ca.lFilter_motor = LN_filter_motor;

% save variable with the custom name
eval([pSLM.varname, ' = sti_mot_ln_ca']);
save(filename, pSLM.varname, '-append')

end

function plot_sig_motor_mod(eVar, eVar_shuffle, filename, oDir)
% plot_sig_motor_mod: plot significance of motor modulation
%
% Usage:
%   plot_sig_motor_mod(eVar, eVar_shuffle, filename, oDir)
%
% Args:
%   eVar: explained variance of raw data
%   eVar_shuffle: explained variance of shuffle data
%   filename: file name
%   oDir: output directory

motpar.mccor_method = 'dep';
motpar.fdr = 0.01;
motpar.hbins = 0:0.001:1;
motpar.hbinsp = 0:0.001:1.1;
motpar.hrange = [0 0.2];
motpar.irange = [0 .001];

% correct explained variance
[eVar_stat, eVar_raw, eVar_shuffle] = ...
    correct_corrcoef(eVar, eVar_shuffle, 1);

% get p-values
pval_raw = calculate_pval_2(eVar_raw, ...
    eVar_shuffle, motpar.fdr, 'raw');

pval_cor = calculate_pval_2(eVar_raw, ...
    eVar_shuffle, motpar.fdr, motpar.mccor_method);

selIdx = pval_cor <= motpar.fdr;

% get histogram of stim mod vs non-stim mod
eVar_hist = hist(flat_matrix(eVar_raw(selIdx, :)), ...
    motpar.hbins);
null_eVar = hist(flat_matrix(eVar_raw(~selIdx, :)), ...
    motpar.hbins);
eVar_hist = bsxfun(@rdivide, ...
    eVar_hist, sum(eVar_hist, 2));
null_eVar = bsxfun(@rdivide, ...
    null_eVar, sum(null_eVar, 2));

% generate histograms per ROI
roi_n = size(eVar_raw, 1);
for roi_i = 1:roi_n
    
    eVar_hist(roi_i, :) = ...
        hist(eVar_raw(roi_i, :), motpar.hbins);
    eVar_hist_null(roi_i, :) = ...
        hist(eVar_shuffle(roi_i, :), motpar.hbins);
    
end
fprintf('\n')

% generate histograms
roi_n = size(eVar_stat, 1);

% plot histogram
figH = figure('Position', ...
    genfigpos(1, 'center', [1600 700])); 
axH(1) = subplot(2, 4, 1);
axH(2) = subplot(2, 4, 5);

axH(3) = subplot(2, 4, 2);
axH(4) = subplot(2, 4, 6);

axH(5) = subplot(2, 4, 3);
axH(6) = subplot(2, 4, 7);

axH(7) = subplot(2, 4, 4);
axH(8) = subplot(2, 4, 8);

% collect all
eVar_hist_all = sum(eVar_hist, 1);
eVar_hist_null_all = sum(eVar_hist_null, 1);

% normalize per column
eVar_hist_all = bsxfun(@rdivide, ...
    eVar_hist_all, sum(eVar_hist_all, 2));
eVar_hist_null_all = bsxfun(@rdivide, ...
    eVar_hist_null_all, sum(eVar_hist_null_all, 2));

eVar_hist = bsxfun(@rdivide, ...
    eVar_hist, sum(eVar_hist, 2));
eVar_hist_null = bsxfun(@rdivide, ...
    eVar_hist_null, sum(eVar_hist_null, 2));

% plot histogram of corrcoef from raw and null
icolormap = [colorGradient([1 1 1], [1 0 1], 15); ...
    colorGradient([1 0 1], [0 1 1], 15)];

[~, idx_order] = sort(eVar_stat);

imagesc(motpar.hbins, 1:roi_n, ...
    eVar_hist_null(idx_order, :), ...
    'Parent', axH(1))
colormap(axH(1), icolormap)
caxis(axH(1), motpar.irange)

imagesc(motpar.hbins, 1:roi_n, ...
    eVar_hist(idx_order, :), ...
    'Parent', axH(2))
colormap(axH(2), icolormap)
caxis(axH(2), [0 1])

axH(1).YDir = 'normal';
axH(2).YDir = 'normal';
axH(1).Title.String = 'CC Shuffle data';
axH(2).Title.String = 'CC Raw data';
axH(1).YLabel.String = 'ROI number';
axH(2).YLabel.String = 'ROI number'; 
axH(2).XLabel.String = 'R^2';
axH(1).XLim = motpar.hrange;
axH(2).XLim = motpar.hrange;

cbH(1) = colorbar('peer', axH(1)); 
cbH(1).Label.String = 'Probability'; 
cbH(1).Label.FontSize = 10;

cbH(2) = colorbar('peer', axH(2)); 
cbH(2).Label.String = 'Probability'; 
cbH(2).Label.FontSize = 10;

% plot histogram modulated vs non-modulated
lineH(1) = plot(motpar.hbins, eVar_hist, ...
    'k', 'Linewidth', 2, 'Parent', axH(3));
hold(axH(3), 'on')
lineH(2) = plot(motpar.hbins, null_eVar, ...
    'Color', [0.5 0.5 0.5], ...
    'Linewidth', 2, 'Parent', axH(3));

axH(3).Title.String = 'mod vs ~mod ROIs';
axH(3).YLabel.String = 'Probability'; 
axH(3).XLim = motpar.hrange;

legend(axH(3), lineH, {'mod', '~mod'}, ...
    'location', 'northwest')

% plot histogram shuffle vs raw
lineH(1) = plot(motpar.hbins, eVar_hist_null_all, ...
    'k', 'Linewidth', 2, 'Parent', axH(4));
hold(axH(4), 'on')
lineH(2) = plot(motpar.hbins, eVar_hist_all, ...
    'Color', [0.5 0.5 0.5], ...
    'Linewidth', 2, 'Parent', axH(4));

axH(4).Title.String = 'shuffle vs raw';
axH(4).XLabel.String = 'R^2'; 
axH(4).YLabel.String = 'Probability'; 
axH(4).XLim = motpar.hrange;

legend(axH(4), lineH, {'shuffle', 'raw'}, ...
    'location', 'northwest')

% plot histogram shuffle vs raw
pval_raw_hist = hist(pval_raw, motpar.hbinsp);
pval_cor_hist = hist(pval_cor, motpar.hbinsp);
pval_raw_hist = bsxfun(@rdivide, ...
    pval_raw_hist, sum(pval_raw_hist, 2));
pval_cor_hist = bsxfun(@rdivide, ...
    pval_cor_hist, sum(pval_cor_hist, 2));

lineH(1) = plot(motpar.hbinsp, pval_raw_hist, ...
    'k', 'Linewidth', 2, 'Parent', axH(5));
hold(axH(5), 'on')
lineH(2) = plot(motpar.hbinsp, pval_cor_hist, ...
    'Color', [0.5 0.5 0.5], ...
    'Linewidth', 2, 'Parent', axH(5));

legend(axH(5), lineH, {'pval', 'corpval'}, ...
    'location', 'northeast')

plot(motpar.hbinsp, pval_raw_hist, ...
    'k', 'Linewidth', 2, 'Parent', axH(6));
hold(axH(6), 'on')
plot(motpar.hbinsp, pval_cor_hist, ...
    'Color', [0.5 0.5 0.5], ...
    'Linewidth', 2, 'Parent', axH(6));

axH(5).Title.String = 'pvalue';
axH(5).XLim = [0 0.9];
axH(5).XLabel.String = 'pvalue';
axH(5).YLabel.String = 'Probability';

axH(6).Title.String = 'pvalue zoom';
axH(6).XLim = [0 0.1];
axH(6).XLabel.String = 'pvalue';
axH(6).YLabel.String = 'Probability';

legend(axH(6), lineH, {'pval', 'corpval'}, ...
    'location', 'northeast')

plot(pval_raw, eVar_stat, ...
    '.k', 'MarkerSize', 20, 'Parent', axH(7));
hold(axH(7), 'on')
plot(pval_raw(pval_cor <= motpar.fdr), ...
    eVar_stat(pval_cor <= motpar.fdr), ...
    '.c', 'MarkerSize', 20, 'Parent', axH(7));

plot(pval_cor, eVar_stat, ...
    '.k', 'MarkerSize', 20, 'Parent', axH(8));

axH(7).Title.String = ['pvalue vs cc (', ...
    num2str(sum(pval_raw <= motpar.fdr)), ', ', ...
    num2str(sum(pval_raw <= motpar.fdr)*100/roi_n), ') (#, %)'];
axH(7).XLim = [0 0.05];
axH(7).XLabel.String = 'pvalue';
axH(7).YLabel.String = 'R^2';
line(.01*ones(1, 2), ylim(axH(7)), 'color', 'r', 'Parent', axH(7))

axH(8).Title.String = ['pvalue-cor vs cc (', ...
    num2str(sum(pval_cor <= motpar.fdr)), ', ', ...
    num2str(sum(pval_cor <= motpar.fdr)*100/roi_n), ') (#, %)'];
axH(8).XLim = [0 0.05];
axH(8).XLabel.String = 'pvalue';
axH(8).YLabel.String = 'R^2';
line(.01*ones(1, 2), ylim(axH(8)), 'color', 'r', 'Parent', axH(8))

fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;

figformat = [1 0 0 0 0 0 0 0 1 0 1];
save_edit_fig_int(axH, figH, oDir, ...
    [filename, '_sig_motor_mod'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH)

end

function plot_evar_stim_vs_motor(...
    var_stim_mean, evar_motor, var_total, motor_prediction, rel_time, ...
    stimvect_, roi_per_trial, roi_per_trial_res_mean, filenamefig, oDir)

% plot summary results
figH = figure('Position', genfigpos(1, 'center', [965 881]));

axH(1) = subplot(6, 3, 1);
axH(2) = subplot(6, 3, 2);
axH(3) = subplot(6, 3, 4);
axH(4) = subplot(6, 3, 5);
axH(5) = subplot(6, 3, 7);
axH(6) = subplot(6, 3, 8);
axH(7) = subplot(6, 3, 10);
axH(8) = subplot(6, 3, 11);
axH(9) = subplot(6, 3, [13 14]);
axH(10) = subplot(6, 3, [16 17]);
axH(11) = subplot(6, 3, [15 18]);
axH(12) = subplot(6, 3, 3);
axH(13) = subplot(6, 3, 6);
axH(14) = subplot(6, 3, 9);
axH(15) = subplot(6, 3, 12);

axH(11).Position = [0.6916 0.0700 0.2134 0.2700];
axH(10).Position = [0.1300 0.0700 0.4942 0.1026];

% find examples
[~, stim_idx] = sort(var_stim_mean);
[~, motor_idx] = sort(evar_motor);
[~, mixed_idx] = sort(var_stim_mean.*evar_motor, 'descend');
mixed_idx = mixed_idx(1);
yrange = [-2.5 4];

plot(rel_time, stimvect_, 'Parent', axH(1))
axH(1).XLim = rel_time([1 end]);
axH(1).Title.String = 'auditory stimuli (example auditory unit)';

plot(rel_time, stimvect_, 'Parent', axH(2))
axH(2).XLim = rel_time([1 end]);
axH(2).Title.String = 'auditory stimuli (example motor unit)';

plot(rel_time, stimvect_, 'Parent', axH(12))
axH(12).XLim = rel_time([1 end]);
axH(12).Title.String = 'auditory stimuli (example mixed unit)';

plot(rel_time, roi_per_trial{stim_idx(end)}', 'Parent', axH(3))
hold(axH(3), 'on');
plot(rel_time, mean(roi_per_trial{stim_idx(end)}, 1), 'k', 'Parent', axH(3))
axH(3).XLim = rel_time([1 end]);
axH(3).YLim = yrange;
axH(3).Title.String = [];

plot(rel_time, roi_per_trial{motor_idx(end)}', 'Parent', axH(4))
hold(axH(4), 'on');
plot(rel_time, mean(roi_per_trial{motor_idx(end)}, 1), 'k', 'Parent', axH(4))
axH(4).XLim = rel_time([1 end]);
axH(4).YLim = yrange;
axH(4).Title.String = 'Ca response across trials';

plot(rel_time, roi_per_trial{mixed_idx}', 'Parent', axH(13))
hold(axH(13), 'on');
plot(rel_time, mean(roi_per_trial{mixed_idx}, 1), 'k', 'Parent', axH(13))
axH(13).XLim = rel_time([1 end]);
axH(13).YLim = yrange;
axH(13).Title.String = [];

plot(rel_time, (roi_per_trial{stim_idx(end)} - mean(roi_per_trial{stim_idx(end)}, 1))', ...
    'Parent', axH(5))
axH(5).XLim = rel_time([1 end]);
axH(5).YLim = yrange;
axH(5).Title.String = [];

plot(rel_time, (roi_per_trial{motor_idx(end)} - mean(roi_per_trial{motor_idx(end)}, 1))', ...
    'Parent', axH(6))
axH(6).XLim = rel_time([1 end]);
axH(6).YLim = yrange;
axH(6).Title.String = 'residual response across trials (- avg trial response)';

plot(rel_time, (roi_per_trial{mixed_idx} - mean(roi_per_trial{mixed_idx}, 1))', ...
    'Parent', axH(14))
axH(14).XLim = rel_time([1 end]);
axH(14).YLim = yrange;
axH(14).Title.String = [];

timestamps = 1:size(roi_per_trial_res_mean, 2);
timestamps = (timestamps*diff(rel_time([1 2]))) + rel_time(1);
plot(timestamps, roi_per_trial_res_mean(stim_idx(end), :), 'k', 'Parent', axH(7))
hold(axH(7), 'on');
plot(timestamps, motor_prediction(stim_idx(end), :), 'r', 'Parent', axH(7))
axH(7).XLim = rel_time([1 end]);
axH(7).YLim = yrange;
axH(7).Title.String = [];

plot(timestamps, roi_per_trial_res_mean(motor_idx(end), :), 'k', 'Parent', axH(8))
hold(axH(8), 'on');
plot(timestamps, motor_prediction(motor_idx(end), :), 'r', 'Parent', axH(8))
axH(8).XLim = rel_time([1 end]);
axH(8).YLim = yrange;
axH(8).Title.String = ['residual and \color{red} ', ...
    'motor fit \color{black}(single trial)'];

plot(timestamps, roi_per_trial_res_mean(mixed_idx, :), 'k', 'Parent', axH(15))
hold(axH(15), 'on');
plot(timestamps, motor_prediction(mixed_idx, :), 'r', 'Parent', axH(15))
axH(15).XLim = rel_time([1 end]);
axH(15).YLim = yrange;
axH(15).Title.String = [];

plot(var_total(stim_idx), '.k', 'Parent', axH(9))
hold(axH(9), 'on');
plot(evar_motor(stim_idx) + var_stim_mean(stim_idx), '.r', 'Parent', axH(9))
plot(var_stim_mean(stim_idx), '.b', 'Parent', axH(9))
axH(9).XLim = [1 length(var_total)];
axH(9).Title.String = ['variance \color{black}(total, ', ...
    '\color{blue}stimuli, \color{red}motor)'];
axH(9).XLabel.String = 'ROIs';
axH(9).YLabel.String = 'variance';

var_res = var_total - var_stim_mean;
plot(evar_motor(stim_idx), 'r', 'Parent', axH(10))
hold(axH(10), 'on');
plot(var_res(stim_idx), 'm', 'Parent', axH(10))
axH(10).XLim = [1 length(var_total)];
axH(10).Title.String = ['variance not explained by stimuli ', ...
    '(\color{red}motor, \color{magenta}unexplained)'];
axH(10).XLabel.String = 'ROIs';
axH(10).YLabel.String = 'variance';

plot(var_stim_mean, evar_motor, '.k', 'Parent', axH(11))
hold(axH(11), 'on');
plot(var_stim_mean(stim_idx(end)), evar_motor(stim_idx(end)), '*g', 'Parent', axH(11))
plot(var_stim_mean(motor_idx(end)), evar_motor(motor_idx(end)), '+g', 'Parent', axH(11))
plot(var_stim_mean(mixed_idx), evar_motor(mixed_idx), '+g', 'Parent', axH(11))
axH(11).XLim = [0 1];
axH(11).YLim = [0 1];
axH(11).Title.String = 'stimuli vs motor variance';
axH(11).XLabel.String = 'stimuli variance';
axH(11).YLabel.String = 'motor variance';

fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;

figformat = [1 0 0 0 0 0 0 0 1 0 1];
save_edit_fig_int(axH, figH, oDir, ...
    [filenamefig, '_stim_vs_motor_var'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH)

end

function plot_filters_full_model(LN_filter_motor, ...
    sti_mot_ln_ca, evar_motor, var_stim_mean, var_total, filenamefig, oDir)
% plot_filters_full_model: plots variance and explained variance
%   of full vs single var type model
%
% Usage:
%   plot_filters_full_model(LN_filter_motor, ...
%       sti_mot_ln_ca, evar_motor, var_stim_mean, filenamefig, oDir)
%
% Args:
%   LN_filter_motor: filters from LN model
%   sti_mot_ln_ca: internal variable
%   evar_motor: explained variance of full model [x, 1]
%   var_stim_mean: variance of stimulus [x, 1]
%   var_total: total variance [x, 1]
%   filenamefig: filename
%   oDir: output directory

im_range = [-0.02 0.02];

figH = figure('Position', genfigpos(1, 'center', [1500 900]));
axH(1) = subplot(4, 3, 1);
axH(2) = subplot(4, 3, 2);
axH(3) = subplot(4, 3, 3);
axH(4) = subplot(4, 3, 4);
axH(1).Position = [0.070 0.1100 0.1 0.8150];
axH(2).Position = [0.2 0.1100 0.55 0.8150];
axH(3).Position = [0.78 0.1100 0.1 0.8150];
axH(4).Position = [0.89 0.1100 0.1 0.8150];

vartype_idx = sti_mot_ln_ca.stim_type(sti_mot_ln_ca.stim_type ~= 1);
vartype_u = unique(vartype_idx);
mean_LN_filter = squeeze(mean(LN_filter_motor, 2))';
[clus] = hierarchicalClus(zscorebigmem(mean_LN_filter), ...
    1, 'euclidean', 'ward', 3);

imagesc(mean_LN_filter(clus.lorder, vartype_idx == vartype_u(1)), 'Parent', axH(1))
caxis(axH(1), im_range)
axH(1).YDir = 'normal';
axH(1).Title.String = 'fictrac-derived speeds';

imagesc(mean_LN_filter(clus.lorder, vartype_idx ~= vartype_u(1)), 'Parent', axH(2))
caxis(axH(2), im_range)
axH(2).YDir = 'normal';
axH(2).Title.String = 'SVD components';

plot(evar_motor(clus.lorder) + var_stim_mean(clus.lorder), ...
    1:numel(clus.lorder), '.r', 'Parent', axH(3))
hold(axH(3), 'on');
plot(var_stim_mean(clus.lorder), 1:numel(clus.lorder), '.b', 'Parent', axH(3))
axH(3).YLim = [1 length(var_total)];
axH(3).YLabel.String = 'ROIs';
axH(3).XLabel.String = 'variance';
axH(3).Title.String = 'variance (\color{blue}auditory, \color{red}auditory+emotor)';

plot(evar_motor(clus.lorder), 1:numel(clus.lorder), ...
    '.r', 'Parent', axH(4))
axH(4).YLim = [1 length(var_total)];
axH(4).XLabel.String = 'variance';
axH(4).Title.String = 'variance (\color{red}emotor)';

fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;

figformat = [1 0 0 0 0 0 0 0 1 0 1];
save_edit_fig_int(axH, figH, oDir, ...
    [filenamefig, '_full_model_filters'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH)

end

function plot_full_vs_single_vartype_model(...
    var_motor_single, var_motor, evar_motor_single, ...
    evar_motor, filenamefig, oDir)
% plot_full_vs_single_vartype_model: plots variance and explained variance
%   of full vs single var type model
%
% Usage:
%   plot_full_vs_single_vartype_model(...
%       var_motor_single, var_motor, evar_motor_single, ...
%       evar_motor, filenamefig, oDir)
%
% Args:
%   var_motor_single: variance of each single model [x, model_i]
%   var_motor: variance of full model [x, 1]
%   evar_motor_single: explained variance of each single model [x, model_i]
%   evar_motor: explained variance of full model [x, 1]
%   filenamefig: filename
%   oDir: output directory

range_ = [0 max([max(var_motor_single(:)) max(var_motor(:))])*1.1];

figH(1) = figure('Position', genfigpos(1, 'center', [800 400]));
for i = 1:size(var_motor_single, 2)
    axH(i) = subplot(1, size(var_motor_single, 2), i);
    
    plot(var_motor_single(:, i), var_motor, '.k', 'Parent', axH(i))
    hold(axH(i), 'on')
    plot(range_, range_, 'r', 'Parent', axH(i))
    axH(i).XLabel.String = 'motor variance single-var type model';
    axH(i).YLabel.String = 'motor variance full model';
    
end

figH(2) = figure('Position', genfigpos(1, 'center', [800 400]));
for i = 1:size(var_motor_single, 2)
    axH_(i) = subplot(1, size(var_motor_single, 2), i);
    
    plot(evar_motor_single(:, i), evar_motor, '.k', 'Parent', axH_(i))
    hold(axH_(i), 'on')
    plot(range_, range_, 'r', 'Parent', axH_(i))
    axH_(i).XLabel.String = 'motor evariance single-var type model';
    axH_(i).YLabel.String = 'motor evariance full model';
    
end

fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;

figformat = [1 0 0 0 0 0 0 0 1 0 1];
save_edit_fig_int(axH, figH(1), oDir, ...
    [filenamefig, '_full_vs_singlevartype_model_var'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH(1))

save_edit_fig_int(axH_, figH(2), oDir, ...
    [filenamefig, '_full_vs_singlevartype_model_eVar'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH(2))

end

function plot_orthogonal_test(fullQRR, filenamefig, oDir)
% plot_orthogonal_test: function that plots qr of regressors
%
% Usage:
%   plot_orthogonal_test(fullQRR, filenamefig, oDir)
%
% Args:
%   fullQRR: ortho triangular decomposition of design regressor matrix
%   filenamefig: filename
%   oDir: output directory

figH = figure('Position', genfigpos(1, 'center', [400 400]));
axH = subplot(1, 1, 1);
plot(abs(diag(fullQRR)), 'linewidth', 2, 'Parent', axH);
title('Regressor orthogonality');
drawnow;
axis square;
ylim([0 1.1]);
ylabel('Norm. vector angle');
xlabel('Regressors');

fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;

figformat = [1 0 0 0 0 0 0 0 1 0 1];
save_edit_fig_int(axH, figH, oDir, ...
    [filenamefig, '_regressors_orthogonality'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH)

end
