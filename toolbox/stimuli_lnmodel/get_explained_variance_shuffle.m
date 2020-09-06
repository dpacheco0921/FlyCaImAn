function [eVar_cs, corrcoef_cs, lFilter_mean, ...
    lFilter_med, lFilter_sd] = ...
    get_explained_variance_shuffle(...
    Y, Y_pred, lgate, stim_bin, redo, ...
    filename, chunk_minsize, chunk_splitn, ...
    chunksiz, corenum, btn, shuffle2use, ...
    add_stim_lag, add_inverse, train_idx, ...
    lnfit_flag, regress_matrix, ridge_algorithm)
% runridgeshuffle_chunks: Run Ridge regression to all 
%   the possible combinations of train and test data
%
% Usage:
%   eVar_cs = get_explained_variance_shuffle(...
%   	Y, Y_pred, lgate, stim_bin, redo, ...
%   	filename, chunk_minsize, chunk_splitn, ...
%   	chunksiz, corenum, btn, shuffle2use, ...
%       add_stim_lag, add_inverse)
%
% Args:
%   Y: time traces of variables [n, T]
%   Y_pred: time traces of predicted variables [n, T]
%   lgate: flag to generate permuted data
%   stim_bin: binary vector of stimuli
%   redo: redo flag
%   filename: name of temporary file to generate
%   chunk_minsize: minimun size of time chunk
%   chunk_splitn: number of times to split catrace_in in the time domain
%   chunksiz: number of chunks for parpool
%   corenum: number of cores
%   btn: number of permutations
%   shuffle2use: type of shuffling to use
%       (default, 0, random chunks, requires chunk_minsize & chunk_splitn & btn);
%       (1, random chunks + circular permutation);
%       (2, circular shuffling);
%   add_stim_lag: aditional lags of stimuli pattern 
%       to avoid (in timestamps) when doing circular permutation
%       (default, [-10 10]);
%   add_inverse: include reverse order when doing circular permutation
%       (default, 1);
%   train_idx: indeces of training timepoints (assumes test_idx = ~train_idx)
%   lnfit_flag: flag to do linear fitting on shuffle data
%   regress_matrix: matrix of regressors [T, m]
%       m: number of regressors.
%   ridge_algorithm: algorithm to use ('bayes', 'MML')
%       (default, 'bayes')
%
% Outputs:
%   eVar_cs: explained variance per permutation for each ROI
%   corrcoef_cs: correlation coefficient per permutation for each ROI
%   lFilter_mean: mean estimated filters
%   lFilter_med: median estimated filters
%   lFilter_sd: std of estimated filters

if ~exist('stim_bin', 'var')
    stim_bin = [];
end

if ~exist('shuffle2use', 'var')
    shuffle2use = 0;
end

if ~exist('add_stim_lag', 'var')
    add_stim_lag = [-10 10];
end

if ~exist('filename', 'var') || isempty(filename)
    filename = [pwd, 'temp_file.mat'];
end

if ~exist('add_inverse', 'var') || isempty(add_inverse)
    add_inverse = 1;
end

if ~exist('train_idx', 'var')
    train_idx = [];
end

if ~exist('lnfit_flag', 'var')
    lnfit_flag = false;
end

if ~exist('ridge_algorithm', 'var') || ...
        isempty(ridge_algorithm)
    ridge_algorithm = 'bayes';
end

if ~exist('regress_matrix', 'var')
    regress_matrix = [];
end

siz1 = size(Y, 1);
dataObj = matfile(filename, 'Writable', true);

if lgate || redo
    
    % initialize parameters
    dataObj.eVar_cs = [];

    % generate shuffle in chunks over time
    if shuffle2use ~= 2
        [~, rperm_chunkshufIdx] = randchunkper(...
        Y(1, :), chunk_splitn, ...
        round(chunk_minsize), btn, stim_bin);
    end
    
    % add circular permutations
    if shuffle2use ~= 0
        rperm_circIdx = randcirshuffleper(...
            numel(Y(1, :)), 10, [], stim_bin, ...
            add_stim_lag, add_inverse);
    end

    if shuffle2use == 0
        rperm_chunkIdx = rperm_chunkshufIdx;
    elseif shuffle2use == 2
        rperm_chunkIdx = rperm_circIdx;
    else
        rperm_chunkIdx = [rperm_chunkshufIdx; rperm_circIdx];
    end
    
    dataObj.rperm_chunkIdx = rperm_chunkIdx;

    % make chunks to run
    [~, ~, chunk_idx] = ppool_makechunks(...
        chunksiz, corenum, siz1);
    
else
    
    % update initial index
    vect_init = size(dataObj.eVar_cs, 1) + 1;

    % Get shuffle time
    rperm_chunkIdx = dataObj.rperm_chunkIdx;

    % make chunks to run
    [~, ~, chunk_idx] = ppool_makechunks(...
        chunksiz, corenum, siz1, vect_init);  
    
end

btn = size(rperm_chunkIdx, 1);

for i = 1:numel(chunk_idx)
    
    if i == 1
        t0 = stic;
    end
    
    batch2run = chunk_idx{i};
    
    parfor ii = 1:numel(batch2run)
        
        idx_i{ii, 1} = (batch2run(ii):min(batch2run(ii)...
            + chunksiz - 1, siz1));
        k = 1;
        t_eVar = [];
        t_corrcoef = [];
        t_fmean_cs_i = [];
        t_fmed_cs_i = [];
        t_fsd_cs_i = [];
        
        for iii = idx_i{ii, 1}
            
            if i == 1 && ii == 1 && k == 1
                t0_ = stic;
            end

            [t_eVar(k, :), t_corrcoef(k, :), t_fmean_cs_i(:, k), ...
                t_fmed_cs_i(:, k), t_fsd_cs_i(:, k)] = ...
                get_explained_variance_per_row(...
                Y(iii, :), Y_pred(iii, :), rperm_chunkIdx, ...
                train_idx, lnfit_flag, regress_matrix, ...
                ridge_algorithm);

            if i == 1 && ii == 1 && k == 1
                fprintf(['time per roi ', ...
                    num2str(stoc(t0_)), ' seconds\n']); 
                fprintf(['Estimated time ', ...
                    num2str(numel(chunk_idx)*chunksiz*stoc(t0_)/3600), ...
                    ' hours\n']); 
            end
        
            k = k + 1;
       end
        
        tb_eVar_cs_i{ii, 1} = t_eVar;
        tb_corrcoef_cs_i{ii, 1} = t_corrcoef;
        tb_fmean_cs_i{ii, 1} = t_fmean_cs_i; 
        tb_fmed_cs_i{ii, 1} = t_fmed_cs_i;
        tb_fsd_cs_i{ii, 1} = t_fsd_cs_i;
        
    end
    
    dataObj.eVar_cs(cat(2, idx_i{:}), 1:btn) = ...
        cat(1, tb_eVar_cs_i{:});
    dataObj.corrcoef_cs(cat(2, idx_i{:}), 1:btn) = ...
        cat(1, tb_corrcoef_cs_i{:});
    
    if lnfit_flag
        dataObj.lFilter_cs_mean(1:size(regress_matrix, 2), cat(2, idx_i{:})) = ...
            cat(2, tb_fmean_cs_i{:});
        dataObj.lFilter_cs_med(1:size(regress_matrix, 2), cat(2, idx_i{:})) = ...
            cat(2, tb_fmed_cs_i{:}); 
        dataObj.lFilter_cs_sd(1:size(regress_matrix, 2), cat(2, idx_i{:})) = ...
            cat(2, tb_fsd_cs_i{:});
    end
    
    clear tb_eVar_cs_i tb_corrcoef_cs_i ...
        idx_i tb_fmean_cs_i ...
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

eVar_cs = dataObj.eVar_cs;
corrcoef_cs = dataObj.corrcoef_cs;

if lnfit_flag
    lFilter_mean = dataObj.lFilter_cs_mean;
    lFilter_med = dataObj.lFilter_cs_med;
    lFilter_sd = dataObj.lFilter_cs_sd;
else
    lFilter_mean = [];
    lFilter_med = [];
    lFilter_sd = [];   
end

end

function [eVar, corrcoef, LN_filter_mean, ...
    LN_filter_med, LN_filter_sd] = get_explained_variance_per_row(...
    Y, Y_pred, rperm_tIdx, train_idx, lnfit_flag, regress_matrix, ...
    ridge_algorithm)
% get_explained_variance_per_row: get explained variance from many shuffle
%   iterations (as defined in rperm_tIdx)
%
% Usage:
%   eVar = get_explained_variance_per_row(...
%      Y, Y_pred, rperm_tIdx)
%
% Args:
%   Y: time traces [1, T]
%   Y_pred: predicted time traces [1, T]
%   rperm_tIdx: set of permutations to apply to Y [n, T],
%       n: permutations
%   train_idx: indeces of training timepoints (assumes test_idx = ~train_idx)
%   lnfit_flag: flag to do linear fitting on shuffle data
%   regress_matrix: matrix of regressors [T, m]
%       m: number of regressors.
%   ridge_algorithm: algorithm to use ('bayes', 'MML')
%       (default, 'bayes')
%
% Outputs:
%   eVar: explained variance per permutation [n, T]
%   corrcoef: correlation coefficient [n, T]
%   LN_filter_mean: mean estimated filters
%   LN_filter_med: median estimated filters
%   LN_filter_sd: std of estimated filters

if ~exist('train_idx', 'var') || isempty(train_idx)
    train_idx = false([1 size(Y, 2)]);
end

if ~exist('lnfit_flag', 'var')
    lnfit_flag = false;
end

if ~exist('ridge_algorithm', 'var') || ...
        isempty(ridge_algorithm)
    ridge_algorithm = 'bayes';
end

if ~exist('regress_matrix', 'var')
    regress_matrix = [];
end

test_idx = ~train_idx;

T = size(Y, 2);

if size(rperm_tIdx, 1) == T
    rperm_tIdx = rperm_tIdx';
end

% apply shuffle to Y
zs_Y_shuffle = zscorebigmem(Y(rperm_tIdx));
Y_pred = zscorebigmem(Y_pred);

% estimate explained variance and correlation per permutation
corrcoef = corr(zscorebigmem(zs_Y_shuffle(:, test_idx))', Y_pred(:, test_idx)');
eVar = (corrcoef').^2;

% estimate filters from shuffle data
if lnfit_flag
    
    zs_Y_shuffle = zscorebigmem(zs_Y_shuffle(:, train_idx))';
    zs_regressM = zscorebigmem(regress_matrix(train_idx, :)')';
    perm_n = size(rperm_tIdx, 1);

    for i = 1:perm_n
        
        if contains(ridge_algorithm, 'bayes')
            LN_filter_shuffle(:, i) = runRidgeOnly(zs_regressM, ...
                zs_Y_shuffle(:, i), size(zs_regressM, 2), 1);
        elseif contains(ridge_algorithm, 'mml')
            [~, LN_filter_shuffle(:, i)] = ...
                ridgeMML(zs_Y_shuffle(:, i), zs_regressM, true);        
        end
        
    end

    % collect filter stats
    LN_filter_mean = mean(LN_filter_shuffle, 2);
    LN_filter_med = prctile(LN_filter_shuffle, 50, 2);
    LN_filter_sd = std(LN_filter_shuffle, [], 2);

else
    
    LN_filter_mean = nan([1 10]);
    LN_filter_med = nan([1 10]);
    LN_filter_sd = nan([1 10]);
    
end

end
