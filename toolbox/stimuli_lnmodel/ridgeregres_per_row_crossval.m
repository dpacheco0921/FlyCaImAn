function [pcor_raw, lFilter, Y_pred, Y_pred_single, Y_pred_unique] = ...
    ridgeregres_per_row_crossval(...
    training_idx, Y, regress_matrix, ...
    chunksiz, corenum, ridge_algorithm, regress_type)
% ridgeregres_per_row_crossval: Run Ridge regression per raw of Y, for each
%   row of train indeces (train_idx), test indeces would be ~train_idx.
%
% Usage:
%   [pcor_raw, lFilter] = ridgeregres_per_row_crossval(...
%       training_idx, Y, regress_matrix, ...
%       chunksiz, corenum, ridge_algorithm)
%
% Args:
%   training_idx: indexes of training timepoints [n, T]
%       n: different train arragaments, T: time.
%   Y: time traces of input variables [n, T]
%   regress_matrix: matrix of regressors [T, m]
%       m: number of regressors.
%   chunksiz: number of chunks for parpool
%   corenum: number of cores
%   ridge_algorithm: algorithm to use ('bayes', 'MML')
%       (default, 'bayes')
%   regress_type: vector with labels (as integers) for different types of
%       regressors.
%
% Outputs:
%   pcor_raw: pearson correlation predicted vs raw for test indeces
%   lFilter: filter per training arragement
%   Y_pred: predicted signal (same size as Y)

if ~exist('ridge_algorithm', 'var') || isempty(ridge_algorithm)
    ridge_algorithm = 'bayes';
end

if ~exist('regress_type', 'var')
    regress_type = [];
end

siz1 = size(Y, 1);

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
        t_Y = [];
        t_Y_single = [];
        t_Y_unique = [];
        
        for iii = roi_i_l{ii, 1}
            
            if i == 1 && ii == 1 && k == 1
                t0_ = stic;
            end
            
            [t_pcor_raw(k, :), t_lFilter(:, :, k), t_Y(k, :), ...
                t_Y_single(:, :, k), t_Y_unique(:, :, k)] = ...
                ridgeregres_per_training_row(training_idx, ...
                Y(iii, :), regress_matrix, ridge_algorithm, regress_type);

            if i == 1 && ii == 1 && k == 1
                fprintf(['time per roi ', ...
                    num2str(stoc(t0_)), ' seconds\n']);
                fprintf(['Estimated time ', ...
                    num2str(numel(chunk_idx)*chunksiz*stoc(t0_)/3600), ...
                    ' hours\n']); 
            end
            
            k = k + 1;
            
        end
        
        tb_pcor_raw{ii, 1} = t_pcor_raw; 
        tb_lFilter{ii, 1} = t_lFilter;
        tb_Y{ii, 1} = t_Y;
        tb_Y_single{ii, 1} = t_Y_single;
        tb_Y_unique{ii, 1} = t_Y_unique;
        
    end
    
    pcor_raw(cat(2, roi_i_l{:}), :) = ...
        cat(1, tb_pcor_raw{:});
    lFilter(:, :, cat(2, roi_i_l{:})) = ...
        cat(3, tb_lFilter{:});
    Y_pred(cat(2, roi_i_l{:}), :) = ...
        cat(1, tb_Y{:});
    Y_pred_single(:, :, cat(2, roi_i_l{:})) = ...
        cat(3, tb_Y_single{:});
    Y_pred_unique(:, :, cat(2, roi_i_l{:})) = ...
        cat(3, tb_Y_unique{:});
    
    clear tb_pcor_raw tb_lFilter tb_Y roi_i_l ...
        t_Y tb_Y_single tb_Y_unique
    
    if mod(i, 1) == 0
        fprintf('%2.1f%% of chunks completed \n', ...
            i*100/numel(chunk_idx));
    end
    
end

end

function [LN_pcor, LN_filter, LN_Y_pred, ...
    LN_Y_pred_perm1, LN_Y_pred_perm2] = ...
    ridgeregres_per_training_row(...
    training_idx, Y, regress_matrix, ...
    ridge_algorithm, regress_type)
% ridgeregres_per_training_row: Run Ridge regression to all 
%   the possible combinations of train and test data
%
% Usage:
%   [LN_filter, LN_pcor_periter, LN_Y_pred] = ...
%       ridgeregres_per_training_row(train_idx, Y, stimM)
%
% Args:
%   training_idx: indexes of training timepoints [n, T]
%       n: different train arragaments, T: time.
%   Y: time traces [1, T]
%   regress_matrix: matrix of regressors [T, m]
%       m: number of regressors.
%   ridge_algorithm: algorithm to use ('bayes', 'MML')
%       (default, 'bayes')
%   regress_type: vector with labels (as integers) for different types of
%       regressors.
%
% Outputs:
%   LN_pcor: pearson correlation predicted vs raw for test indeces
%   LN_filter: filter per train arragement
%   LN_Y_pred: predicted trace for test indeces

if ~exist('ridge_algorithm', 'var') || isempty(ridge_algorithm)
    ridge_algorithm = 'bayes';
end

test_idx = ~training_idx; 
LN_Y_pred = zeros(size(Y)); 
LN_filter = []; 
LN_pcor = zeros(size(training_idx, 1), 1); 

if ~isempty(regress_type)
    regress_u = unique(regress_type)';
end

% run per training row
for iter_i = 1:size(training_idx, 1)
    
    % estimate filter from training indeces
    zs_Y_i = zscorebigmem(Y(1, training_idx(iter_i, :)));
    zs_regressM_i = zscorebigmem(regress_matrix(training_idx(iter_i, :), :)')';
    
    if contains(ridge_algorithm, 'bayes')
        LN_filter(:, iter_i) = runRidgeOnly(zs_regressM_i, zs_Y_i', ...
            size(zs_regressM_i, 2), 1);
    elseif contains(ridge_algorithm, 'mml')
        [~, LN_filter(:, iter_i)] = ridgeMML(zs_Y_i', zs_regressM_i, true);        
    end
        
    % calculate predicted trace and correlation 
    %   coefficient on test indeces
    zs_Y = zscorebigmem(Y(1, test_idx(iter_i, :)));
    zs_regressM = zscorebigmem(regress_matrix(test_idx(iter_i, :), :)')';
    Y_pred = (zs_regressM*LN_filter(:, iter_i))';
    LN_pcor(iter_i, 1) = corr(zs_Y', Y_pred');
    
    % generate predicted trace from test data for full model
    LN_Y_pred(1, test_idx(iter_i, :)) = Y_pred;
    
    % generate predicted trace from test data for imcomplete models
    if ~isempty(regress_type)
        
        % permute selected stim_type
        for i = 1:numel(regress_u)
            
            % shuffle
            regress2perm = zs_regressM(:, regress_type == regress_u(i))';
            siz = size(regress2perm);
            row_idx = randshuffleper(siz(2), siz(1));
            col_idx = repmat((1:siz(1))', 1, siz(2));
            regress2perm = regress2perm(sub2ind(siz, col_idx, row_idx))';
            
            % replace
            zs_regressM_2 = zs_regressM;
            zs_regressM_2(:, regress_type == regress_u(i)) = regress2perm;
            
            % get predicted
            LN_Y_pred_perm1(i, test_idx(iter_i, :)) = (zs_regressM_2*LN_filter(:, iter_i))';
            
            clear regress2perm siz row_idx col_idx zs_regressM_2
        end
        
        % permute unselected stim_type
        if numel(regress_u) > 2
            for i = 1:numel(regress_u)

                % shuffle
                regress2perm = zs_regressM(:, regress_type ~= regress_u(i))';
                siz = size(regress2perm);
                row_idx = randshuffleper(siz(2), siz(1));
                col_idx = repmat((1:siz(1))', 1, siz(2));
                regress2perm = regress2perm(sub2ind(siz, col_idx, row_idx))';

                % replace
                zs_regressM_2 = zs_regressM;
                zs_regressM_2(:, regress_type ~= regress_u(i)) = regress2perm;

                % get predicted
                LN_Y_pred_perm2(i, test_idx(iter_i, :)) = (zs_regressM_2*LN_filter(:, iter_i))';

                clear regress2perm siz row_idx col_idx zs_regressM_2
            end
        end
        
    end
       
    clear zs_Y_i zs_regressM_i zs_Y zs_regressM Y_pred
    
end

if ~exist('LN_Y_pred_perm1', 'var')
    LN_Y_pred_perm1 = nan(size(LN_Y_pred));
end

if ~exist('LN_Y_pred_perm2', 'var')
    LN_Y_pred_perm2 = nan(size(LN_Y_pred));
end

end
