function [LN_filter, Y_pred] = ridgeregres_per_row(...
    Y, regress_matrix, ridge_algorithm, verbose)
% ridgeregres_per_row: Run Ridge regression for each row independently
%
% Usage:
%   [LN_filter, Y_pred] = ridgeregres_per_row(...
%       Y, regress_matrix, ridge_algorithm, verbose)
%
% Args:
%   Y: time traces [n, T]
%   regress_matrix: matrix of regressors [T, m]
%       (m: number of regressors)
%   ridge_algorithm: algorithm to use ('bayes', 'MML'; case insensitive)
%       (default, 'bayes')
%   verbose: 1 to print timing and progress (default, 0)
%
% Outputs:
%   LN_filter: ridge filter (one weight per regressor, applied to
%       z-scored regressors) for each row of Y [m, n]
%   Y_pred: predicted trace per row of Y, in z-score units (the filter is
%       fit on all timepoints, no train/test split) [n, T]

if ~exist('verbose', 'var') || isempty(verbose)
    verbose = 0;
end

if ~exist('ridge_algorithm', 'var') || isempty(ridge_algorithm)
    ridge_algorithm = 'bayes';
end

LN_filter = []; 
chunk_n = 400;

for i = 1:size(Y, 1)
    
    if i == 1
    	t0_ = stic;
    end
    
    % estimate filter
    zs_Y_i = zscorebigmem(Y(i, :));
    zs_regressM_i = zscorebigmem(regress_matrix')';
    
    if contains(lower(ridge_algorithm), 'bayes')
        LN_filter(:, i) = runRidgeOnly(zs_regressM_i, zs_Y_i', ...
            size(zs_regressM_i, 2), 1);
    elseif contains(lower(ridge_algorithm), 'mml')
        [~, LN_filter(:, i)] = ridgeMML(zs_Y_i', zs_regressM_i, true);
    else
        error('Unknown ridge_algorithm: %s (use bayes or MML)', ...
            ridge_algorithm);
    end
        
    clear zs_Y_i zs_stim_i
    
    if i == 1 && verbose
        fprintf(['time per roi ', ...
            num2str(stoc(t0_)), ' seconds\n']); 
    end
    
    if mod(i, chunk_n) == 0 && verbose
        fprintf('%2.1f%% of chunks completed \n', ...
            i*100/size(Y, 1));
    end
    
end

% generate predicted trace
Y_pred = (zscorebigmem(regress_matrix')'*LN_filter)';

end
