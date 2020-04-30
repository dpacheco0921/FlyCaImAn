function Y_pred = get_predicted_signal(training_idx, LN_filter, regress_matrix)
% get_predicted_signal: get predicted signal from training/testing idx, and
%   calculated filters
%
% Usage:
%   Y_pred = get_predicted_signal(training_idx, LN_filter, regress_matrix)
%
% Args:
%   training_idx: indexes of training timepoints [n, T]
%   LN_filter: filter per train arragement [m, n, x]
%   regress_matrix: matrix of regressors [T, m]
%       m: number of regressors.
%       n: different train arragaments
%       T: time.
%       x: output variables to generate.

test_idx = ~training_idx; 

for iter_i = 1:size(training_idx, 1)
    zs_regressM = zscorebigmem(regress_matrix(test_idx(iter_i, :), :)')';
    Y_pred(:, test_idx(iter_i, :)) = (zs_regressM*LN_filter(:, iter_i, :))';
end

end
