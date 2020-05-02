function Y_pred = get_predicted_signal(training_idx, ...
    LN_filter, regress_matrix)
% get_predicted_signal: get predicted signal from training/testing idx, and
%   calculated filters
%
% Usage:
%   Y_pred = get_predicted_signal(training_idx, ...
%       LN_filter, regress_matrix)
%
% Args:
%   training_idx: indexes of training timepoints [n, T]
%   LN_filter: filter per train arragement [m, n, x]
%   regress_matrix: matrix of regressors [T, m]
%       m: number of regressors.
%       n: different train arragaments
%       T: time.
%       x: output variables to generate.

testing_idx = ~training_idx; 
Y_pred = [];

for iter_i = 1:size(training_idx, 1)
    zs_regressM = zscorebigmem(regress_matrix(testing_idx(iter_i, :), :)')';
    filter_in = squeeze(LN_filter(:, iter_i, :));
    Y_pred(:, testing_idx(iter_i, :)) = ...
        (zs_regressM*filter_in)';
end

end
