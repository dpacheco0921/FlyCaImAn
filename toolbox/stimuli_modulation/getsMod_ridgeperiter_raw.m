function [LN_pcor, LN_filter] = ...
    getsMod_ridgeperiter_raw(...
    train_idx, Y, stimM, stimSiz)
% getsMod_ridgeperiter_raw: Run Ridge regression to all 
%   the possible combinations of train and test data
%
% Usage:
%   [LN_filter, LN_pcor_periter] = ...
%       getsMod_ridgeperiter_raw(...
%       train_idx, Y, stimM, stimSiz)
%
% Args:
%   train_idx: indexes of train timepoints [n, T]
%       n: different train arragaments, T: time.
%   Y: time traces [1, T]
%   stimM: stimuli to use for prediction (at different lags) [T, m]
%       m: number of weights.
%   stimSiz: stimuli size (in case stimM is composed of many stimuli types)
%
% Outputs:
%   LN_pcor: pearson correlation predicted vs raw for test indeces
%   LN_filter: filter per train arragement

test_idx = ~train_idx; 
LN_filter = []; 
LN_pcor = zeros(size(train_idx, 1), 1); 
zs_stim = [];

% Estimate model from train indeces
for iter_i = 1:size(train_idx, 1)
    
    % estimate filter
    zs_Y_i = zscorebigmem(Y(1, train_idx(iter_i, :)));
    zs_stim_i = zscorebigmem(stimM(train_idx(iter_i, :), :)')';
    LN_filter(:, iter_i) = runRidgeOnly(zs_stim_i, zs_Y_i', stimSiz, 1);
    
    % Test model on test indeces
    
    zs_Y = zscorebigmem(Y(1, test_idx(iter_i, :)));
    zs_stim = zscorebigmem(stimM(test_idx(iter_i, :), :)')';
    LN_pcor(iter_i, 1) = corr(zs_Y', (zs_stim*LN_filter(:, iter_i)));
    
    clear zs_Y_i zs_stim_i zs_Y zs_stim
    
end

end
