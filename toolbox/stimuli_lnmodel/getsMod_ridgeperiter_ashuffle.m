function [LN_pcor, LN_filter_mean, ...
    LN_filter_med, LN_filter_sd] = ...
    getsMod_ridgeperiter_ashuffle(...
    train_idx, Y, stimM, stimSiz, ...
    rperm_tIdx, ifilter, filterebound)
% getsMod_ridgeperiter_ashuffle: Run Ridge regression to 
%   all *rperm_tIdx* shuffled combinations 
%   (this could be random or chunk-random, depends on how rperm_tIdx is generated)
%
% Usage:
%   [LN_pcor, LN_filter_mean, ...
%       LN_filter_med, LN_filter_sd] = ...
%       getsMod_ridgeperiter_ashuffle(...
%       train_idx, Y, stimM, stimSiz, ...
%       rperm_tIdx, ifilter, filterebound)
%
% Args:
%   train_idx: indexes of train timepoints [1, T],
%       train arragements, T: time
%   Y: time traces [1, T]
%   stimM: stimuli to use for prediction (at different lags)
%       [T, m], m: number of weights
%   stimSiz: stimuli size
%       (in case stimM is composed of many stimuli types)
%   rperm_tIdx: set of permutations to apply to Y [n, T],
%       n: permutations
%   ifilter: input filter to use for prediction (LN_pcor) [m, 1],
%       m: number of weights
%   filterebound: filter error bounds 
%       (get filter for shuffle data, only collects mean, median and sd)
%
% Outputs:
%   LN_pcor: pearson correlation predicted vs raw for test indecex
%   LN_filter_mean: mean estimated filter
%   LN_filter_med: median estimated filter
%   LN_filter_sd: std of estimated filters

if ~exist('filterebound', 'var') || ...
        isempty(filterebound)
    filterebound = 0;
end

T = size(Y, 2);

if size(rperm_tIdx, 1) == T
    rperm_tIdx = rperm_tIdx';
end
perm_n = size(rperm_tIdx, 1);

test_idx = ~train_idx;

% apply shuffle to Y
Y_shuffle = Y(rperm_tIdx);

% correct Y for train and test data (zscore)
zs_Y_shuffle_train = zscorebigmem(Y_shuffle(:, train_idx))';
zs_Y_shuffle_test = zscorebigmem(Y_shuffle(:, test_idx))';

% get stim for train and test data
zs_stim_train = zscorebigmem(stimM(train_idx, :)')';
zs_stim_test = zscorebigmem(stimM(test_idx, :)')';

% initialize variables
LN_pcor = zeros(perm_n, 1); 
LN_filter_shuffle = zeros(stimSiz, perm_n);

% estimate filters and calculate LN_pcor of shuffle to input filter
for i = 1:perm_n
    
    % get estimated filter from shuffle data
    if filterebound
        LN_filter_shuffle(:, i) = ...
            runRidgeOnly(zs_stim_train, ...
            zs_Y_shuffle_train(:, i), stimSiz, 1);
    end
    
    % get correlation coefficient
    LN_pcor(i) = corr(zs_Y_shuffle_test(:, i), ...
        (zs_stim_test*ifilter));
    
end

% collect filter stats
LN_filter_mean = mean(LN_filter_shuffle, 2);
LN_filter_med = prctile(LN_filter_shuffle, 50, 2);
LN_filter_sd = std(LN_filter_shuffle, [], 2);
clear LN_filter_shuffle

end
