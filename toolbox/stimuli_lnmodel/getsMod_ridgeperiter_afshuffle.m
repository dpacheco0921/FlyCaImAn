function [LN_pcor, LN_filter_mean, ...
    LN_filter_med, LN_filter_sd] = ...
    getsMod_ridgeperiter_afshuffle(...
    train_idx, Y, stimM, stimSiz, ...
    ifilter, perm_n, stim, filterebound)
% getsMod_ridgeperiter_afshuffle: Run Ridge regression 
%   to phase shuffled combinations
%
% Usage:
%   [LN_pcor, LN_filter_mean, ...
%       LN_filter_med, LN_filter_sd] = ...
%       getsMod_ridgeperiter_afshuffle(...
%       train_idx, Y, stimM, stimSiz, ...
%       ifilter, perm_n, stim, filterebound)
%
% Args:
%   train_idx: indexes of train timepoints [1, T],
%       train arragements, T: time
%   Y: time traces [1, T]
%   stimM: stimuli to use for prediction (at different lags)
%       [T, m], m: number of weights
%   stimSiz: stimuli size
%       (in case stimM is composed of many stimuli types)
%   ifilter: input filter to use for prediction (LN_pcor) [m, 1],
%       m: number of weights 
%   perm_n: number of permutations
%   stim: vector with stimuli structure information (binary) [1, T]
%   filterebound: get filter error bounds
%       (get filter for shuffle data, but only collects mean, median and sd)
%   
% Outputs:
%   LN_pcor: pearson correlation predicted vs raw for test indecex
%   LN_filter_mean: mean estimated filter
%   LN_filter_med: median estimated filter
%   LN_filter_sd: std of estimated filters

if ~exist('filterebound', 'var') || isempty(filterebound)
    filterebound = 0;
end

% perform aaft shuffle (shuffle phases but keeps amplitude & frequency distribution)
% it deletes arragements where stim_struct is the same
Y_shuffle = aaft(Y, perm_n, stim)';

% get Y for train and test data
test_idx = ~train_idx;
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
    
    % get filter
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
