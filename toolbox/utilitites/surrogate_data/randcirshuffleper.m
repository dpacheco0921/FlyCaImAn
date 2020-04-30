function rperm_tIdx = randcirshuffleper(...
    timepoints_n, min_shift, repn, stim, stim_lags2avoid)
% randcirshuffleper: generate a set of circular permutations
%
% Usage:
%   rperm_tIdx = randcirshuffleper(...
%      timepoints_n, min_shift, repn, stim)
%
% Args:
%   timepoints_n: # of timepoints
%   min_shift: min shift length
%   repn: number of permutations
%   stim: vector with stimuli information
%   stim_lags2avoid: set of lags to remove from circ shuffle
%       (default, [-4 4])
%
% Output:
%   rperm_tIdx: matrix with indeces of random permutation per rows

if ~exist('min_shift', 'var') || isempty(min_shift)
    min_shift = 10;
end

if ~exist('stim', 'var') || isempty(stim)
    stim = [];
end

if ~exist('stim_lags2avoid', 'var') || isempty(stim_lags2avoid)
    stim_lags2avoid = [-4 4];
end

rng('shuffle');

max_circ_perm = randperm(timepoints_n);
max_circ_perm = max_circ_perm(max_circ_perm > min_shift);

rperm_tIdx = zeros(numel(max_circ_perm), timepoints_n);

for s_i = 1:numel(max_circ_perm)
    rperm_tIdx(s_i, :) = circshift(1:timepoints_n, max_circ_perm(s_i));
end

% remove permutations similar to stimulus
if ~isempty(stim)
    
    if ~isempty(stim_lags2avoid)
        stim_ = [];
        shifts2use = [1:stim_lags2avoid(2), ...
            numel(max_circ_perm)+stim_lags2avoid(1):numel(max_circ_perm)];
        for i = shifts2use
            stim_(shifts2use == i, :) = stim(circshift(1:timepoints_n, i));
        end
    else
        stim_ = stim;
    end
    
    perm2rem_ = find(ismember(rperm_tIdx(s_i, :), ...
            1:timepoints_n, 'rows'));
    rperm_tIdx(perm2rem_, :) = [];    
        
    perm2rem_ = find(ismember(stim(rperm_tIdx), stim_, 'rows'));
    rperm_tIdx(perm2rem_, :) = [];    
    
end

% subsampling permutations
if exist('repn', 'var') && ~isempty(repn)
    
    if repn > size(rperm_tIdx, 1)
        fprintf('requested number of permutations above max possible')
    else
        rperm_tIdx = rperm_tIdx(1:repn, :);
    end
    
end

end
