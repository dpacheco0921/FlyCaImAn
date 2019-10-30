function rperm_tIdx = randshuffleper(...
    timepoints_n, repn, stim)
% randshuffleper: Generate a Random permutation 
%   excluding the arragements that are the same as the
%   original or have similar stim structure as stim
%
% Usage:
%   rperm_tIdx = randshuffleper(...
%      timepoints_n, repn, stim)
%
% Args:
%   timepoints_n: # of timepoints
%   repn: number of permutations
%   stim: vector with stimuli information
%
% Output:
%   rperm_tIdx: matrix with indeces of random permutation per rows

if ~exist('stim', 'var') || isempty(stim)
    stim = [];
end

rng('shuffle');

rperm_tIdx = zeros(repn, timepoints_n);

for s_i = 1:repn
    
    i_gate = 0;
    
    while i_gate == 0
        
        rperm_tIdx(s_i, :) = randperm(timepoints_n);
        
        if ~isempty(stim)
            
            if ~ismember(rperm_tIdx(s_i, :), ...
                    1:timepoints_n, 'rows') ...
                    && ~ismember(stim(s_i), stim, 'rows')
                i_gate = 1; 
            end
            
        else
            
            if ~ismember(rperm_tIdx(s_i, :), ...
                    1:timepoints_n, 'rows')
                i_gate = 1; 
            end
            
        end
        
    end
    
end

end
