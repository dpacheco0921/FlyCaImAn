function [pval_raw, pvalc_dep, ...
    pvalc_pdep, pvalc_bh, pvalks] = ...
    calculate_pval(pcor_stat, ...
    pcor_raw, pcor_shuffle, fdr, ...
    cortype, chunksiz, corenum)
% calculate_pval: generate raw and corrected
%   (for multiple comparisons) pvals from 2 sets of correlation
%
% Usage:
%   [pval_raw, pvalc_dep, ...
%      pvalc_pdep, pvalc_bh, pvalks] = ...
%      calculate_pval(pcor_stat, ...
%      pcor_raw, pcor_shuffle, fdr, ...
%      cortype, chunksiz, corenum)
%
% Args:
%   pcor_stat: percentile of pcor_raw 
%   pcor_raw: correlation of raw data
%   pcor_shuffle: correlation of shuffle data
%   fdr: false discovery rate
%   cortype: type of corrections to perform
%   chunksiz: number of chunks for parpool
%   	(default, 2*10^3)
%   corenum: number of cores
%   	(default, 4)

if ~exist('cortype', 'var') ...
        || isempty(cortype)
    cortype = zeros(1, 4);
end

if ~exist('fdr', 'var') ...
        || isempty(fdr)
    fdr = 0.01;
end

if ~exist('chunksiz', 'var') ...
        || isempty(chunksiz)
    chunksiz = 2*10^3;
end

if ~exist('corenum', 'var') ...
        || isempty(corenum)
    corenum = 4;
end

pval_raw = [];
pvalc_dep = [];
pvalc_pdep = [];
pvalc_bh = [];
pvalks = [];

% 1) compute raw p-val
%   interpretation of negative CC,
%   during the regression the filter is
%   suppose to have a positive CC,
%   so any negative value is taken as 0.
idx_pos = pcor_stat >= 0;
pval_raw = ones(numel(pcor_stat), 1);

if sum(idx_pos) > 0
    pval_raw(idx_pos, 1) = sum(pcor_shuffle(idx_pos, :) ...
        >= pcor_stat(idx_pos), 2)./...
        sum(~isnan(pcor_shuffle(idx_pos, :)), 2);
end

% so far negative values do not have 
%   any meaning so I dont know what to do about them.

% 2) Correct pvalues
[pvalc_dep, pvalc_pdep, pvalc_bh] = ...
    pval_corr_multi_com(pval_raw, [], ...
    fdr, cortype);

% 3) alternative pvalue, compute ks test
if cortype(4)
    
    pvalks = [];
    K = size(pcor_raw, 1);
    
    [~, ~, chunk_idx] = ...
        ppool_makechunks(chunksiz, corenum, K);
    chunk_idx = chunk_idx{1};
    
    parfor i = 1:numel(chunk_idx)
        
        i_idx{i, 1} = ...
            (chunk_idx(i):min(chunk_idx(i) ...
            + chunksiz - 1, K));
        k = 1;
        pvalks_i = [];
        
        for ii = i_idx{i, 1}
            [~, pvalks_i(k, :)] = ...
                kstest2(pcor_shuffle(ii, :), ...
                pcor_raw(ii, :));
            k = k + 1;
        end
        
        pvalks_i_c{i, 1} = pvalks_i;
        
    end
    
    pvalks(cell2mat(i_idx')) = cell2mat(pvalks_i_c);
    
end

end
