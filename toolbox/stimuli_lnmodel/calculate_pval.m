function [pval_raw, pvalc_dep, ...
    pvalc_pdep, pvalc_bh, pvalks] = ...
    calculate_pval(coef_stat, ...
    coef_raw, coef_shuffle, fdr, ...
    pval_cortype, chunksiz, corenum)
% calculate_pval: generate raw and corrected
%   (for multiple comparisons) pvals from 2 sets of correlation
%
% Usage:
%   [pval_raw, pvalc_dep, ...
%      pvalc_pdep, pvalc_bh, pvalks] = ...
%      calculate_pval(coef_stat, ...
%      coef_raw, coef_shuffle, fdr, ...
%      pval_cortype, chunksiz, corenum)
%
% Args:
%   coef_stat: statistic per ROI (percentile of coef_raw) [K, 1]
%   coef_raw: coefficient from raw data [K, folds]
%       (correlation coefficient or explained variance, etc)
%   coef_shuffle: coefficient from shuffle data [K, permutations]
%       (correlation coefficient or explained variance, etc)
%   fdr: false discovery rate (default, 0.01)
%   pval_cortype: logical vector [1, 4] choosing the corrections to run:
%       [dep pdep bh ks] (default, zeros(1, 4)). dep: Benjamini &
%       Yekutieli, pdep: Benjamini & Hochberg, bh: Bonferroni-Holm,
%       ks: KS test p-value (only the first corenum*chunksiz ROIs are
%       computed, see Notes)
%   chunksiz: number of chunks for parpool
%   	(default, 2*10^3)
%   corenum: number of cores
%   	(default, 4)
%
% Returns:
%   pval_raw: uncorrected p-value per ROI [K, 1] (fraction of shuffle
%       values >= coef_stat; negative coef_stat gives 1)
%   pvalc_dep, pvalc_pdep, pvalc_bh: corrected p-values [K, 1]
%       (NaN for corrections not requested in pval_cortype)
%   pvalks: KS-test p-values, empty unless pval_cortype(4) is set
%
% Notes:
%   The KS branch only processes chunk_idx{1} (the first batch of
%   corenum chunks), so ROIs beyond corenum*chunksiz are left unfilled.

if ~exist('pval_cortype', 'var') ...
        || isempty(pval_cortype)
    pval_cortype = zeros(1, 4);
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
%   interpretation of negative coefficient values:
%   any negative value is taken as 0 (no correlation/zero explained variance).
idx_pos = coef_stat >= 0;
pval_raw = ones(numel(coef_stat), 1);

if sum(idx_pos) > 0
    pval_raw(idx_pos, 1) = sum(coef_shuffle(idx_pos, :) ...
        >= coef_stat(idx_pos), 2)./...
        sum(~isnan(coef_shuffle(idx_pos, :)), 2);
end

% 2) Correct pvalues
[pvalc_dep, pvalc_pdep, pvalc_bh] = ...
    pval_corr_multi_com(pval_raw, [], ...
    fdr, pval_cortype);

% 3) alternative pvalue, compute ks test
if pval_cortype(4)
    
    pvalks = [];
    K = size(coef_raw, 1);
    
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
                kstest2(coef_shuffle(ii, :), ...
                coef_raw(ii, :));
            k = k + 1;
        end
        
        pvalks_i_c{i, 1} = pvalks_i;
        
    end
    
    pvalks(cell2mat(i_idx')) = cell2mat(pvalks_i_c);
    
end

end
