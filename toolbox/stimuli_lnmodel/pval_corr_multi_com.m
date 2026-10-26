function [pvalc_dep, pvalc_pdep, pvalc_bh] = ...
    pval_corr_multi_com(pval_input, file_idx, ...
    fdr, cortype)
% pval_corr_multi_com: correct pvalues to 
%   account for multiple comparisons
%
% Usage:
%   [pvalc_dep, pvalc_pdep, pvalc_bh] = ...
%       pval_corr_multi_com(pval_input, file_idx, ...
%       fdr, cortype)
%
% Args:
%   pval_input: input pvalues
%   file_idx: indeces of files where data comes from
%   fdr: false discovery rate
%   cortype: logical/numeric vector [1, 3] choosing which corrections to
%       compute (default, zeros(1, 3)):
%       cortype(1): Benjamini & Yekutieli (dep)
%       cortype(2): Benjamini & Hochberg FDR (pdep)
%       cortype(3): bonferroni-holm test (bh)
%
% Returns:
%   pvalc_dep, pvalc_pdep, pvalc_bh: corrected p-values, same size as
%       file_idx (pval_input). NaN for corrections that were not requested.

if ~exist('file_idx', 'var') || isempty(file_idx)
    file_idx = ones(size(pval_input));
end

if ~exist('cortype', 'var') || isempty(cortype)
    cortype = zeros(1, 3);
end

if ~exist('fdr', 'var') || isempty(fdr)
    fdr = 0.01;
end

pvalc_dep = nan(size(file_idx));
pvalc_pdep = nan(size(file_idx));
pvalc_bh = nan(size(file_idx));

for i = unique(file_idx)'
    
    idx2use = file_idx == i;
    
    % 1) compute Benjamini & Yekutieli (dep)
    if cortype(1)
        [~, ~, ~, pvalc_dep(idx2use)] = ...
            fdr_bh(pval_input(idx2use), fdr, 'dep', 'no');
    end

    % 2) compute Bejnamini & Hochberg FDR (pdep)
    if cortype(2)
        [~, ~, ~, pvalc_pdep(idx2use)] = ...
            fdr_bh(pval_input(idx2use), fdr, 'pdep', 'no');
    end

    % 3) compute bonferroni-holm test
    if cortype(3)
        [pvalc_bh(idx2use), ~] = ...
            bonf_holm(pval_input(idx2use), fdr);
    end
    
end

end
