function pval = calculate_pval_2(...
    coef_raw, coef_shuffle, fdr, ...
    pval_cortype, prct2use)
% calculate_pval_2: generate pvals from a set of raw and shuffle
%   coefficients
%
% Usage:
%   pval = calculate_pval_2(...
%      coef_raw, coef_shuffle, fdr, ...
%       pval_cortype, prct2use)
%
% Args:
%   coef_raw: coefficient of raw data [K, 1], or [K, n] if prct2use is
%       given (one row per ROI/variable)
%   coef_shuffle: coefficient from shuffle data [K, permutations]
%       (correlation coefficient or explained variance, etc)
%   fdr: false discovery rate
%       (default, 0.01)
%   pval_cortype: type of multiple comparison correction to use, as a
%       string: 'raw' (none), 'dep' (Benjamini & Yekutieli), 'pdep'
%       (Benjamini & Hochberg) or 'bh' (Bonferroni-Holm)
%       (default, 'dep')
%   prct2use: percentile (0-100) of coef_raw across columns to use as the
%       statistic (default, [], use coef_raw as is)
%
% Returns:
%   pval: one p-value per ROI/variable [K, 1], corrected as requested.
%       Fraction of shuffle values >= the raw coefficient; negative or NaN
%       coefficients are treated as 0.

if ~exist('pval_cortype', 'var') ...
        || isempty(pval_cortype)
    pval_cortype = 'dep';
end

if ~exist('fdr', 'var') ...
        || isempty(fdr)
    fdr = 0.01;
end

if ~exist('prct2use', 'var') ...
        || isempty(prct2use)
    prct2use = [];
end

pval = [];

% 1) deal with nan (their are interpreted as failures) so
%   they are set to 0
coef_shuffle(isnan(coef_shuffle)) = 0;
coef_raw(isnan(coef_raw)) = 0;

% 2) get percentile
if ~isempty(prct2use)
    coef_raw = prctile(coef_raw, prct2use, 2);
end

% 3) compute raw p-val
%   interpretation of negative coefficient values:
%   any negative value is taken as 0 (no correlation/zero explained variance).
% one raw coefficient per row and a non-empty null distribution are
%   required, otherwise the p-values below would be meaningless or NaN
assert(numel(coef_raw) == size(coef_shuffle, 1), ...
    ['coef_raw must have one value per row of coef_shuffle (use prct2use ', ...
    'to collapse several columns)'])
assert(size(coef_shuffle, 2) > 0, 'coef_shuffle is empty')

idx_pos = coef_raw >= 0;
pval_raw = ones(numel(coef_raw), 1);

if sum(idx_pos) > 0
    pval_raw(idx_pos, 1) = sum(coef_shuffle(idx_pos, :) ...
        >= coef_raw(idx_pos), 2)./...
        sum(~isnan(coef_shuffle(idx_pos, :)), 2);
end

% 4) Correct pvalues
%   pval_corr_multi_com expects a logical vector [dep pdep bh], not a
%   string (indexing the string 'bh' as cortype(3) would error), so
%   convert the requested type into that vector.
cortype_flags = [strcmp(pval_cortype, 'dep'), ...
    strcmp(pval_cortype, 'pdep'), strcmp(pval_cortype, 'bh')];

if ~any(cortype_flags) && ~strcmp(pval_cortype, 'raw')
    error('calculate_pval_2:badCortype', ...
        'pval_cortype must be raw, dep, pdep or bh (got %s)', ...
        num2str(pval_cortype));
end

[pvalc_dep, pvalc_pdep, pvalc_bh] = ...
    pval_corr_multi_com(pval_raw, [], ...
    fdr, cortype_flags);

if strcmp(pval_cortype, 'raw')
    pval = pval_raw;
end

if strcmp(pval_cortype, 'dep')
    pval = pvalc_dep;
end

if strcmp(pval_cortype, 'pdep')
    pval = pvalc_pdep;
end

if strcmp(pval_cortype, 'bh')
    pval = pvalc_bh;
end

end
