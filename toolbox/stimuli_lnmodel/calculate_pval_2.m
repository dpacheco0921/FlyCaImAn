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
%   coef_raw: coefficient of raw data 
%   coef_shuffle: coefficient from shuffle data
%       (correlation coefficient or explained variance, etc)
%   fdr: false discovery rate
%       (default, 0.01)
%   pval_cortype: type of multiple comparison correction to use: dep, pdep, bh)
%       (default, 'dep')
%   prct2use: percentile of coef_raw to use

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
idx_pos = coef_raw >= 0;
pval_raw = ones(numel(coef_raw), 1);

if sum(idx_pos) > 0
    pval_raw(idx_pos, 1) = sum(coef_shuffle(idx_pos, :) ...
        >= coef_raw(idx_pos), 2)./...
        sum(~isnan(coef_shuffle(idx_pos, :)), 2);
end

% 2) Correct pvalues
[pvalc_dep, pvalc_pdep, pvalc_bh] = ...
    pval_corr_multi_com(pval_raw, [], ...
    fdr, pval_cortype);

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
