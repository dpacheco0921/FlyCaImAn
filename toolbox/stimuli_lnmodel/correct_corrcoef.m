function [coef_stat, coef_raw, coef_shuffle] = ...
    correct_corrcoef(coef_raw, coef_shuffle, prct2use)
% correct_corrcoef: correct coefficients (correlation/explained variance/etc)
%   and get percentile
%
% Usage:
%   [coef_stat, coef_raw, coef_shuffle] = ...
%       correct_corrcoef(coef_raw, coef_shuffle, prct2use)
%
% Args:
%   coef_raw: coefficient from raw data [K, folds]
%       (correlation coefficient or explained variance, etc)
%   coef_shuffle: coefficient from shuffle data [K, permutations]
%       (correlation coefficient or explained variance, etc)
%   prct2use: percentile (0-100) of coef_raw across columns used as the
%       statistic (required)
%
% Returns:
%   coef_stat: percentile of coef_raw per row [K, n_percentiles]
%   coef_raw: input with NaN replaced by 0
%   coef_shuffle: input with NaN replaced by 0
%
% Notes: 
% interpretation of nans (nans are zeroed): it assumes that nan means that
%   the regression did not find a suitable filter so the coefficient
%   measured of the predicted to raw responses should be 0.

coef_raw(isnan(coef_raw)) = 0;
coef_shuffle(isnan(coef_shuffle)) = 0;

% get a statistic (coef_stat)
coef_stat = prctile(coef_raw, prct2use, 2);

end
