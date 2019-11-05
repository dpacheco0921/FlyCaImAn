function [pcor_stat, pcor_raw, pcor_shuffle] = ...
    correct_corrcoef(pcor_raw, pcor_shuffle, prct2use)
% correct_corrcoef: correct rank correlations and get percentile
%
% Usage:
%   [pcor_stat, pcor_raw, pcor_shuffle] = ...
%       correct_corrcoef(pcor_raw, pcor_shuffle, prct2use)
%
% Args:
%   pcor_raw: correlation of raw data
%   pcor_shuffle: correlation of shuffle data
%   prct2use: percentile to use
%
% Notes: 
% Zero nans
% interpration of nans: intuitively one would think that a nan means that
% the regression did not find a suitable filter so the rank correlation of
% the predicted to raw responses should be 0.

pcor_raw(isnan(pcor_raw)) = 0;
pcor_shuffle(isnan(pcor_shuffle)) = 0;

% get a CC statistic (prct_cc)
pcor_stat = prctile(pcor_raw, prct2use, 2);

end
