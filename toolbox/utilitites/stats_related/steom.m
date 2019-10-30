function sem_ = steom(Y, dim)
% steom: calculates the standart error of the mean
%
% Usage:
%   sem_ = steom(Y, dim)
%
% Args:
%   Y: n-by-m matrix
%   dim: dimension of observations

if ~exist('dim', 'var') || isempty(Y)
    dim = 1;
end

sem_ = nanstd(Y, [], dim)/sqrt(size(Y, dim));

end