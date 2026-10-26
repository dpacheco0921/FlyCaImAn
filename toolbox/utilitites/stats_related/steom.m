function sem_ = steom(Y, dim)
% steom: calculates the standart error of the mean
%
% Usage:
%   sem_ = steom(Y, dim)
%
% Args:
%   Y: n-by-m matrix
%   dim: dimension of observations (default, 1)
%
% Returns:
%   sem_: std across dim divided by sqrt(size(Y, dim)); NaNs are ignored in the
%       std but still counted in size(Y, dim) (left unchanged)

if ~exist('dim', 'var') || isempty(dim)
    dim = 1;
end

sem_ = nanstd(Y, [], dim)/sqrt(size(Y, dim));

end