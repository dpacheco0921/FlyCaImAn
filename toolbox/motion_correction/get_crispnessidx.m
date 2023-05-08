function crispnessidx = get_crispnessidx(Y)
% get_crispnessidx: calculate crispness index as in
%   (http://dx.doi.org/10.1016/j.jneumeth.2017.07.031)
%
% Usage:
%   crisp_idx = get_crispnessidx(Y)
%
% Args:
%   Y: input 3D or 2D matrix
%
% Returns:
%   crispnessidx: crispness index

crispnessidx = [];

siz_ = size(Y);
pre_matrix = abs(reshape(squeeze(mean(Y, length(siz_))),...
    [prod(siz_(1:end-1)) 1]));

crispnessidx = norm(pre_matrix, 'fro');

end