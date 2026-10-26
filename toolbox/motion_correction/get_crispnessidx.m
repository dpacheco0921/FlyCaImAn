function crispnessidx = get_crispnessidx(Y)
% get_crispnessidx: calculate crispness index as in
%   (http://dx.doi.org/10.1016/j.jneumeth.2017.07.031)
%
% Usage:
%   crispnessidx = get_crispnessidx(Y)
%
% Args:
%   Y: input image over time, last dimension is time
%       ([height width time] or [height width depth time]), no NaNs
%
% Returns:
%   crispnessidx: scalar (arbitrary fluorescence units), the norm of
%       the time-averaged image; larger = sharper mean image.
%       Note: it is the norm of the mean image itself, not of its gradient
%       (kept as is so values stay comparable with earlier runs).

crispnessidx = [];

siz_ = size(Y);
pre_matrix = abs(reshape(squeeze(mean(Y, length(siz_))),...
    [prod(siz_(1:end-1)) 1]));

crispnessidx = norm(pre_matrix, 'fro');

end