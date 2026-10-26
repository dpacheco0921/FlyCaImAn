function [Y_norm] = normsdata(Y, t_span)
% normsdata: normalize each row to unit norm
%   and smooth the squared values
% 
% Usage:
%   Y_norm = normsdata(Y, t_span)
%
% Args:
%   Y: 2DxT matrix
%   t_span: time span (number of samples) for smoothing
%
% Returns:
%   Y_norm: 2DxT matrix (double); each row is x/norm(x), squared, then smoothed
%
% Notes: see norm

% normalize data to unity and smooth the square value
fnorm = @(x) x/norm(x);

Y_norm = double(Y);

nSamp = size(Y_norm, 1);

for j = 1:nSamp
    Y_norm(j, :) = smooth(fnorm(Y_norm(j, :)).^2, t_span);
end

end
