function [Y_norm] = normsdata(Y, t_span)
% normdata: normalize data to unity
%   and smooth the square value
% 
% Usage:
%   Y_norm = normdata(Y, p, dim)
%
% Args:
%   Y: 2DxT matrix
%   t_span: time span for smoothing
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
