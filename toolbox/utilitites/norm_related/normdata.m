function Y_norm = normdata(Y, p, dim)
% normdata: normalize Y to unity 
%   (generalized norm function to work with matrices of group of vectors)
% 
% Usage:
%   Y_norm = normdata(Y, p, dim)
%
% Args:
%   Y: 2DxT matrix
%   p: norm type, see norm (default, 2)
%   dim: dimension across which to normalize (default, 1)
%
% Returns:
%   Y_norm: same size as Y (double); each vector along dim is divided by its p-norm
%       (all-zero vectors give NaN)
%
% Notes: see norm

if ~exist('p', 'var') || isempty(p)
    p = 2;
end

if ~exist('dim', 'var') || isempty(dim)
    dim = 1;
end

% normalize data to unity

% write norm function (euclidean norm of vector)
fnorm = eval(['@(x) x/norm(x, ', num2str(p), ')']);

% edit Y to run on a cellfun
if dim ~= 1
    Y = Y';
end

Y = double(Y);
Y_cell = chunk2cell(Y, 1, 1);
Y_norm = cellfun(fnorm, Y_cell, 'UniformOutput', false);
Y_norm = cell2mat(Y_norm);

% undo edit
if dim ~= 1
    Y_norm = Y_norm';
end

end
