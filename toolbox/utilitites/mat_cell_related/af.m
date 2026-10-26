function out = af(func, varargin)
% af: Convenience wrapper for arrayfun with non-uniform output.
%
% Usage:
%   out = af(func, A);
%
% Args:
%   func: function handle
%   varargin: array(s) passed to arrayfun
%
% Returns:
%   out: cell array of the (non-uniform) outputs
%
% Notes:
%   equivalent to out = arrayfun(func, A, 'UniformOutput', false)

out = arrayfun(func, varargin{:}, 'UniformOutput', false);

end

