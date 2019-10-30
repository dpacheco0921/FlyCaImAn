function out = cf(func, varargin)
% cf: Convenience wrapper for cellfun with non-uniform output.
%
% Usage:
%   out = cf(func, C);
%
% Args:
%   func: function
%   varargin: variables
%
% Notes:
%   equivalent to out = cellfun(func, C, 'UniformOutput, false)

out = cellfun(func, varargin{:}, 'UniformOutput', false);

end