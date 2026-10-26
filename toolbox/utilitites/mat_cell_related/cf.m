function out = cf(func, varargin)
% cf: Convenience wrapper for cellfun with non-uniform output.
%
% Usage:
%   out = cf(func, C);
%
% Args:
%   func: function handle
%   varargin: cell array(s) passed to cellfun
%
% Returns:
%   out: cell array of the (non-uniform) outputs
%
% Notes:
%   equivalent to out = cellfun(func, C, 'UniformOutput', false)

out = cellfun(func, varargin{:}, 'UniformOutput', false);

end