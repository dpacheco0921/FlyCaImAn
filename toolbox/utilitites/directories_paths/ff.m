function varargout = ff(varargin)
% ff Description: get full file of each varargin{i}
% Usage:
%   varargout = ff(varargin)
% 
% Args:
%   varargin: path parts to join (strings), as in fullfile
%
% Returns:
%   varargout: the full file path (a single output, since fullfile
%       returns one string)
%
% See also:

N = max(nargout,1);
varargout{1:N} = fullfile(varargin{:});

end
