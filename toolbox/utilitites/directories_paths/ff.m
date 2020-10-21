function varargout = ff(varargin)
% ff Description: get full file of each varargin{i}
% Usage:
%   varargout = ff(varargin)
% 
% Args:
%   varargin: files to run
% 
% See also:

N = max(nargout,1);
varargout{1:N} = fullfile(varargin{:});

end
