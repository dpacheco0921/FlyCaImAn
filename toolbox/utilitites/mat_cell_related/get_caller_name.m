function caller = get_caller_name(varargin)
%GET_CALLER_NAME Returns the name of the caller function.
% Usage:
%   caller = get_caller_name()
%   caller = get_caller_name('Path', true)
%   caller = get_caller_name(..., 'Warn', true) % warning if called from workspace
%
% Args (name-value):
%   Path: if true, return the file path instead of the function name (default, false)
%   Warn: warn when called from the workspace (default, false)
%
% Returns:
%   caller: name (or path) of the function that CALLED the function calling
%       get_caller_name (two levels up the stack; used by parse_params to name the
%       function being parsed); empty if there is no such caller
%
% See also: mfilename, funpath, dbstack

% Parse inputs
defaults.Path = false;
defaults.Warn = false;
params = defaults;
if nargin > 0
    params = parse_params(varargin, defaults);
end

% Get function call stack
if params.Path
    [ST, I] = dbstack('-completenames');
else
    [ST, I] = dbstack();
end

% This function or the function that called it were called directly from
% the workspace
if numel(ST) < I + 2
    caller = '';
    if params.Warn
        if numel(ST) < I + 1
            warning('%s called directly from workspace.', mfilename)
        else
            warning('%s called directly from workspace.', ST(I + 1).name)
        end
    end
    return
end

% Return caller
if params.Path
    caller = ST(I + 2).file;
else
    caller = ST(I + 2).name;
end


end
