function [dt, timer_id] = stocf(timer_id, str, varargin)
% stocf: Report elapsed time with print formatting.
%
% Usage:
%   stocf(str, ...)
%   stocf(timer_id, str, ...)
%   [dt, timer_id] = stocf(_)
%
% Args:
%   timer_id: id returned by stic (optional; default, the most recent timer)
%   str: sprintf-style message; elapsed seconds are appended as [x.xxs]
%       (or h/m/s text after 5 min) unless an empty [] is passed in varargin,
%       in which case dt replaces the last [] there
%   varargin: values for the sprintf-style message
%
% Returns (only if requested):
%   dt: elapsed seconds since the timer was started
%   timer_id: id of the timer that was read
%
% See also: stic, stoc

% Get elapsed time
if ischar(timer_id)
    
    if nargin > 1
        varargin = [{str} varargin];
    end
    
    % No timer specified, so just use last
    str = timer_id;
    [dt, timer_id] = stoc;
    
else
    dt = stoc(timer_id);
end

% Figure out where to plug in the elapsed time
dt_idx = find(areempty(varargin),1,'last');

% Use default formatting if no empty arrays specified
if isempty(dt_idx)
    
    dt_idx = numel(varargin) + 1;
    
    % pad with space
    if str(end) ~= ' '
        str(end+1) = ' ';
    end
    
    if dt > 5*60
        
        % after 5 mins, use string representation
        str = [str '[' secsf(dt) ']'];
        
    else
        
        % Put elapsed time in argument list
        str = [str '[%.2fs]'];
        varargin{dt_idx} = dt;
        
    end
    
else
    
    % an empty [] was passed in varargin: plug the elapsed time in there
    varargin{dt_idx} = dt;
    
end

% Print!
printf(str, varargin{:})

% No returns if no output requested
if nargout == 0
    clear dt timer_id
end

end
