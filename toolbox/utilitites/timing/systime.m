function secs = systime
% systime: Returns precise system time in seconds.
%
% Usage:
%   secs = systime
%
% Returns:
%   secs: time in seconds (double). Only differences between two calls are
%       meaningful: on non-Windows systems the constant below also encodes a
%       hardcoded 4 h local-time offset, so absolute values differ from the
%       Windows call.
%
% Note: This is just a wrapper for GetSystemTimePreciseAsFileTime
%
% See also: GetSystemTimePreciseAsFileTime

if ispc
    secs = GetSystemTimePreciseAsFileTime;
else
    secs = now * 86400 - 50522817600; % same units
end

end

