function [h, m, s] = secs2hms(numSecs)
% secs2hms: Converts a number of seconds to hours, minutes and fractional seconds.
%
% Usage:
%   hms = secs2hms(numSecs) % numeric vector
%   [h, m, s] = secs2hms(numSecs)
%
% Args:
%   numSecs: number of seconds (scalar)
%
% Returns:
%   h, m: whole hours and minutes
%   s: remaining seconds, including the fractional part
%
% See also: duration, hms, secs2str

[h, m, s] = hms(duration([0 0 numSecs]));

end
