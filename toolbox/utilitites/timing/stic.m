function timer_id = stic
% stic: tic equivalent using precise system time.
%
% Usage:
%   stic
%   timer_id = stic
%
% Returns (only if requested):
%   timer_id: integer id of this timer, to pass to stoc
%
% See also: stoc, systime

global stic_timers;

timer_id = numel(stic_timers) + 1;
stic_timers(timer_id) = systime;

if nargout < 1
    clear timer_id
end

end

