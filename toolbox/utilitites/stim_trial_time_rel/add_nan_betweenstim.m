function [Yn, time, stim_trial_nn, stim_on_nn] = ...
    add_nan_betweenstim(Y, stim_trial, stim_on, nan_n, time)
% add_nan_betweenstim: group the time points of Y by stimulus (trial) number and
%   append nan_n NaN columns after each group, so plots show a gap between them
%
% Usage:
%   [Yn, time, stim_trial_nn, stim_on_nn] = ...
%       add_nan_betweenstim(Y, stim_trial, stim_on, nan_n, time)
%
% Args:
%   Y: ROIs x T matrix
%   stim_trial: 1 x T vector, stimulus (trial) number of each time point
%   stim_on: 1 x T vector, stimulus-on marker of each time point
%   nan_n: number of NaN columns to append after each stimulus group
%   time: optional 1 x T time vector (seconds, uniform spacing)
%
% Returns:
%   Yn: ROIs x T2 matrix, time points regrouped by stimulus with NaN gaps
%   time: 1 x T2 time vector rebuilt from the first sample and dt (as given if empty)
%   stim_trial_nn, stim_on_nn: 1 x T2 vectors regrouped like Yn, 0 in the gaps

if ~exist('time', 'var'); time = []; end % so the output is always defined

stim_n = unique(stim_trial); 
Yn = [];

% groups sit on the cell diagonal (empty off-diagonals), so cell2mat below
% joins them side by side in stimulus order
for i = 1:numel(stim_n)
    
    Yn{1, i} = Y(:, stim_trial == stim_n(i));
    Yn{1, i}(:, end + 1:end + nan_n) = nan;
    
    stim_trial_nn{i, i} = stim_trial(stim_trial == stim_n(i));
    stim_trial_nn{i, i}(end + 1:end + nan_n) = 0;
    
    stim_on_nn{i, i} = stim_on(stim_trial == stim_n(i));
    stim_on_nn{i, i}(end + 1:end + nan_n) = 0;
    
end

Yn = cell2mat(Yn);
stim_trial_nn = cell2mat(stim_trial_nn);
stim_on_nn = cell2mat(stim_on_nn);

if exist('time', 'var') && ~isempty(time)
    dt = time(2)-time(1);
    % build from integer steps to avoid the length drift of a 0:dt:end range
    time = time(1) + (0:numel(stim_trial_nn) - 1)*dt;
end

end
