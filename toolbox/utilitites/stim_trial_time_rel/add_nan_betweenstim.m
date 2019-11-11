function [Yn, time, stim_trial_nn, stim_on_nn] = ...
    add_nan_betweenstim(Y, stim_trial, stim_on, nan_n, time)

stim_n = unique(stim_trial); 
Yn = [];

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
    time = (0:dt:(numel(stim_trial_nn) - 1)*dt) + time(1);
end

end
