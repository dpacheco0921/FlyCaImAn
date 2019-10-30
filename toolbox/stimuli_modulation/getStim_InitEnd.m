function stEn = getStim_InitEnd(...
    frame_time, stimuli_time)
% getStim_InitEnd: generates a [n, 2] matrix with 
%   stimuli start and end in frame indeces
%
% Usage:
%   stEn = getStim_InitEnd(...
%      frame_time, stimuli_time)
%
% Args:
%   frame_time: frame timestamps (seconds)
%   stimuli_time: stimuli start and end (seconds)

% round to tens of miliseconds
frame_time = round(frame_time*100)/100;
stimuli_time = round(stimuli_time*100)/100;

dt = frame_time(2) - frame_time(1);

for trial_i = 1:size(stimuli_time, 1)
    
    % dt stim
    dt_stim = stimuli_time(trial_i, 2) - stimuli_time(trial_i, 1);
    
    % compute difference to get the timepoint
    %   closer to stimuli init or end
    dt_i = abs(frame_time - stimuli_time(trial_i, 1));
    dt_e = abs(frame_time - stimuli_time(trial_i, 2));
    
    % maybe systematically choose the timepoint
    %   before rather than after stim star?
    stEn(trial_i, 1) = find(dt_i == min(dt_i), 1);
    stEn(trial_i, 2) = find(dt_e == min(dt_e), 1);
    
    % correct length
    d_stEn = stEn(trial_i, 2) - stEn(trial_i, 1);
    
    if d_stEn > round(dt_stim/dt)
       stEn(trial_i, 2) = stEn(trial_i, 2) - 1;
    elseif d_stEn < round(dt_stim/dt) && d_stEn < numel(frame_time)
       stEn(trial_i, 2) = stEn(trial_i, 2) + 1;
    end
    
    clear dt_i dt_e d_stEn
    
end

end
