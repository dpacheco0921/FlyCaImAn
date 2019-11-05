function [stimvect, stEn, stimuli_name] = ...
    getStimVect(wDat, delta, trial2use, stim2use)
% getStimVect: generates a binary vector of the stimuli in 
%   imaging time resolution
%
% Usage:
%   [stimvect, stEn, stimuli_name] = ...
%       getStimVect(wDat, delta, trial2use, stim2use)
%
% Args:
%   wDat: input main metadata variable
%   delta: delta to use for binary vector
%   trial2use: specific trials to use 
%   stim2use: specific stimuli (by ordered idx)
%
% Notes:
% 1) Uses wDat.sTime, wDat.fTime, wDat.sPars.int, 
%   wDat.sPars.name, wDat.sPars.order
%   to generate the binary vector of the stimuli
% 2) Rounds time to hundreds of miliseconds.
% 3) get stimuli names and indexes (assumes stim is the same (number and
%   order) across subvolumes).
% 4) any stimuli shorter than dt will be assigned a dt width.

if ~exist('trial2use', 'var'); trial2use = []; end
if ~exist('stim2use', 'var'); stim2use = []; end

% round to tens of miliseconds
wDat.fTime = round(wDat.fTime*100)/100;
wDat.sTime = round(wDat.sTime*100)/100;

% remove stimuli after end of recording
stim2del = wDat.sTime(:, 1) > wDat.fTime(end) | ...
    wDat.sTime(:, 2) > wDat.fTime(end);
wDat.sTime = wDat.sTime(~stim2del, :);

% update name using stimuli name plus intensity and stimuli idx
stimuli_name = cellfun(@(x, y) [x '_AMP_' num2str(y)], ...
    wDat.sPars.name(1, :), chunk2cell(wDat.sPars.int(1, :), 1), ...
    'UniformOutput', false);

[~, ~, stim_u_idx] = unique(stimuli_name, 'stable');
stim_all_idx = stim_u_idx(wDat.sPars.order(1, :));

if size(stim_all_idx, 1) == 1
    stim_all_idx = stim_all_idx';
end

stim_all_idx = stim_all_idx(1:size(wDat.sTime, 1), 1);

if ~exist('delta', 'var') || isempty(delta)
    delta = 1;
end

% select stim and trials to use:
if ~isempty(trial2use) && isempty(stim2use)
    
    % if only trials2use is provided consider all trials across stimuli
    wDat.sTime = wDat.sTime(trial2use, :);
    
elseif ~isempty(stim2use) && isempty(trial2use)
    
    % if only stim2use consider all trials per selected stimuli
    if size(stim2use, 1) > 1
        stim2use = stim2use';
    end
    
    % collect all trials that correspond to selected stimuli
    for i = stim2use
        trial2use = [find(stim_all_idx == i)', trial2use];
    end
    
    trial2use = sort(trial2use);
    wDat.sTime = wDat.sTime(trial2use, :);
    
elseif ~isempty(stim2use) && ~isempty(trial2use)
    
    % if stim2use & trial2use consider selected trials per selected stimuli
    if size(stim2use, 1) > 1
        stim2use = stim2use';
    end
    
    trial2use_across_stim = [];
    
    % collect all trials that correspond to selected stimuli
    for i = stim2use
        trialsperstim = find(stim_all_idx == i)';
        trial2use_across_stim = ...
            [trialsperstim(trial2use), ...
            trial2use_across_stim];
    end
    
    trial2use_across_stim = sort(trial2use_across_stim);
    wDat.sTime = wDat.sTime(trial2use_across_stim, :);    
    
    % replace trial2use (relative to stimuli)
    %   to abs trials to use (across stimuli)
    trial2use = trial2use_across_stim;
    
else
    
    % include all
    trial2use = 1:size(wDat.sTime, 1);
    
end

maxtrials = numel(trial2use);

% build stimvect
stimvect = zeros(1, numel(wDat.fTime));

if ~isempty(wDat.sTime)
    
    % get stimuli start and end in frame indeces
    stEn = getStim_InitEnd(wDat.fTime, wDat.sTime);

    for trial_i = 1:maxtrials
        stimvect(1, stEn(trial_i, 1):stEn(trial_i, 2)) = delta;
    end

else
    
    fprintf('No trials with this stimulus found\n')
    stEn = [];
    
end

end
