function [trialvect, trial_init, trial_end] = ...
    getTrialVect(wDat, stim2merge, customtrial_init, ...
    customtrial_end, trials2use)
% getTrialVect: generates a binary vector of the stimuli
%
% Usage:
%   [trialvect, trial_init] = ...
%       getTrialVect(wDat, stim2merge, customtrial_init, ...
%      customtrial_end, trials2use)
%
% Args:
%   wDat: input main metadata variable
%   stim2merge: stimuli to merge for index vector
%   customtrial_init: number of seconds before stimuli 
%       start to count trial start (seconds) abs number
%   customtrial_end: number of seconds after stimuli
%       end to count trial end (seconds) abs number
%   trials2use: which trials to use
%
% Notes:
% Use wDat.sTime, wDat.fTime, wDat.sPars.basPre, 
%   wDat.sPars.order, and wDat.sPars.name, 
%   wDat.sPars.name, wDat.sPars.int
% to generate a trial idx vector

if ~exist('customtrial_init', 'var'); customtrial_init = []; end
if ~exist('customtrial_end', 'var'); customtrial_end = []; end
if ~exist('trials2use', 'var'); trials2use = []; end

% round to tens of miliseconds
wDat.fTime = round(wDat.fTime*100)/100;
wDat.sTime = round(wDat.sTime*100)/100;

% get start of a trial (basPre + stim-on + basPost)
trial_init = zeros(size(wDat.sTime, 1), 1);
trial_end = zeros(size(wDat.sTime, 1), 1);

% Get init of all stims & trials
for trial_i = 1:size(wDat.sTime, 1)
   
    % compute difference to get the timepoint closer to stimuli init or end
    if isempty(customtrial_init)
        dt_i = abs(wDat.fTime - wDat.sTime(trial_i, 1) ...
            + wDat.sPars.basPre(1, wDat.sPars.order(1, trial_i)));
    else
        dt_i = abs(wDat.fTime - wDat.sTime(trial_i, 1) + customtrial_init);
    end
    
    % compute difference to get the timepoint closer to stimuli init or end
    if ~isempty(customtrial_end)
        dt_e = abs(wDat.fTime - wDat.sTime(trial_i, 2) - customtrial_end);
    end
    
    % Get begining of trial
    if trial_i == 1
        trial_init(trial_i, 1) = max([find(dt_i == min(dt_i), 1) 1]);
    else
        trial_init(trial_i, 1) = find(dt_i == min(dt_i), 1);
    end
    
    % Get end of trial
    if ~isempty(customtrial_end)
        trial_end(trial_i, 1) = ...
            min([find(dt_e == min(dt_e), 1), numel(wDat.fTime)]);
    end
    
    clear dt_i dt_e
    
end

if isempty(customtrial_end)
    trial_end(:, 1) = [trial_init(2:end, 1) - 1; numel(wDat.fTime)];
end

% update name using stimuli name plus intensity and stimuli idx
stimuli_name = cellfun(@(x, y) [x '_AMP_' num2str(y)], ...
    wDat.sPars.name(1, :), chunk2cell(wDat.sPars.int(1, :), 1), ...
    'UniformOutput', false);

[~, ~, stim_u_idx] = unique(stimuli_name, 'stable');
stim_all_idx = stim_u_idx(wDat.sPars.order(1, :));
stim_all_idx = stim_all_idx(1:size(wDat.sTime, 1), 1);

% make trial vect per stim
trialvect = nan(numel(unique(stim_u_idx)), numel(wDat.fTime));

for i = sort(unique(stim_u_idx))'
    
    k = 1;
    
    for trial_i = find(stim_all_idx == i)'
        trialvect(i, trial_init(trial_i):trial_end(trial_i)) = k;
        k = k + 1;
    end
    
end

% Select stims to use
trialvect = nansum(trialvect(stim2merge, :), 1);
if ~isempty(trials2use)
    trialvect(~ismember(trialvect, trials2use)) = 0;
end

end
