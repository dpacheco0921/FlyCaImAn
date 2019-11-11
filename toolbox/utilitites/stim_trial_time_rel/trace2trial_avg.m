function [wDat, CaRaw_zs_tl, CaPred_zs_tl, CaRef_zs_tl, ...
    stim_idx, stim_vect, time_tl, time_tl_rel] = ...
    trace2trial_avg(tRmode, stim2use, time_range, ...
    CaRaw, CaPred, CaRef, wDat)
% st_pros_roi: split trace into trials relative to stimuli type
%
% Usage:
%   [wDat, CaRaw_zs_tl, CaPred_zs_tl, CaRef_zs_tl, ...
%       stim_idx, stim_vect, time_tl, time_tl_rel] = ...
%       trace2trial_avg(tRmode, stim2use, time_range, ...
%       CaRaw, CaPred, CaRef, wDat)
%
% Args:
%   tRmode: use time_range to be relative to 
%       stimuli start time_range(1) and stimuli end time_range(2)
%   stim2use: stimuli to use
%   time_range: range of time to use relative to stimulus onset [a b].
%       (for example [-5 10] will load chunks from -5 seconds to 
%       plus 10 seconds relative to stimumus onset)
%   CaRaw, CaPred, CaRef: raw data, predicted data, and reference data
%   wDat
% 
% Notes: 
% 1) Chop trace into trials [init end] defined by "time_range" (in seconds)
%   tRmode defines is it used relative to stimuli start or to both start
%   and end
% 2) Compile them into trials for each stimuli used
% 3) do this for raw and predicted data

if ~exist('tRmode', 'var')
    tRmode = 0;
end

if ~exist('stim2use', 'var')
    stim2use = [];
end

if ~exist('time_range', 'var')
    time_range = [];
end

if isempty(time_range)
    fprintf('Error needs: time_range');
    return;
end

% collect stimuli info
stimuli_name = cellfun(@(x, y) [x '_AMP_' num2str(y)], ...
wDat.sPars.name(1, :), chunk2cell(wDat.sPars.int(1, :), 1), ...
'UniformOutput', false);

% assumes that stim order is the same across all segments
[wDat.sName, ~, stim_u_idx] = ...
    unique(stimuli_name, 'stable');
stim_all_idx = stim_u_idx(wDat.sPars.order(1, :));

if size(stim_all_idx, 1) == 1
    stim_all_idx = stim_all_idx';
end

wDat.sIdx = stim_all_idx(1:size(wDat.sTime, 1), 1);

clear stim_all_idx stim_u_idx stimuli_name
  
% zscore signal
CaRaw_zs = zscorebigmem(CaRaw);

[~, a3, ~, ~, a1, ~, ~, a4, a5, a6] = ...
    trace2trials(wDat, CaRaw_zs, ...
    time_range, wDat.sIdx, tRmode, stim2use);

% store splitted trials:
% raw zs-signal per ROI (all trials) (stitched along stim type)
CaRaw_zs_tl = a1;

% store stim structure
% time indeces of trial per stimuli type
stim_idx = a5;
% time indeces of stimuli occurrence (coded by stimuli type)
stim_vect = a6;

% store time info
% absolute increase in time
time_tl = a3;
% increase in time relative to each stim (reset per stim)
time_tl_rel = a4;

[~, ~, ~, ~, a2] = trace2trials(wDat, ...
    CaPred, time_range, wDat.sIdx, tRmode, stim2use);

Ca_temp = zscorebigmem(CaRef);
[~, ~, ~, ~, a2_Ref] = trace2trials(wDat, Ca_temp, ...
    time_range, wDat.sIdx, tRmode, stim2use);

% predicted signal per ROI (stitched along stim type)
CaPred_zs_tl = a2;

% signal from reference channel (red channel/structural channel)
CaRef_zs_tl = a2_Ref;
    
end
