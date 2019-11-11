function [stim_idx, bs_idx] = get_stimrel_timeidx(wDat, bs_range)
% get_stimrel_timeidx: generates roi density
% Usage:
%   get_stimrel_timeidx(wDat, tRange)
%
% Args:
%   wDat: wDat structure with fTime and sTime fields
%   bs_range: indeces of files to group
%
% Notes: saves data in obj.temp.roidensity

if ~exist('bs_range', 'var'); bs_range = [-20 -1]; end

% Round time (assumes total duration of stimuli is close to an integer)
wDat.fTime = round(wDat.fTime*10)/10; 
wDat.sTime = round(wDat.sTime*10)/10;

% Get time index, parsing trials per stim type
stim_idx = []; bs_idx = [];

for i = 1:size(wDat.sTime, 1)
    
    t_i = find(wDat.fTime == wDat.sTime(i, 1));
    t_e = find(wDat.fTime == wDat.sTime(i, 2));
    stim_idx{i, 1} = t_i:t_e;
    bs_idx{i, 1} = t_i + bs_range(1):t_i + bs_range(2);
    clear t_i t_e
    
end

end