function stimuli_onset_offset = ...
    find_stim_int(itrace, minwidth, threshold_val)
% find_stim_int: function reads stimuli trace and outputs 
%   the start and end of stimuli.
%
% Usage:
%   stimuli_onset_offset = ...
%       find_stim_int(itrace, minwidth, threshold_val)
%
% Args:
%   itrace: input stimuli trace (vector, volts, one value per sample)
%   minwidth: minimun gap between two stimuli (samples), closer
%       stimuli are merged into one
%       (default, 0, nothing is merged)
%   threshold_val: voltage threshold, |trace| above it counts as stimulus
%       (default, 0)
%
% Returns:
%   stimuli_onset_offset: [n-stimuli x 2] sample index of the first
%       sample above threshold (onset) and of the first sample back
%       below threshold (offset)
%
% Notes:
% onsets and offsets are paired by order, so a trace that starts
%   already above threshold will pair the wrong edges.

% an empty minwidth would make the comparison below fail (and an empty
%   threshold would select nothing), so fall back to the defaults
if ~exist('minwidth', 'var') || isempty(minwidth)
    minwidth = 0;
end

if ~exist('threshold_val', 'var') || isempty(threshold_val)
   threshold_val = 0;
end

itrace = abs(itrace) > threshold_val;

if size(itrace, 1) > 1
    itrace = itrace';
end

% get onset
onset = (find(diff(itrace) == 1) + 1)';
% get offset
offset = (find(diff(itrace) == -1) + 1)';

min_num = min([numel(onset) numel(offset)]);

stimuli_onset_offset = [onset(1:min_num), offset(1:min_num)];

if size(stimuli_onset_offset, 1) > 1
    
    i = 1;
    k = 0;
    
    while k == 0
        
        if (stimuli_onset_offset(i + 1, 1) - ...
                stimuli_onset_offset(i, 2) < minwidth) ...
                && (i < size(stimuli_onset_offset, 1))
            
            % update upper limit
            stimuli_onset_offset(i, 2) = stimuli_onset_offset(i + 1, 2);
            % delete contiguous one
            stimuli_onset_offset(i + 1, :) = []; 
            
        else
            
            i = i + 1; 
            
        end
        
        k = double(i == size(stimuli_onset_offset, 1));
        
    end
    
end

end
