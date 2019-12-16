function [stimuli_onset_offset, stimuli_cat_trace] = ...
    songstim_init_end(rDat, start_end, stimCha, minSize)
% songstim_init_end: function reads rDat from song stim and outputs 
% the start and end of a song stimuli.
%
% Usage:
%   [segY, stimuli_cat_trace] = ...
%       songstim_init_end(rDat, init_end, stimCha, minSize)
%
% Args:
%   rDat: prv metadata structure variable
%   start_end: start and end of trace
%   stimCha: analog output channel
%   minSize: is the min inter pulse interval to consider two pulse as a separate bout
%       (default, 10000 or 1000ms)
%
% Returns:
%   stimuli_onset_offset: stimuli onset and offset [onset, offset]
%
% Notes:

if ~exist('stimCha', 'var') || ...
        isempty(stimCha)
    stimCha = 1;
end

if ~exist('minSize', 'var') || ...
        isempty(minSize)
    minSize = 10000;
end

% raw input auditory stim
stimuli_cat_trace = [];

for j = rDat.sti
    stimuli_cat_trace = ...
        cat(1, stimuli_cat_trace, rDat.stimAll{j}(:, stimCha));
    % need to check if I require to reorder them when the order is not the
    % same as in the txt file (when using randomize option)
end

stimuli_cat_trace = ...
    stimuli_cat_trace(start_end(1):start_end(2));

% binarize the trace and find bout init and end
inputTrace_bi = abs(stimuli_cat_trace) > 0;
stimuli_onset_offset = ...
    [(find(diff(inputTrace_bi) == 1) + 1), ...
    find(diff(inputTrace_bi) == -1) + 1];

% plot hist
%figure()
%diffsegY = segY(2:end,1) - segY(1:end-1,1);
%hist(diffsegY, min(diffsegY):1:max(diffsegY))

% concatenate small pulses
if size(stimuli_onset_offset, 1) > 1
    
    i = 1;
    k = 0;
    
    while k == 0
        
        if (stimuli_onset_offset(i + 1, 1) - stimuli_onset_offset(i, 2) < minSize) && ...
                (i < size(stimuli_onset_offset, 1))
            
            % update upper limit
            stimuli_onset_offset(i, 2) = ...
                stimuli_onset_offset(i + 1, 2);
            % delete contiguous one
            stimuli_onset_offset(i + 1, :) = [];
            
        else
            
            i = i + 1;
            
        end
        
        k = double(i == size(stimuli_onset_offset, 1));
        
    end
end

end
