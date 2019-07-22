function [stimuli_onset_offset, stimuli_cat_trace] = ...
    songstim_init_end_rDat(rDat, start_end, stimCha, findStim, minwidth)
% songstim_init_end_rDat: function reads rDat from song stim and outputs 
% the start and end of a song stimuli.
%
% Usage:
%   [stimuli_onset_offset, stimuli_cat_trace] = songstim_init_end_rDat(rDat, init_end, stimCha)
%
% Args:
%   rDat: prv metadata structure variable
%   start_end: start and end of trace
%   stimCha: analog output channel
%   findStim: gate to use stimuli trace to find stimuli start and end
%   minwidth: minimun width of stimuli (ms)
%
% Returns:
%   stimuli_onset_offset: stimuli onset and offset [onset, offset]
%   stimuli_cat_trace: vector of stimuli delivered througout experiment
%
% Note:
% rDat.log.silencePre and rDat.log.silencePost are in ms

if ~exist('stimCha', 'var') || isempty(stimCha); stimCha = 1; end
if ~exist('findStim', 'var') || isempty(findStim); findStim = 0; end
if ~exist('minwidth', 'var') || isempty(minwidth); minwidth = 100; end

% Raw input auditory stim
stimuli_cat_trace = [];

% get actual number of trials
sti_played = min([numel(rDat.log.silencePre), ...
    numel(rDat.stiStartSample), numel(rDat.log.silencePost)]);

for j = 1:sti_played
    
    stimuli_onset_offset(j, 1) = rDat.stiStartSample(j) + ...
        rDat.log.silencePre(j)*10 + 1;
    stimuli_onset_offset(j, 2) = rDat.stiStartSample(j + 1) - ...
        rDat.log.silencePost(j)*10;    
    stimuli_cat_trace = cat(1, stimuli_cat_trace, ...
        rDat.stimAll{rDat.sti(j)}(:, stimCha));
    % need to check if I require to reorder them when the order is not the
    % same as in the txt file (when using randomize option)
    
end

% Chop stimuli_cat_trace and stimuli_onset_offset
stimuli_cat_trace = stimuli_cat_trace(start_end(1):start_end(2))';

if findStim
    stimuli_onset_offset = find_stim_int(stimuli_cat_trace, minwidth);
else
    stimuli_onset_offset = stimuli_onset_offset - start_end(1); 
end

l2chop = find(stimuli_onset_offset(:, 1) > start_end(2));
if ~isempty(l2chop)
    stimuli_onset_offset = stimuli_onset_offset(1:(l2chop - 1), :);
end

end

function stimuli_onset_offset = find_stim_int(itrace, minwidth)
% find_stim_int: function reads stimuli trace and outputs 
% the start and end of stimuli.
%
% Usage:
%   stimuli_onset_offset = find_stim_int(itrace, minwidth)
%
% Args:
%   itrace: input stimuli trace
%   minwidth: minimun width of stimuli (ms)

itrace = abs(itrace) > 0;

stimuli_onset_offset = [(find(diff(itrace) == 1) + 1)', find(diff(itrace) == -1)' + 1];

if size(stimuli_onset_offset, 1) > 1
    
    i = 1;
    k = 0;
    
    while k == 0
        
        if (stimuli_onset_offset(i + 1, 1) - stimuli_onset_offset(i, 2) < minwidth) ...
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