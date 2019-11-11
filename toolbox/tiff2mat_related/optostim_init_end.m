function stimuli_onset_offset = ...
    optostim_init_end(trace, lStim)
% optostim_init_end: Using the stim trace, 
%   it gets the beginning and end of stimulation(LED)
%
% Usage:
%   stimuli_onset_offset = ...
%       optostim_init_end(trace, lStim)
%
% Args:
%   trace: stimuli trace
%   lStim: stimuli variable from LEDcontroler (only requires .freq and .fs)
%
% Returns:
%   stimuli_onset_offset: stimuli onset and offset [onset, offset]
%
% Notes:
% assumes any change from 0 is an stimulus

% Using the stim trace, it gets the beginning and end of stimulation (LED)
pulseIdx = diff(trace);  
pulsenInit = find(pulseIdx > 0) + 1;
pulseEnd = find(pulseIdx < 0);

% find trains
try
    % in time units (a.u.)
    Ths = (1/lStim.freq(1))*lStim.fs(1);
catch
    % in time units (a.u.)
    Ths = (1/lStim.freq(1))*lStim.rate(1);
end

initDis = diff(pulsenInit);
peaks = find(initDis > Ths*2);
stimuli_onset_offset = zeros(numel(peaks) + 1, 2);
stimuli_onset_offset(1,1) = pulsenInit(1);
stimuli_onset_offset(end, end) = pulseEnd(end);

for i = 2:size(stimuli_onset_offset, 1)
    
    stimuli_onset_offset(i, 1) = pulsenInit(peaks(i-1) + 1);
    stimuli_onset_offset(i-1, 2) = pulseEnd(peaks(i-1));
    
end

end
