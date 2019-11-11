function [CI, pVal] = corrpertrial(iTrace, sTrace, tInitEnd, btn)
% get CI per trial, assumes that iTrace and sTrace are XxT
if ~exist('btn', 'var') || isempty(btn); btn = 10^4; end % bootstrap n
iTrace = zscorebigmem(iTrace);
sTrace = zscorebigmem(sTrace);
pVal = []; CI = [];
for roi_i = 1:size(iTrace, 1)
    for trial_i = 1:size(tInitEnd, 1)
        chop_trace = iTrace(roi_i, tInitEnd(trial_i, 1):tInitEnd(trial_i, 2));
        chop_stim = sTrace(roi_i, tInitEnd(trial_i, 1):tInitEnd(trial_i, 2));
        [pVal(roi_i, trial_i), CI(roi_i, trial_i)] =  getCI_bootstrap(chop_trace, chop_stim, btn);
        clear chop_trace chop_stim
    end
end
end