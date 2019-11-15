function idx2keep = threshold_C(C, A, YrA, iparams)
% threshold_C: function to do dF thresholding based
%   on the max delta F induced by the LED bleedthrough.
%   The idea is that the 90 - 10 percentile of the dF for each component
%   should be greater than ths.
%
% Usage:
%   idx2keep = threshold_C(C, A, YrA, iparams)
%
% Args:
%   C: predicted temporal response
%   A: spatial component
%   YrA: residual per component
%   iparams: internal parameters
%       (hperc: upper percentile)
%           (default, 90)
%       (lperc: lower percentile)
%           (default, 10)
%       (t_init: time to use as start, in secs)
%           (default, 7.5)
%       (ths_ratio: threshold ratio)
%           (default, 0.9)
%       (sr: sampling rate (Hz))
%           (default, 2)
%       (d_prct: percentile)
%           (default, 20)
%       (d_window: window over which to compute the percentile)
%           (default, 100)
%       (d_shift: length of window shifting)
%           (default, 10)
%       (l_df: threshold for delta in intensity)
%           (default, 0)

% interal default parameters
% high percentile 
lopars.hperc = 90;
lopars.lperc = 10;
lopars.t_init = 7.5;
lopars.ths_ratio = 0.9;
lopars.sr = 2;
lopars.d_prct = 20;
lopars.d_window = 100;
lopars.d_shift = 10;
lopars.l_df = 0;
idx2keep = [];

if ~exist('iparams', 'var'); iparams = []; end
lopars = loparam_updater(lopars, iparams);

if ~exist('YrA', 'var') || ...
    isempty(YrA)
    YrA = [];
end

if ~isempty(C)
    
    if ~isempty(YrA)
        C = C + YrA;
    end
    
    % detrend C
    C = prctfilt(C, lopars.d_prct, ...
        lopars.d_window, lopars.d_shift);
    
    % only use timepoints after the buffer time:
    tIdx = round(lopars.t_init*lopars.sr);
    
    % get lower and upper precentile
    C_lp = prctile(C(:, tIdx:end), lopars.lperc, 2);
    C_hp = prctile(C(:, tIdx:end), lopars.hperc, 2);
    
    % threshold raw trace
    Cdf_ths = lopars.l_df ./ full(max(A, [], 1))';
    idx2keep = find((C_hp - C_lp) ./ Cdf_ths > lopars.ths_ratio);
    
end

end
