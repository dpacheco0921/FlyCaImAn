function idx2keep = threshold_C(C, A, YrA, iparams)
% idx2keep = threshold_C(C, A, ths, YrA, tIdx, iparams)
% function to do dF thresholding based on the max delta F induced by the LED bleedthrough.
% The idea is that the 90 - 10 percentile of the dF for each component
% should be greater than ths.

lopars.hperc = 90;
lopars.lperc = 10;
lopars.t_init = 7.5; % timepoints or 7.5 sec
lopars.ths_ratio = 0.9;
lopars.dtype = 1; % 0 = lowpass filter, 1 = percentile filter
lopars.sr = 2;
lopars.lfreq = 0.002;
lopars.d_prct = 20;
lopars.d_window = 100;
lopars.d_shift = 10;
lopars.l_df = 0;
idx2keep = [];

if ~exist('iparams', 'var'); iparams = []; end
lopars = loparam_updater(lopars, iparams);

if ~exist('YrA', 'var') || ...
    isempty(YrA); YrA = [];
end

if ~isempty(C)
    
    if ~isempty(YrA)
        C = C + YrA;
    end
    
    % Detrend C
    if lopars.dtype == 0
        [filt_b, filt_a] = butter(2, ...
            lopars.lfreq /(lopars.sr/2), 'low');
        C = C - filtfilt(filt_b, filt_a, C')';
    else
        C = prctfilt(C, lopars.d_prct, ...
            lopars.d_window, lopars.d_shift);
    end
    
    % Only use timepoints after the buffer time:
    tIdx = round(lopars.t_init*lopars.sr);
    
    % Get lower and upper precentile
    C_lp = prctile(C(:, tIdx:end), lopars.lperc, 2);
    C_hp = prctile(C(:, tIdx:end), lopars.hperc, 2);
    
    % Threshold raw trace
    Cdf_ths = lopars.l_df ./ full(max(A, [], 1))';
    idx2keep = find((C_hp - C_lp) ./ Cdf_ths > lopars.ths_ratio);
    
end

end
