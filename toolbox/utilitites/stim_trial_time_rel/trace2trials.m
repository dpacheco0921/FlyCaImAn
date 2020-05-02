function [oTraces, oTime, tInitEnd, lInitEnd, ...
    oMean, oMed, oSEM, oTimerel, oStimidx, ...
    oStimvect] = trace2trials(wDat, iTrace, ...
    tRange, sTypeIdx, tRmode, sType)
% trace2trials: splits iTrace into trials per stimuli
%
% Usage:
%   [oTraces, oTime, tInitEnd, lInitEnd, ...
%       oMean, oMed, oSEM, oTimerel, oStimidx, ...
%       oStimvect] = trace2trials(wDat, iTrace, ...
%       tRange, tRmode, sType)
%
% Args:
%   wDat: input main metadata variable
%   iTrace: signal matrix, size: KxT (K, roi number, T, time)
%   tRange: timepoints relative to stimuli start to include in trials (in seconds)
%   sTypeIdx: vector with stimuli indexes (in numbers)
%   tRmode: use tRange to be relative to stimuli start tRange(1) and stimuli end tRange(2)
%   sType: stimuli to use
%
% Outputs:
%   oTraces: all trials per K
%   oTime: absolute time used
%   tInitEnd: time indexces of start-end trial
%   lInitEnd: time indexces of start-end of stim
%   oMean, oMed, oSEM: main statistics of trials per K

if ~exist('tRange', 'var')
    tRange = [-min(wDat.sPars.basPre(:)), ...
        max(wDat.sPars.basPost(:))];
end

if ~exist('sTypeIdx', 'var')
    sTypeIdx = ones(1, size(wDat.sTime, 1));
end

if ~exist('tRmode', 'var')
    tRmode = 0;
end

if ~exist('sType', 'var')
    [sType, i_idx] = unique(sTypeIdx);
end

[K, T] = size(iTrace);

% round time to tens of miliseconds 
%   (assumes duration of stimuli is close to a integer)
wDat.fTime = round(wDat.fTime*100)/100;
wDat.sTime = round(wDat.sTime*100)/100;

% get info of stim: type of stim, and duration per stim
if size(sTypeIdx, 1) > 1
    sTypeIdx = sTypeIdx';
end

if isempty(sType)
    
    [sType, i_idx] = unique(sTypeIdx);
    ii_idx = (1:numel(i_idx))';
    
else
    
    [~, i_idx] = unique(sTypeIdx);
    [~, ii_idx, i_intsc] = ...
        intersect(sType, unique(sTypeIdx));
    i_idx = i_idx(i_intsc);
    clear i_intsc
    
end

sCount = sum(cell2mat(arrayfun(@(x) x == sType, ...
    sTypeIdx, 'UniformOutput', false)'));

if size(sType, 1) > 1
    sType = sType';
end

sWidth = nan(numel(sType), 2);
sWidth(ii_idx, :) = wDat.sTime(i_idx, :);
sWidth = bsxfun(@minus, sWidth(:, 2), sWidth(:, 1));

% Get time index, parsing trials per stim type
lInitEnd = nan(max(sCount), numel(sType)*2);
ti = nan(max(sCount), numel(sType));
tend = nan(max(sCount), numel(sType));

dt = wDat.fTime(2) - wDat.fTime(1);

ii = 1;

for st_i = sType
    
    i = 1;
    t_idx = find(sTypeIdx(:) == st_i)';
    
    for trial_i = t_idx
        
        % 1) get time before and after stim defined by tRange
        
        % compute difference to get the timepoint
        %   closer to stimuli init or end
        if size(tRange, 1) == 1
            dt_stim = tRange(1, 2) - tRange(1, 1);
            dt_i = abs(wDat.fTime - wDat.sTime(trial_i, 1)...
                - tRange(1, 1));
            dt_e = abs(wDat.fTime - wDat.sTime(trial_i, 1)...
                - tRange(1, 2));
        else
            dt_stim = tRange(st_i, 2) - tRange(st_i, 1);
            dt_i = abs(wDat.fTime - wDat.sTime(trial_i, 1)...
                - tRange(st_i, 1));
            dt_e = abs(wDat.fTime - wDat.sTime(trial_i, 1)...
                - tRange(st_i, 2));            
        end
        
        % maybe systematically choose the 
        %   timepoint before rather than after
        % stim star?
        ti(i, ii) = find(dt_i == min(dt_i), 1);
        tend(i, ii) = find(dt_e == min(dt_e), 1);
                
        if tRmode
            % add stimuli length
            stimtp = round(sWidth(ii, 1)/dt);
            tend(i, ii) = tend(i, ii) + stimtp;
            
            % correct length
            
            d_stEn = tend(i, ii) - ti(i, ii);

            if d_stEn > round((dt_stim + sWidth(ii, 1))/dt) ...
                    && tend(i, ii) < T
               tend(i, ii) = tend(i, ii) - 1;
            elseif d_stEn < round((dt_stim + sWidth(ii, 1))/dt) ...
                    && tend(i, ii) < T
               tend(i, ii) = tend(i, ii) + 1;
            elseif tend(i, ii) > T
               tend(i, ii) = T; 
            end
        
            clear stimtp
            
        else
            % correct length
            
            d_stEn = tend(i, ii) - ti(i, ii);

            if d_stEn > round(dt_stim/dt) && tend(i, ii) < T
               tend(i, ii) = tend(i, ii) - 1;
            elseif d_stEn < round(dt_stim/dt) && tend(i, ii) < T
               tend(i, ii) = tend(i, ii) + 1;
            end
            
        end
        
        % 1) get stimuli start and end
        stEn = getStim_InitEnd(wDat.fTime, wDat.sTime(trial_i, :));
        lInitEnd(i, 1 + 2*(ii-1)) = stEn(1, 1);
        lInitEnd(i, 2 + 2*(ii-1)) = stEn(1, 2);
        
        %lInitEnd(i, 1 + 2*(st_i-1)) = ...
        %   find(wDat.fTime == wDat.sTime(trial_i, 1));
        %lInitEnd(i, 2 + 2*(st_i-1)) = ...
        %   find(wDat.fTime == wDat.sTime(trial_i, 2));
        
        clear dt_i dt_e
        
        i = i + 1;
        
    end
    
    ii = ii + 1;
    
end

% get max oTrace length (count timepoints per stimuli trial)
maxlength = max(nansum(tend - ti + 1, 2), [], 1);

% get stimuli trial start index
initidx = [0 cumsum(nanmax(tend - ti + 1, [], 1))];

% allocate variables
oStimidx = ones(1, maxlength);
oMean = nan(K, maxlength);
oTime = nan(K, maxlength);

% compile all trials per stimuli type
% all stimuli are concatenated along columns, and trials are different rows
tInitEnd{1} = [];
tInitEnd{2} = initidx + 1;

for roi_i = 1:K
    
    % generate matrix trials x all-stimuli-type (concatenated)
    oTraces{roi_i, 1} = ...
        nan(size(tend, 1), maxlength);
    
    % find and allocate each trial per stimuli type to oTraces 
    for trial_j = 1:size(tend, 1)
        
        for stim_j = 1:size(tend, 2) % collect all indeces
            
            if ~isnan(ti(trial_j, stim_j)) ...
                    && ~isnan(tend(trial_j, stim_j))
                
                try
                    
                    ts_idx = ti(trial_j, stim_j)...
                        :tend(trial_j, stim_j);

                    if trial_j == 1
                       tInitEnd{1} = [tInitEnd{1}, ts_idx]; 
                    end

                    oTraces{roi_i, 1}(trial_j, initidx(stim_j) ...
                        + (1:numel(ts_idx))) = iTrace(roi_i, ts_idx);
                    clear ts_idx

                    if trial_j == 1
                        oStimidx(1, (initidx(stim_j) + 1)...
                            :initidx(stim_j + 1)) = sType(stim_j);
                    end
                    
                catch
                    
                    keyboard
                    
                end
            end
            
        end
    end
    
    % calculate basic stats of responses across trials (disregard nans)
    oMean(roi_i, :) = ...
        nanmean(oTraces{roi_i, 1}, 1);
    oMed(roi_i, :) = ...
        nanmedian(oTraces{roi_i, 1}, 1);
    oSEM(roi_i, :) = ...
        nanstd(oTraces{roi_i, 1}, [], 1)/numel(tend)^2;
    
    clear pre_oTrace
    
    if mod(roi_i, 300) == 0
        fprintf('%2.1f%% of ROIs completed \n', ...
            roi_i*100/K);
    end
    
end

% generate time-related variables
oTime = tRange(1, 1):dt:((maxlength - 1)*dt ...
    + tRange(1, 1));
oTimerel = zeros(1, maxlength);

% generate stimuli-related variables
oStimvect = zeros(1, maxlength);

for stim_ = 1:size(tend, 2)
    
    if ~isnan(initidx(stim_ + 1))
        time_ = (initidx(stim_)+1):initidx(stim_ + 1);
        time_ = time_ - time_(1);
        time_ = time_*dt;

        oTimerel(1, (initidx(stim_) ...
            + 1):initidx(stim_ + 1)) = time_;

        if size(tRange, 1) == 1
            t_init_ = find(oTimerel ...
                + tRange(1, 1) == 0);
        else
            t_init_ = find(oTimerel ...
                + tRange(sType(stim_), 1) == 0);
        end

        t_init(stim_) = t_init_(stim_);
        
    else
        
        t_init(stim_) = nan;
        
    end
    
end

for i = 1:numel(sType)
    
    if ~isnan(t_init(i))
        oStimvect(t_init(i):(t_init(i) ...
            + round(sWidth(i, 1)/dt))) = sType(i);
    end
    
end

end