function [CorrData, cDat, lStim] = ...
    LEDdenoiser(Data, Baseline, iDat, ...
    lStim, axH, cDat, plotgate)
% LEDdenoiser: This function perfomrs LED noise correction on imaging data
%
% Usage:
%   [CorrData, cDat, lStim] = ...
%       LEDdenoiser(Data, Baseline, iDat, ...
%       lStim, axH, cDat, plotgate)
%
% Args:
%   Data: 2DxT or 3DxT dataset
%   Baseline: time without opto stimulation
%   iDat: image metadata structure
%   lStim: stimuli metadata structure
%   axH: axes handles [1 x 3] used to plot this iteration (only used if plotgate)
%   cDat: LED correction metadata structure (needs Iter, buffer, CorType)
%   plotgate: flag to plot process
%
% Returns:
%   CorrData: LED-corrected data, same dimensions as Data except the
%       height, which is 2 pixels shorter (last rows are dropped)
%   cDat: updated correction structure (light triggered median / average
%       before and after correction, per iteration, in fluorescence units)
%   lStim: updated stimuli structure (lStim2D = stimulus interpolated to
%       pixel times, pstEn = pulse start / end in pixel-time samples)
%
% Notes:
% updates cDat with LTA, LTM per iteration
% LightTriggeredMedian = LTM, LightTriggeredAverage = LTA

if ~exist('Baseline','var') || isempty(Baseline); Baseline = 1:40; end
if ~exist('plotgate','var') || isempty(plotgate); plotgate = 0; end
if exist('flip', 'builtin')
    str2use = 'flip';
else
    str2use = 'flipdim';
end

% normalizing / smoothing data and correcting video shape
cDat.minInit = min(Data(:));
[DataNorm2DxT, DataMedian, oDim] = ...
    preprocessingIm(Data, Baseline, iDat, lStim, plotgate);
clear Data
wDim = size(DataNorm2DxT);

% Get Stim info (only first iteration) and chopping frame
if ~isfield(lStim, 'lStim2D') && cDat.Iter == 1
    [lStim.lStim2D, DataNorm2DxT, DataMedian] = ...
        lStimInterp(wDim, oDim, iDat, lStim, DataNorm2DxT, DataMedian, plotgate); 
end

if ~isempty(iDat.DelFrames) 
    
    % Discarding frames that are empty
    %   (lStim2D is only chopped if it still has all frames: after the
    %   first iteration it was already chopped, and chopping it again
    %   would remove the wrong frames)
    PrelStim2D = lStim.lStim2D;
    if size(lStim.lStim2D, 3) == size(DataNorm2DxT, 3)
        lStim.lStim2D(:, :, iDat.DelFrames) = [];
    end
    PreNorm2D = DataNorm2DxT;
    DataNorm2DxT(:, :, iDat.DelFrames) = [];
    
end

% updating size, in case it was changed
oDim(1) = size(DataNorm2DxT, 2);
wDim = size(DataNorm2DxT);

lDim = size(lStim.lStim2D);
lStim1D = reshape(lStim.lStim2D, [prod(lDim), 1]);

% vectorize DataNorm2DxT
DataNorm1D = reshape(DataNorm2DxT, [prod(wDim), 1]);

% Getting light triggered average (lta)
lPulse = lStim1D > 0; lpstEn = diff(lPulse);
lStim.pstEn = []; lStim.pstEn = [find(lpstEn == 1) + 1, find(lpstEn == -1)];
clear lpstEn lPulse
LT_matrix = ltagen(lStim, DataNorm1D, cDat.buffer);
cDat.LTM(cDat.Iter, :) = median(LT_matrix);
cDat.LTM(cDat.Iter, :) = cDat.LTM(cDat.Iter, :) - mean(cDat.LTM(cDat.Iter, 1:cDat.buffer));
cDat.LTA(cDat.Iter, :) = mean(LT_matrix);
cDat.LTA(cDat.Iter, :) = cDat.LTA(cDat.Iter, :) - mean(cDat.LTA(cDat.Iter, 1:cDat.buffer));
clear lta

% Correcting video
eval(['[CorDataNorm2DxT, LT_matrix] = ', ...
    'videoCorr(DataNorm1D, lStim1D, lStim, cDat.', ...
    cDat.CorType, '(cDat.Iter, :), cDat.buffer, wDim);'])
plottingCorrVideo(CorDataNorm2DxT, DataNorm2DxT, lStim.lStim2D, plotgate);
cDat.LTMCor(cDat.Iter, :) = median(LT_matrix);
cDat.LTMCor(cDat.Iter, :) = cDat.LTMCor(cDat.Iter, :) - mean(cDat.LTMCor(cDat.Iter, 1:cDat.buffer));
cDat.LTACor(cDat.Iter, :) = mean(LT_matrix);
cDat.LTACor(cDat.Iter, :) = cDat.LTACor(cDat.Iter, :) - mean(cDat.LTACor(cDat.Iter, 1:cDat.buffer));
clear CorDataNorm1D  lta

% Storing images to plot
cDat.XYprojRaw{cDat.Iter} = squeeze(median(DataNorm2DxT, 1)); clear DataNorm2DxT
cDat.XYprojRawCor{cDat.Iter} = squeeze(median(CorDataNorm2DxT, 1));

% Plotting results
% plot mean and median
if plotgate; resultplotter(cDat, axH); end

% Adding discarded frames to volume
if  ~isempty(iDat.DelFrames)
    CorrFrames = ones(1, size(PreNorm2D, 3));
    CorrFrames(1, iDat.DelFrames) = 0;
    PreNorm2D = zeros([wDim(1:2) prod(oDim(3:end))]);
    PreNorm2D(:, :, CorrFrames == 1) = CorDataNorm2DxT;
    clear CorDataNorm2DxT
    CorDataNorm2DxT = PreNorm2D;
    clear PreNorm2D CorrFrames
end

% reshaping data
if length(oDim) == 3
    FramesPerTrial = floor(iDat.StackN/lStim.trialn);
    eval(['CorDataNorm2DxT(:, 1:2:end, :) = ', ...
        str2use, '(CorDataNorm2DxT(: , 1:2:end, :), 1);'])
    DataNorm4D = reshape(CorDataNorm2DxT, [oDim(2), oDim(1), oDim(3:end)]);
    
    % add background
    for trialn = 1:lStim.trialn
        if trialn == lStim.trialn
            CorrData(:, :, (FramesPerTrial*(trialn-1) + 1):iDat.StackN) = ...
                bsxfun(@plus, DataNorm4D(:, :, ...
                (FramesPerTrial*(trialn-1) + 1):iDat.StackN), DataMedian{trialn});
        else
            CorrData(:, :, (FramesPerTrial*(trialn-1) + 1):FramesPerTrial*trialn) = ...
                bsxfun(@plus, DataNorm4D(:, :, ...
                (FramesPerTrial*(trialn-1) + 1):FramesPerTrial*trialn), DataMedian{trialn});
        end
    end
    
    eval(['CorrData = permute(',str2use,'(CorrData, 1), [2 1 3]);'])
    
    % replacing empty frames
    if  ~isempty(iDat.DelFrames)
    	CorrData = framegapfill(iDat.DelFrames, CorrData);
    end
    
    cDat.minEnd = min(CorrData(:));
    
else
    
    StacksPerTrial = iDat.StackN/lStim.trialn;
    eval(['CorDataNorm2DxT(:, 1:2:end, :) = ', ...
        str2use, '(CorDataNorm2DxT(: , 1:2:end, :), 1);'])
    DataNorm4D = reshape(CorDataNorm2DxT, [oDim(2), oDim(1), oDim(3:end)]);
    
    % add background
    for trialn = 1:lStim.trialn
        CorrData(:, :, :, (StacksPerTrial*(trialn-1) + 1):StacksPerTrial*trialn) = ...
            bsxfun(@plus, DataNorm4D(:, :, :, ...
            (StacksPerTrial*(trialn-1) + 1):StacksPerTrial*trialn), DataMedian{trialn});
    end
    
    eval(['CorrData = permute(',str2use,'(CorrData, 1), [2 1 3:length(oDim)]);'])
    
    % replacing empty frames
    if ~isempty(iDat.DelFrames)
    	CorrData = OptVolFiller(iDat.DelFrames, CorrData);
    end
    
    cDat.minEnd = min(CorrData(:));
    
end

end

function Data = OptVolFiller(frame2corr, Data)
% OptVolFiller: fill whole volumes that contain empty frames
%
% Args:
%   frame2corr: indexes of empty frames, counted over all planes and
%       timepoints of the raw data (plane fastest)
%   Data: 4D [Y X Z T] data
%
% Returns:
%   Data: same size, volumes with any empty frame replaced (framegapfill)

% assumes that Data has the same planes and timepoints as the original (raw)
Idx_1D = zeros(1, size(Data, 3)*size(Data, 4));
Idx_1D(frame2corr) = 1;
Idx_1D = sum(reshape(Idx_1D, [size(Data, 3) size(Data, 4)]), 1);
vol2corr = find(Idx_1D > 0); clear Idx_1D
Data = framegapfill(vol2corr, Data);
fprintf('\n')

end

function [DataNorm2DxT, DataMedian, oDim] = preprocessingIm(Data, Baseline, iDat, lStim, plotgate)
% preprocessingIm: subtract the baseline median image (per trial), undo
%   the bidirectional scanning flip and reshape to [X Y frames]
%
% Args:
%   Data: [Y X (Z) T] imaging data; Baseline: timepoints without light
%   iDat, lStim: metadata (StackN, trialn); plotgate: flag to plot
%
% Returns:
%   DataNorm2DxT: baseline-subtracted data [X Y frames]
%   DataMedian: cell (one per trial) with the subtracted median images
%   oDim: size of the input Data

% Smoothing / reshaping Data
Data = double(Data);
oDim = size(Data);
Data = permute(Data, [2 1 3:length(oDim)]);
if exist('flip', 'builtin')
    str2use = 'flip';
else
    str2use = 'flipdim';
end
eval(['Data = ',str2use,'(Data, 1);']);

% whitening video
if length(oDim) == 3
    
    FramesPerTrial = floor(iDat.StackN/lStim.trialn);
    
    % median substract
    for trialn = 1:lStim.trialn
        DataMedian{trialn} = median(Data(:, :, FramesPerTrial*(trialn-1) + Baseline), 3);
        if trialn == lStim.trialn
            DataNorm2DxT(:, :, (FramesPerTrial*(trialn-1) + 1):iDat.StackN) ...
                = bsxfun(@minus, Data(:, :, ...
                (FramesPerTrial*(trialn-1) + 1):iDat.StackN), DataMedian{trialn});
        else
            DataNorm2DxT(:, :, (FramesPerTrial*(trialn-1) + 1):FramesPerTrial*trialn) ...
                = bsxfun(@minus, Data(:, :, ...
                (FramesPerTrial*(trialn-1) + 1):FramesPerTrial*trialn), DataMedian{trialn});            
        end
    end
    
    % Correct for bidirectionality
    eval(['DataNorm2DxT(:, 1:2:end, :) = ', ...
        str2use, '(DataNorm2DxT(: , 1:2:end, :), 1);'])
    
else
    
    StacksPerTrial = iDat.StackN/lStim.trialn;
    
    % median substract
    for trialn = 1:lStim.trialn
        DataMedian{trialn} = median(Data(:, :, :, StacksPerTrial*(trialn-1) + Baseline), 4);
        DataNorm4D(:, :, :, (StacksPerTrial*(trialn-1) + 1):StacksPerTrial*trialn) ...
            = bsxfun(@minus, Data(:, :, :, ...
            (StacksPerTrial*(trialn-1) + 1):StacksPerTrial*trialn), DataMedian{trialn});        
    end
    DataNorm2DxT = reshape(DataNorm4D, [oDim([2 1]), prod(oDim(3:end))]);
    clear DataNorm4D
    
    % Correct for bidirectionality
    eval(['DataNorm2DxT(:, 1:2:end, :) = ', ...
        str2use, '(DataNorm2DxT(: , 1:2:end, :), 1);'])
    
end

% Plot image
if plotgate == 1
    
    figure('position', [397 694 1006 374])
    
    for i = 1:size(Data, 3)
        
       subplot(1, 2, 1);
       imagesc(Data(:, :, i));
       caxis([0 50]);
       title(num2str(i));
       subplot(1, 2, 2);
       imagesc(DataNorm2DxT(:, :, i));
       caxis([0 50]);
       title(num2str(i));
       pause(0.05)
       
    end
    
end

end

function lta = ltagen(lStim, DataNorm1D, buffer)
% ltagen: light triggered traces, one row per light pulse of the most
%   common width, each with 'buffer' samples before and after the pulse
%
% Args:
%   lStim: stimuli structure (pstEn: pulse start / end samples)
%   DataNorm1D: data as one long vector in pixel-time order
%   buffer: samples to include before and after each pulse
%
% Returns:
%   lta: [pulses x (pulse width + 2*buffer + 1)] smoothed traces

% Selecting the median pulse width and use those for pulse lta
ModPulse = mode(lStim.pstEn(:, 2) - lStim.pstEn(:, 1));
MedIdx = find((lStim.pstEn(:, 2) - lStim.pstEn(:, 1)) == ModPulse);
MedpstEn = lStim.pstEn(MedIdx, :);

% it seems that without smoothing is working fine so far
lta = zeros(numel(MedIdx), (ModPulse + buffer*2 + 1));

for p_idx = 1:size(MedpstEn, 1)
    LocalStore = DataNorm1D((MedpstEn(p_idx, 1)-buffer):(MedpstEn(p_idx, 2)+buffer));
    lta(p_idx, 1:numel(LocalStore)) = smooth(LocalStore, 6);
end

end

function [lStim2DxT, DataNorm2DxT, DataMedian] = ...
    lStimInterp(dDim, oDim, iDat, lStim, DataNorm2DxT, DataMedian, plotgate)
% lStimInterp: interpolate the light stimulus trace to the time each pixel
%   was scanned, and drop the last 2 columns of data / stimulus / medians
%   (same pruning is applied to the red channel in batch_LEDdenoise)
%
% Args:
%   dDim: size of DataNorm2DxT; oDim: original data size
%   iDat: image metadata (fstEn: frame start / end samples)
%   lStim: stimuli structure (trace); DataMedian: cell of median images
%
% Returns:
%   lStim2DxT: stimulus value per pixel [X-2 Y frames]
%   DataNorm2DxT, DataMedian: inputs without the last 2 columns

% Interpolate LED pulses to frame times and reduces frame size (x-2 pixels)
PixelSpacing = prod(dDim(1:2));
lStim2DxT = [];

for fIdx = 1:dDim(3)
    %fprintf([num2str(fIdx), '\n'])
    TimePerFrame = diff(iDat.fstEn(fIdx,:)) + 1;
    % Interpolate stim
    LinearStim = interp1((0:1:(TimePerFrame-1))/((TimePerFrame-1)), ...
        lStim.trace(iDat.fstEn(fIdx, 1):iDat.fstEn(fIdx, 2)), (0:1:(PixelSpacing-1))/(PixelSpacing-1));
    lStim2DxT(:, :, fIdx) = reshape(LinearStim, [dDim(1), dDim(2)])';
    clear preStim2D
end

lStim2DxT = permute(lStim2DxT, [2 1 3]);

% Plot image
if plotgate == 1
    for i = 1:dDim(3)
       subplot(1, 2, 1); imagesc(DataNorm2DxT(:, :, i));
       caxis([0 50])
       subplot(1, 2, 2); imagesc(lStim2DxT(:, :, i));
       caxis([0 2]); title(num2str(i));
       pause(0.05)
    end
end

lStim2DxT = lStim2DxT(:, 1:end-2, :);
DataNorm2DxT = DataNorm2DxT(:, 1:end-2, :);

for i = 1:numel(DataMedian)
    DataMedian{i} = DataMedian{i}(:, 1:end-2, :);
end

end

function [CorDataNorm2DxT, lta] = ...
    videoCorr(DataNorm1D, lStim1D, lStim, corrVector, buffer, dDim)
% videoCorr: subtract the light triggered median (corrVector) around
%   every light pulse; pulses cut by the start / end of the recording or
%   of different width use the matching part of corrVector.
%   The hard-coded 300 and 50 sample values trim the stitching region
%   between the two ends of corrVector and the pulse edges (set to 0).
%
% Returns:
%   CorDataNorm2DxT: corrected data [dDim]; lta: light triggered traces
%   of the corrected data

% This function perfomrs video correction (substraction)
CorDataNorm1D = DataNorm1D;

% Getting single ended corrVectors, for stitching
corrVector_a = corrVector(300:end); corrVector_b = corrVector(1:end-300);
maxlength = numel(corrVector_b)+ numel(corrVector_a);

for p_idx = 1:size(lStim.pstEn, 1)
    
    % getting correct corrVector size
    timeIdx = (lStim.pstEn(p_idx, 1) - buffer):(lStim.pstEn(p_idx, 2) + buffer);
    Interval = numel(timeIdx);
    halfWay = ceil(Interval/2);
    if Interval == numel(corrVector)
        lta2sub = corrVector(1:Interval);
        lta2sub([1:50 end-50:end]) = 0;
    else
        if min(lStim1D(lStim.pstEn(p_idx, 1):lStim.pstEn(p_idx, 1)+30)) == max(lStim1D)
            
            % When pulse is chopped at the begining
            timeIdx = lStim.pstEn(p_idx, 1):(lStim.pstEn(p_idx, 2) + buffer);
            Interval = numel(timeIdx);
            lta2sub = corrVector((end - Interval + 1):end);
            lta2sub(end-50:end) = 0;
            
        elseif min(lStim1D(lStim.pstEn(p_idx, 2)-30:lStim.pstEn(p_idx, 2))) == max(lStim1D)
            
            % When pulse is chopped at the end
            timeIdx = (lStim.pstEn(p_idx, 1) - buffer):lStim.pstEn(p_idx, 2);
            Interval = numel(timeIdx);
            lta2sub = corrVector(1:Interval);
            lta2sub(1:50) = 0;
            
        else % pulse is slightly chopped or lenghtend
            
            try
                
                if Interval > maxlength
                    lta2sub = corrVector_b; 
                    lta2sub(numel(corrVector_b)+ (1:(Interval - maxlength))) = ...
                        corrVector(halfWay + (1:(Interval - maxlength)));
                    lta2sub((numel(corrVector_b)+ Interval - maxlength + 1):Interval) = ...
                        corrVector_a;
                else
                    lta2sub(1:halfWay) = corrVector_b(1:halfWay);
                    lta2sub(halfWay + 1:Interval) = ...
                        corrVector_a((end -(Interval-halfWay) + 1):end);
                end
                
                lta2sub([1:50 end-50:end]) = 0;
                
            catch ME

                % lta2sub is not defined if this failed, so stop here
                error('LEDdenoiser:videoCorr', ...
                    'Pulse %d could not be corrected: %s', p_idx, ME.message)

            end
            
        end
    end
    
    % Substracting
    CorDataNorm1D(timeIdx) = CorDataNorm1D(timeIdx) - lta2sub';
    clear lta2sub timeIdx
    
end

lta = ltagen(lStim, CorDataNorm1D, buffer);
CorDataNorm2DxT = ...
    reshape(CorDataNorm1D, [dDim(1), dDim(2), dDim(3)]);

end

function plottingCorrVideo(CorDataNorm2DxT, DataNorm2DxT, lStim2D, plotgate)
% plottingCorrVideo: (if plotgate == 1) play raw, corrected and stimulus
%   frames side by side; no return value

if plotgate == 1
    
    figure();
    
    for p_idx = 1:size(CorDataNorm2DxT, 3)
        
       subplot(1, 3, 1)
       imagesc(DataNorm2DxT(:, :, p_idx)); caxis([0 50])
       subplot(1, 3, 2)
       imagesc(CorDataNorm2DxT(:, :, p_idx)); caxis([0 50])
       subplot(1, 3, 3)
       imagesc(lStim2D(:, :, p_idx)); caxis([0 2])
       title(num2str(p_idx)); pause(0.05)
       
    end
    
end

end

function resultplotter(cDat, axH)
% resultplotter: plot light triggered average / median before and after
%   correction (axH(1)) and the raw / corrected YX projections (axH(2:3));
%   no return value

% LTA
plot(cDat.LTM(cDat.Iter, :), 'b', 'Parent', axH(1)); hold(axH(1), 'on'); 
plot(cDat.LTA(cDat.Iter, :), 'r', 'Parent', axH(1))
plot(cDat.LTMCor(cDat.Iter, :), 'c', 'Parent', axH(1)); 
plot(cDat.LTACor(cDat.Iter, :), 'm', 'Parent', axH(1))
title(axH(1), 'Light triggered avegare')
xlabel(axH(1), 'Time (pixels units)')
ylabel(axH(1), 'Fluorescence (a.u)')
xlim(axH(1), [1 size(cDat.LTACor, 2)])
ylim(axH(1), [min([min(cDat.LTA)*1.1, -2]) max([max(cDat.LTA)*1.1, 2])])

% 2D Pre
imagesc(cDat.XYprojRaw{cDat.Iter}, 'Parent', axH(2))
title(axH(2), 'Raw movie')
xlabel(axH(2), 'Time (frames)')
ylabel(axH(2), 'YX projection (pixels)')
caxis(axH(2), [0 max(cDat.XYprojRaw{cDat.Iter}(:))])
colorbar(axH(2))

% 2D Post
imagesc(cDat.XYprojRawCor{cDat.Iter}, 'Parent', axH(3))
title(axH(3), 'Corrected movie')
xlabel(axH(3), 'Time (frames)')
ylabel(axH(3), 'YX projection (pixels)')
colorbar(axH(3))
caxis(axH(3), [0 max(cDat.XYprojRawCor{cDat.Iter}(:))])

end
