function wDat_F_hist_per_stack(wDat, ...
    histnorm_flag, nbins, image2replace, figH, axH)
% wDat_F_hist_per_stack: plot histograms of fluorescence 
%   per constituting stacks
%
% Usage:
%   wDat_F_hist_per_stack(wDat, ...
%    histnorm_flag, nbins, image2replace, figH, axH)
%
% Args:
%   wDat: main metadata structure
%   histnorm_flag: flag to normalize histograms (1) or not (0)
%   nbins: number of planes to use for bining (scalar), or a vector of
%       plane edges. If empty, wDat.Zstitch.Zidx defines the stacks.
%   image2replace: input image to use instead of wDat.GreenChaMean
%       (default, [], use wDat.GreenChaMean)
%   figH: figure handle
%   axH: axes handle
%
% Returns:
%   nothing, it only draws one histogram line per stack in axH
%   (x axis: fluorescence in a.u., y axis: voxel counts)

if ~exist('figH', 'var') || isempty(figH)
    figH = figure();
end

if ~exist('axH', 'var') || isempty(axH)
    axH = subplot(1, 1, 1);
end

% only replace the image if one was actually provided (an empty input would
%   otherwise wipe wDat.GreenChaMean)
if exist('image2replace', 'var') && ~isempty(image2replace)
    wDat.GreenChaMean = image2replace;
end

if ~exist('histnorm_flag', 'var') || isempty(histnorm_flag)
    histnorm_flag = 0;
end

hbins = -2*10^2:0.01:3*10^3;
lh = []; lstr = [];

if exist('nbins', 'var') && ~isempty(nbins)
    
    if numel(nbins) == 1
        intrange = 0:nbins:wDat.vSize(3);
        intrange(end) = wDat.vSize(3);
    else
        intrange = nbins;
    end
    
    cormap_int = jet(max(wDat.Zstitch.Zidx));
    
    for i = 1:numel(intrange)-1
        
        subv = wDat.GreenChaMean(:, :, intrange(i)+1: intrange(i+1));
        [lhist, ~] = hist(subv(:), hbins);
        
        if histnorm_flag
            lhist = lhist/max(lhist);
        end
        
        lh(i) = plot(hbins, lhist, ...
            'color', cormap_int(i, :), ...
            'Parent', axH);
        lstr{i} = num2str(i);
        hold(axH, 'on')
        clear subv
        
    end
    
else
    
    if isfield(wDat, 'plane2keep')
        wDat.Zstitch.Zidx = ...
            wDat.Zstitch.Zidx(wDat.plane2keep ~= 0);
    end
    
    cormap_int = jet(max(wDat.Zstitch.Zidx));
    
    for i = 1:max(wDat.Zstitch.Zidx)
        
        subv = ...
            wDat.GreenChaMean(:, :, wDat.Zstitch.Zidx == i);
        [lhist, ~] = hist(subv(:), hbins);
        
        if histnorm_flag
            lhist = lhist/max(lhist);
        end
        
        lh(i) = plot(hbins, lhist, ...
            'color', cormap_int(i, :), ...
            'Parent', axH);
        lstr{i} = num2str(i);
        hold(axH, 'on')
        clear subv
        
    end
    
end

xlim(axH, [min(wDat.GreenChaMean(:)), ...
    max(wDat.GreenChaMean(:))])
legend(axH, lh, lstr)
title(axH, 'Histogram of fluorescent values')

end
