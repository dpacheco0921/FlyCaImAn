function wDat_F_hist_per_stack(wDat, ...
    histnorm_flag, nbins, image2replace, figH, axH)
% wDat_F_hist_per_stack: plot histograms of fluorescence 
%   per constituting stacks
%
% Usage:
%   wDatFluoHistPerSegment(wDat, ...
%    histnorm_flag, nbins, image2replace, figH, axH)
%
% Args:
%   wDat: main metadata structure
%   histnorm_flag: flag to normalize histograms
%   nbins: number of planes to use for bining
%   image2replace: input image to use instead of wDat.GreenChaMean
%   figH: figure handle
%   axH: axes handle

if ~exist('figH', 'var') || isempty(figH)
    figH = figure();
end

if ~exist('axH', 'var') || isempty(axH)
    axH = subplot(1, 1, 1);
end

if exist('image2replace', 'var') || ~isempty(image2replace)
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
    
    cormap_int = ...
        colormapgen(2, max(wDat.Zstitch.Zidx));
    
    for i = 1:numel(intrange)-1
        
        subv = wDat.GreenChaMean(:, :, intrange(i)+1: intrange(i+1));
        [lhist, ~] = histogram(subv(:), hbins);
        
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
    
    cormap_int = colormapgen(2, max(wDat.Zstitch.Zidx));
    
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

xlim(axH, [-20 50])
legend(axH, lh, lstr)
title(axH, 'Histogram of fluorescent values')

end