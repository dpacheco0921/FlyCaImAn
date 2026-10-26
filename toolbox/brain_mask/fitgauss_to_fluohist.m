function wDat = fitgauss_to_fluohist(wDat, ...
    lnbis, gaussnumber, sigmaRed, ...
    gauss2use, image2replace, figH, ...
    axH, plot_flag)
% fitgauss_to_fluohist: fit gaussian to distribution of 
%   fluorescence per stack
%
% Usage:
%   wDat = fitgauss_to_fluohist(wDat, ...
%       lnbis, gaussnumber, sigmaRed, ...
%       gauss2use, image2replace, figH, ...
%       axH, plot_flag)
%
% Args:
%   wDat: main metadata structure
%   lnbis: stack edges in z axis (plane indices, ascending, first is 0)
%   gaussnumber: number of gaussian to fit (gauss1, gauss2, etc)
%       (default, 'gauss5')
%   sigmaRed: sigma, number of widths [left right] around the gaussian
%       center used for the threshold
%       (default, [2 5])
%   gauss2use: gaussian number to retrieve
%       (default, 1)
%   image2replace: input image to use instead of wDat.GreenChaMean
%       (default, [], use wDat.GreenChaMean)
%   figH: figure handle
%   axH: axes handle
%   plot_flag: flag to plot results
%       (default, 1)
%
% Returns:
%   wDat: input structure with
%       wDat.bF_zbin: stack edges (same as lnbis)
%       wDat.bF_ths: 1 x (number of stacks) intensity threshold (a.u.),
%           the upper threshold (center + sigmaRed(2) widths) of the
%           selected gaussian fitted to each stack's fluorescence histogram

if ~exist('image2replace', 'var'); image2replace = []; end

if numel(lnbis) < 2
    error('fitgauss_to_fluohist: lnbis needs at least two stack edges')
end

if ~exist('gaussnumber', 'var') || isempty(gaussnumber)
    gaussnumber = 'gauss5';
end

if ~exist('sigmaRed', 'var') || isempty(sigmaRed)
    sigmaRed = [2 5];
end

if ~exist('gauss2use', 'var') || isempty(gauss2use)
    gauss2use = 1;
end

if ~exist('plot_flag', 'var') || isempty(plot_flag)
    plot_flag = 1;
end

if plot_flag
    
    if ~exist('figH', 'var') || isempty(figH)
        figH = figure();
    end
    
    if ~exist('axH', 'var') || isempty(axH)
        axH = subplot(1, 1, 1);
    end
    
end

% variables to update
wDat.bF_ths = [];
wDat.bF_zbin = lnbis;

for i = 1:numel(wDat.bF_zbin)-1
    
    % decide if using internal or external mean image
    if ~isempty(image2replace)
        subv = image2replace(:, :, ...
            wDat.bF_zbin(i)+1: wDat.bF_zbin(i+1));
    else
        subv = wDat.GreenChaMean(:, :, ...
            wDat.bF_zbin(i)+1: wDat.bF_zbin(i+1));
    end
    
    hDat = fitgauss1D(subv(:), ...
        gaussnumber, sigmaRed, gauss2use, 0.01);
    
    % plot histograms
    if plot_flag
        
        fitgaussPlotter(hDat, figH, axH);
        title(axH, 'Gaussian fit to fluorescence distribution')
        
        if ~isempty(image2replace)
            xlim(axH, [min(image2replace(:)), ...
                max(image2replace(:))])
        else
            xlim(axH, [min(wDat.GreenChaMean(:)), ...
                max(wDat.GreenChaMean(:))])
        end
        
    end
    
    wDat.bF_ths(1, i) = hDat.xths(1, 2);
    clear subv hDat
    
end

end
