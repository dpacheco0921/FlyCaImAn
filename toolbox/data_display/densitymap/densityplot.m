function [figH, axH, densitymap, x_pix, y_pix] = ...
    densityplot(x, y, width, height, limits, ...
    sig, siz, axH, log_norm_flag, max_norm_flag)
% densityplot: collects values from fields of a given structure
%
% Usage:
%   [figH, axH, densitymap, x_pix, y_pix] = ...
%      densityplot(x, y, width, height, ...
%      limits, sig, siz, axH, log_norm_flag, max_norm_flag)
%
% Args:
%   x, y: two vectors of equal length providing scatterplot x, y coordinates
%   width, height: dimensions of the data density plot, in pixels
%       (default, [])
%   limits: [xmin xmax ymin ymax]
%       (default, [], it internally will set it to x, & y max/min)
%   sig, siz: size and sd of 2D smoothing kernel
%   axH: input axis handle
%   log_norm_flag: flag to perform log transformation
%       (default, 0)
%   max_norm_flag: flag to perform max normalization
%       (default, 0)
%
% Notes:
% see vect2density

if ~exist('width', 'var'); width = []; end
if ~exist('height', 'var'); height = []; end
if ~exist('limits', 'var'); limits = []; end

if ~exist('log_norm_flag', 'var')
    log_norm_flag = 0;
end

if ~exist('max_norm_flag', 'var')
    max_norm_flag = 0;
end

% get density map
[densitymap, x_pix, y_pix] = ...
    vect2density(x, y, width, height, ...
    limits, sig, siz, log_norm_flag, max_norm_flag);

% plot density map
if ~exist('axH', 'var') || isempty(axH)
    figH = figure();
    axH = subplot(1, 1, 1); 
end

imagesc(x_pix, y_pix, densitymap, 'Parent', axH);

set(axH, 'YDir', 'normal')

end
