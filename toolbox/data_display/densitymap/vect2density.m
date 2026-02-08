function [densitymap, x_pix, y_pix] = ...
    vect2density(x, y, width, height, ...
    limits, sig, siz, max_norm_flag, ...
    log_norm_flag, sum_norm_flag, ...
    n_norm_flag, tn_norm_flag, input_n)
% vect2density: make a density image of x, y coordinates
%
% Usage:
%   [densitymap, x_pix, y_pix] = ...
%       vect2density(x, y, width, height, ...
%       limits, sig, siz, max_norm_flag, ...
%       log_norm_flag, sum_norm_flag, ...
%       n_norm_flag, tn_norm_flag, input_n)
%
% Args:
%   x, y: two vectors of equal length giving scatterplot x, y coordinates
%   width, height: dimensions of the data density plot, in pixels
%       (default, [64 64])
%   limits: [xmin xmax ymin ymax] defaults to data max/min
%       (default, x, & y max/min)
%   sig, siz: size and sd of 2D smoothing kernel
%       (default, [1 1])
%   max_norm_flag: flag to perform max normalization
%       (default, 0)
%   log_norm_flag: flag to perform log transformation
%       (default, 0)
%   sum_norm_flag: flag to perform sum normalization
%       (default, 0)
%   n_norm_flag: flag to perform n normalization (size of x) 
%       (default, 0)
%   tn_norm_flag: flag to perform normalization to inputed sizes (input_n)
%       (default, 0)
%   input_n: set of input sizes
%
% Notes:
% Inspired by:
%   https://www.mathworks.com/matlabcentral/fileexchange/31726-data-density-plot?focused=3853911&tab=function

if ~exist('limits', 'var') || isempty(limits)
    limits(1) = min(x);
    limits(2) = max(x);
    limits(3) = min(y);
    limits(4) = max(y);
end

if ~exist('width', 'var') || isempty(width)
    width = 64;
end

if ~exist('height', 'var') || isempty(height)
    height = 64;
end

if ~exist('sig', 'var') || isempty(sig)
    sig = 1;
end

if ~exist('siz', 'var') || isempty(siz)
    siz = [1 1];
end

if ~exist('log_norm_flag', 'var') || isempty(log_norm_flag)
    log_norm_flag = 0;
end

if ~exist('max_norm_flag', 'var')|| isempty(max_norm_flag)
    max_norm_flag = 0;
end

if ~exist('sum_norm_flag', 'var')|| isempty(sum_norm_flag)
    sum_norm_flag = 0;
end

if ~exist('n_norm_flag', 'var')|| isempty(n_norm_flag)
    n_norm_flag = 0;
end

if ~exist('tn_norm_flag', 'var')|| isempty(tn_norm_flag)
    tn_norm_flag = 0;
end

if ~exist('input_n', 'var') || isempty(input_n) || isnan(input_n)
    input_n = length(x);
end

deltax = (limits(2) - limits(1))/width;
deltay = (limits(4) - limits(3))/height;

% initialize
densitymap = zeros(height, width);

[densitymap, yx_] = hist3([y', x'], ...
    {limits(3):deltay:limits(4), limits(1):deltax:limits(2)});
%[densitymap, yx_] = hist3([y', x'], [height, width]);

% smooth map (2D kernel)
densitymap = imblur(densitymap, sig, siz, 2);

% normalize density max to 1
if max_norm_flag
    densitymap = ...
        (densitymap - min(densitymap(:)))/max(densitymap(:));
end

% normalize density to sum
if sum_norm_flag
    densitymap = densitymap/sum(densitymap(:));
end

% normalize density to total number of samples (same as sum)
if n_norm_flag
    densitymap = densitymap/length(x);
end

% normalize density to n
if tn_norm_flag
    densitymap = densitymap/input_n;
end

% log transform
if log_norm_flag
    sizY = size(densitymap);
    densitymap = log(densitymap(:));
    densitymap = reshape(densitymap, sizY);
end

y_pix = yx_{1};
x_pix = yx_{2};

end
