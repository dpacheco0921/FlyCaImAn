function [figH, axH] = makefigs(y_subplot, x_subplots, ...
    figure_siz, screen_position)
% makefigs: function that generates figure
%   providing figure handles and axes handles.
%
% Usage:
%   [figH, axH] = makefigs(y_subplot, x_subplots, ...
%       figure_siz, screen_position)
%
% Args:
%   y_subplot: number of columns
%       (default, 5)
%   x_subplots: number of rows
%       (default, 5)
%   screen_position: position of figure
%       (center, NW, NE, SW, SE, [] (whole screen))
%       (default, 'center')
%   figure_siz: figure size [width height]
%       (default, [])

if ~exist('y_subplot', 'var') || isempty(y_subplot)
    y_subplot = 5;
end

if ~exist('x_subplots', 'var') || isempty(x_subplots)
    x_subplots = 5;
end

if ~exist('figure_siz', 'var') || isempty(figure_siz)
    figure_siz = [];
end

if ~exist('screen_position', 'var') || isempty(screen_position)
    screen_position = 'center';
end

figH = figure('Position', ...
    genfigpos(1, screen_position, figure_siz));

% generate axes
for i = 1:y_subplot*x_subplots
    axH(i) = subplot(y_subplot, x_subplots, i);
end

end