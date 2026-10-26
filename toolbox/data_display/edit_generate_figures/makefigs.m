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
%   y_subplot: number of subplot rows (first argument of subplot)
%       (default, 5)
%   x_subplots: number of subplot columns
%       (default, 5)
%   screen_position: position of figure
%       ('center', 'nw', 'ne', 'sw', 'se' (lowercase), or other/[] for the whole screen)
%       (default, 'center')
%   figure_siz: figure size [width height]
%       (default, [], genfigpos then uses [300 300])
%
% Returns:
%   figH: figure handle
%   axH: vector of axes handles, one per subplot (row by row)

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
