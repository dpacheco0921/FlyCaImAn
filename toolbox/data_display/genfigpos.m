function figpos = genfigpos(monitor_id, screen_position, figure_siz)
% genfigpos: function that provides with positions for figures relative to
%   current screen(s)
%
% Usage:
%   figpos = genfigpos(monitor_id, screen_position, figure_siz)
%
% Args:
%   monitor_id: which monitor to use
%       (1, main screen, default)
%       (2, secondary screen)
%   screen_position: center, NW, NE, SW, SE, [] (means the whole screen)
%       (default, 'center')
%   figure_siz: figure size [width height]
%       (default, [300 300])
% 
% Returns:
%   figpos: figure position

if ~exist('monitor_id', 'var') || isempty(monitor_id)
    monitor_id = 1;
end
if ~exist('figure_siz', 'var') || isempty(figure_siz)
    figure_siz = [300 300];
end
if ~exist('screen_position', 'var') || isempty(screen_position)
    screen_position = 'center';
end

% get monitor position
mon_pos = get(0, 'MonitorPositions');

if ispc
   mon_pos(:, 1:2) = mon_pos(:, 1:2) + 60;
   mon_pos(:, 3) = mon_pos(:, 3) - 130;
   mon_pos(:, 4) = mon_pos(:, 4) - 150;
end

switch screen_position
    
    case 'center'
        
        mon_center(:, 1) = mon_pos(:, 1) + mon_pos(:, 3)/2;
        mon_center(:, 2) = mon_pos(:, 2) + mon_pos(:, 4)/2;
        figpos = [mon_center(monitor_id, 1) - (figure_siz(1)/2), ...
            mon_center(monitor_id, 2) - (figure_siz(2)/2), ...
            figure_siz];
        
    case 'nw'
        
        mon_center = mon_pos;
        mon_center(:, 2) = mon_center(:, 2) + mon_center(:, 4) - figure_siz(2);
        mon_center(:, 3) = figure_siz(1); mon_center(:, 4) = figure_siz(2);
        figpos = mon_center(monitor_id, :);
        
    case 'ne'
        
        mon_center = mon_pos;
        mon_center = mon_center(:, 1:2) + mon_center(:, 3:4);
        mon_center(:, 1) = mon_center(:, 1) - figure_siz(1);
        mon_center(:, 2) = mon_center(:, 2) - figure_siz(2);
        mon_center(:, 3) = figure_siz(1); mon_center(:, 4) = figure_siz(2);
        figpos = mon_center(monitor_id, :);
        
    case 'se'
        
        mon_center = mon_pos;
        mon_center(:, 1) = mon_center(:, 1) + mon_center(:, 3) - figure_siz(1);
        mon_center(:, 3) = figure_siz(1); mon_center(:, 4) = figure_siz(2);
        figpos = mon_center(monitor_id, :);
        
    case 'sw'
        
        mon_center = mon_pos;
        mon_center(:, 3) = figure_siz(1); mon_center(:, 4) = figure_siz(2);
        figpos = mon_center(monitor_id, :);
        
    otherwise % whole screen
        figpos = mon_pos(monitor_id, :);
        
end

end