function xyz = to_plot_point_3D(obj, i_idx, imatrix)
% to_plot_point_3D: plot selected indeces per trace
%
% Usage:
%   xyz = to_plot_point_3D(obj, i_idx, imatrix)
%
% Args:
%   obj: traceObj obj
%   i_idx: index to use per trace ([n by trees]: row r, column t is the
%       point index to plot in tree t; NaN entries are skipped)
%   imatrix: field to use
%
% Returns:
%   xyz: [n*trees by 3] XYZ coordinates of the plotted points (same units
%       as the traces), also drawn as blue circles in the current axes
%
% Notes: see traceObj

% extract xyz matrix (cell with one [points by 3] matrix per tree)
xyz_trees = to_pick_xyz_object(obj, [], imatrix);

% each column of i_idx belongs to one tree
assert(size(i_idx, 2) == numel(xyz_trees), ...
    'to_plot_point_3D: i_idx needs one column per tree (%d), got %d', ...
    numel(xyz_trees), size(i_idx, 2))

% get single points per matrix
xyz = [];

for i_tree = 1:numel(xyz_trees)

    idx = i_idx(:, i_tree);
    idx(isnan(idx)) = [];
    xyz = [xyz; xyz_trees{i_tree}(idx, :)]; %#ok<AGROW>

end

scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), ...
    'o', 'MarkerEdgeColor', [0 0 1])

end
