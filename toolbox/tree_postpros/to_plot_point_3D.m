function xyz = to_plot_point_3D(obj, i_idx, imatrix)
% to_plot_point_3D: plot selected indeces per trace
%
% Usage:
%   xyz = to_plot_point_3D(obj, i_idx, imatrix)
%
% Args:
%   obj: traceObj obj
%   i_idx: index to use per trace
%   imatrix: field to use

% extract xyz matrix
xyz = to_pick_xyz_object(obj, [], imatrix);

% get single points per matrix
xyz = cellfun(@(x) x(y, :), ...
    xyz, chunk2cell(i_idx, 1));

scatter3(i_idx(:, 1), i_idx(:, 2), i_idx(:, 3), ...
    'o', 'MarkerEdgeColor', [0 0 1])

end
