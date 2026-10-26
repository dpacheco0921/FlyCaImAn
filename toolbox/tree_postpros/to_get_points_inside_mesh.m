function [xyz, idx_inside] = to_get_points_inside_mesh(X, Y, Z, i_mesh, chunksiz)
% to_get_points_inside_mesh: function that finds XYZ coordinates that are
% inside a i_mesh
% Usage:
%   [xyz, idx_inside] = to_get_points_inside_mesh(X, Y, Z, i_mesh, chunksiz)
%
% Args:
%   X, Y, Z: point coordinates (same size, any shape)
%   i_mesh: input mesh (treestoolbox 'elem' format or faces/vertices)
%   chunksiz: number of points tested at a time (default, 1000)
%
% Returns:
%   xyz: [m by 3] coordinates of the points inside the mesh
%   idx_inside: [1 by m] linear index (into X(:)) of those points
%
% Notes: see 'treestoolbox' for format info

if ~exist('chunksiz', 'var') || isempty(chunksiz); chunksiz = 10^3; end

% create xyz matrix and split into cell chunks
xyz = [X(:), Y(:), Z(:)];
xyz = chunk2cell(xyz, chunksiz);

% change mesh format
if isfield(i_mesh, 'elem')
    i_mesh = flip_surf_format(i_mesh);
end

parfor i = 1:numel(xyz)
    pts_inside{i} = inpolyhedron(i_mesh, xyz{i});
end

pts_inside = cell2mat(pts_inside');
clear i xyz i_mesh

idx_inside = find(pts_inside == 1)';

xyz = [X(:), Y(:), Z(:)];
xyz = xyz(idx_inside, :);

end