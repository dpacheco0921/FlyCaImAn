function [i_point, i_normal] = to_get_point_normal(iObject, norm_step)
% to_get_point_normal: walk through the tree and generates point and normal vector
% to generate orthogonal planes to each xyz point
% Usage:
%   to_get_point_normal(iObject, i_xyz)
%
% Args:
%   iObject: tree obj / xyz matrix
%   norm_step: number of points beyond to use for normal stimation
%
% Notes: see 'treestoolbox' for format info

% get_point_normal(obj, f2use, norm_step, imatrix)
% walk through the tree and generates point and normal vector
% to generate orthogonal planes to each xyz point

if isempty(norm_step)
    norm_step = 5;
end

i_normal = [];
i_point = [];

xyz = getObj_xyz(iObject);

% plane equation: (p - i_point).i_normal = 0
for i = 1:(size(xyz, 1) - norm_step)
    
    i_point(i, :) = xyz(i, :);
    i_normal(i, :) = xyz(i, :) - xyz(i + norm_step, :);
    
end
    
end