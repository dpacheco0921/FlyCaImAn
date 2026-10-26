function [axis_jitter, distance_jitter] = to_get_xyz_jitter(i_coord, i_centroid)
% to_get_xyz_jitter: calculate jitter by measuring distance to centroid per axis and euclidean distance
% Usage:
%   [axis_jitter, distance_jitter] = to_get_xyz_jitter(i_coord, i_centroid)
%
% Args:
%   i_coord: [n by 3 by trees] input set of coordinates
%   i_centroid: [n by 3] reference centroid (default, mean across trees;
%       any NaN in a point makes that centroid NaN)
%
% Returns:
%   axis_jitter: [n by 3 by trees] signed X, Y, Z offset from the centroid
%       (same units as the coordinates)
%   distance_jitter: [n by trees] euclidean distance to the centroid

if ~exist('i_centroid', 'var') || isempty(i_centroid); i_centroid = []; end
    
% get centroid
if isempty(i_centroid)
    i_centroid = sum(i_coord, 3)/size(i_coord, 3); 
end

% get jitter on each axis separately
axis_jitter = bsxfun(@minus, i_coord, i_centroid);

% get jitter on distance;
distance_jitter = bsxfun(@minus, i_coord, i_centroid);
distance_jitter = squeeze(sum(distance_jitter.^2, 2).^0.5);

end