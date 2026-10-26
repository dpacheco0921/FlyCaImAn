function dist_to_refpoint = to_dist_to_refpoimt(float_points, ref_point)
% to_dist_to_refpoimt: calculate euclidean distance of a matrix of points [n, 3]  (n, [xyz]) to a reference point
% Usage:
%   dist_to_refpoint = to_dist_to_refpoimt(float_points, ref_point)
%
% Args:
%   float_points: floating points ([n by xyz])
%   ref_point: reference point ([1 by xyz])
%
% Returns:
%   dist_to_refpoint: [n by 1] distance of each point to ref_point (same
%       units as the coordinates, um)
% Notes: see traceObj

dist_to_refpoint = bsxfun(@minus, float_points, ref_point);
dist_to_refpoint = (sum(dist_to_refpoint.^2, 2)).^0.5;

end