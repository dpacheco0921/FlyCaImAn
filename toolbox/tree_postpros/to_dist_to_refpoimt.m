function dist_to_refpoint = to_dist_to_refpoimt(float_points, ref_point)
% to_get_planeintersects: calculate euclidean distance of a matrix of points [n, 3]  (n, [xyz]) to a reference point
% Usage:
%   to_dist_to_refpoimt(float_points, ref_point)
%
% Args:
%   float_points: floating points ([n by xyz])
%   ref_point: reference point ([1 by xyz])
% Notes: see traceObj

dist_to_refpoint = bsxfun(@minus, float_points, ref_point);
dist_to_refpoint = (sum(dist_to_refpoint.^2, 2)).^0.5;

end