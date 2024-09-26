function i_coord = to_get_plane2segline_intersect(f_points2test, float_points, ref_point, ref_normal)
% to_get_plane2segline_intersect: find intersection of point(s) &
% plane(s) to segments of lines (points contained within those segments)
% Usage:
%   to_get_plane2segline_intersect(f_points2test, float_points, ref_point, ref_normal)
%
% Args:
%   f_points2test: indeces to test
%   ref_point: reference point ([1 by xyz])
%   float_points: floating points ([n by xyz])
%   ref_normal: reference normal
%
% Notes: see traceObj

i = 1;
i_coord = [];
f_points_n = size(float_points, 1);
i_dist = [];
p_dist = [];

for k = f_points2test
    
    % for each point the next or previous point to calculate the line equation
    k_points = [k k + 1];
    if k == f_points_n
        k_points = [k k - 1];
    end
    
    % find plane & line intersection 'i_p'
    l = float_points(k_points(1), :) - float_points(k_points(2), :);
    d = (ref_point - float_points(k_points(1), :))*ref_normal'/((l)*ref_normal');
    i_p = float_points(k_points(1), :) + d*l;
    
    % if i_p = Inf line is parallel, if zero: ref point is contained in plane and line
    % determine if 'i_p' is between the points used to calculate it
    i_dist(i, 1) = (sum(l.^2))^0.5;
    p_dist(i, :) = [(sum((float_points(k_points(1), :) - i_p).^2))^0.5, ...
        (sum((float_points(k_points(2), :) - i_p).^2))^0.5];
    
    if round(i_dist(i, 1)*10^4) ~= round(sum(p_dist(i, :))*10^4)
        i_p = [nan nan nan];
    end
    
    % if p = nan, then the intersection is outside line interval
    i_coord(i, :) = i_p;
    i = i + 1;
    
    clear close_points l d p
    
end

end