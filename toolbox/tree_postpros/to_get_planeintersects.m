function [i_idx, i_coord] = to_get_planeintersects(iObject, ref_point, ref_normal, max_dist, verbose)
% to_get_planeintersects: get line plane intersections close to reference point (i_point(i))
% p = dl + lo; line & plane intersection
% d = (i_point - lo).i_normal/l.i_normal
%
% Usage:
%   to_get_planeintersects(iObject)
%
% Args:
%   iObject: tree obj / xyz matrix
%   ref_point: reference point
%   ref_normal: reference normal
%   max_dist: max distance of points to test in a.u units (um)
%   verbose: verbose gate
%
% Notes: see traceObj

if ~exist('verbose', 'var'); verbose = 0; end

i_idx = nan([size(ref_point, 1), numel(iObject)]);
i_coord = nan([size(ref_point, 1), 3, numel(iObject)]);

if verbose; fprintf(['Running ', num2str(numel(iObject)), ' trees \n']); end

for i_tree = 1:numel(iObject)
    
    float_points = getObj_xyz(iObject(i_tree));
    
    i_idx_pp = nan(1, size(ref_point, 1));
    i_coord_pp = cell(1, size(ref_point, 1));
    
    parfor i_ref_p = 1:size(ref_point, 1)
        
        if verbose; fprintf('*'); end
        
        % compute distance of points to ref points
        dist_to_refpoint = to_dist_to_refpoimt(float_points, ref_point(i_ref_p, :));
        k_idx = find(dist_to_refpoint < max_dist)';
        
        % find intersections (plane to line segments)
        pre_i_coord = to_get_plane2segline_intersect(k_idx, float_points, ref_point(i_ref_p, :), ref_normal(i_ref_p, :));
        
        % chose the set of points with the smallest distance
        if sum(~isnan(sum(pre_i_coord, 2))) ~= 0
            
            k_idx2sel = k_idx(~isnan(sum(pre_i_coord, 2)));
            [~, k_idx_i] = sort(dist_to_refpoint(k_idx2sel), 'ascend');
            i_idx_pp(i_ref_p) = k_idx2sel(k_idx_i(1));
            i_coord_pp{i_ref_p} = pre_i_coord(k_idx == i_idx_pp(i_ref_p), :);
            
        else
            
            i_idx_pp(i_ref_p) = nan;
            i_coord_pp{i_ref_p} = [nan nan nan];
            
        end
       
    end
    
    i_idx(:, i_tree) = i_idx_pp';
    i_coord(:, 1:3, i_tree) = cell2mat(i_coord_pp');
    
    if verbose; fprintf('\n'); end
    
end

end