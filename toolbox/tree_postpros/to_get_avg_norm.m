function [sel_idx, sel_coord, sel_point, sel_norm] = ...
    to_get_avg_norm(float_trees, ref_point, ref_normal, ...
    intersect_idx, intersect_coord, max_dist, i_buffer, neigh_dist)
% to_get_avg_norm: get average norm across a tree (from a population of trees and a initial reference tree)
% Usage:
%   to_get_avg_norm(float_trees, ref_point, ref_normal, ...
%      intersect_idx, intersect_coord, max_dist, i_buffer, neigh_dist)
%
% Args:
%   float_trees: floating trees
%   ref_point: reference points
%   ref_normal: reference normals
%   intersect_idx: indeces of intersection
%   intersect_coord: coordinates of intersection
%   max_dist: maximun distance from centroid
%   i_buffer: buffer
%   neigh_dist: distance from neighbor
%
% Output
%   sel_idx: index of points in trees close to intersection with normal plane
%   sel_coord: coordinates of points in trees intersecting normal plane
%   sel_point: points from average trace
%   sel_norm: normal from average trace

if ~exist('i_buffer', 'var') || isempty(i_buffer); i_buffer = 5; end
float_trees = getObj_xyz(float_trees);

% tag with nans all points lower/higher or equal to i_buffer and end - i_buffer

intersect_idx(intersect_idx <= i_buffer) = nan;
max_ = max(intersect_idx, [], 1) - i_buffer;
intersect_idx(intersect_idx - max_ >= 0) = nan;
 
clear max_

% find nan in rows
idx2del = sum(isnan(intersect_idx), 2) > 0;

intersect_idx(idx2del, :) = [];
intersect_coord(idx2del, :, :) = [];
ref_point(idx2del, :) = [];
ref_normal(idx2del, :) = [];

% get avg normals from intercepts
sel_point = nan(size(intersect_idx, 1), 3);
sel_idx = nan([size(intersect_idx, 1), numel(float_trees)]);
sel_coord = nan([size(intersect_idx, 1), 3, numel(float_trees)]);

for i = 1:size(intersect_idx, 1)

    fprintf('*')

    % get array of nearby points
    sel_pnts = intersect_idx(i, :);
    sel_pnts_m = [];
    
    try
        for j = 1:size(intersect_idx, 2)
            sel_pnts_m(:, :, j) = ...
                float_trees{j}(sel_pnts(j) - i_buffer: sel_pnts(j) + i_buffer, :);
        end
    catch
        keyboard
    end
    
    % replace intersect point by actual intersection
    sel_pnts_m(i_buffer + 1, :, :) = intersect_coord(i, :, :);
    sel_points_d  = squeeze((sum(bsxfun(@minus, ...
        sel_pnts_m, sel_pnts_m(i_buffer + 1, :, :)).^2, 2)).^0.5) <= neigh_dist;
    
    for j = 1:size(intersect_idx, 2)
        sel_norm(j, :) = sel_pnts_m(find(sel_points_d(:, j), 1, 'last'), :, j) ...
            - sel_pnts_m(find(sel_points_d(:, j), 1, 'first'), :, j);
    end

    % normalize length
    sel_norm = bsxfun(@rdivide, sel_norm, (sum(sel_norm.^2, 2)).^0.5);
    
    % avg normal
    sel_norm_avg = sum(sel_norm, 1);
    sel_norm(i, :) = sel_norm_avg;
    sel_point_avg = sum(squeeze(sel_pnts_m(i_buffer + 1, :, :))', 1)/size(intersect_idx, 2);
    sel_point(i, :) = sel_point_avg; % centroid

    [sel_idx(i, :), sel_coord(i, :, :)] = ...
        to_get_planeintersects(float_trees, sel_point_avg, sel_norm_avg, max_dist);

    % plot new avg plane
    % if i == 10
    % 
    %     i_delta = 10; i_sel = 1:2; r_ = 20; xy_range = [-r_ r_ -r_ r_];
    %     plotorthoplanes_3D(obj, [], [sel_point_avg; ref_point(i, :)], ...
    %         [sel_norm_avg; ref_normal(i, :)], xy_range, i_delta, imatrix, i_sel);
    %     [figH, axH, pxH] = plotorthoplanes_3D(obj, f2use, i_point, i_normal, xy_range, i_delta, imatrix, i_sel, p_colorvect, t_colorvect)
    % 
    % end

end

fprintf('\n')

end
        