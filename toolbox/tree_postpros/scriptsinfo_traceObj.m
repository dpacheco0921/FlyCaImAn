%% Tools & scripts to edit and perfom vasic operations with trees | neurons (swc files)
%% 1) segment neurons
% use fiji simple neuron trace and save in swc format

%% 2) Using treestoolbox to load and save trees in matlab
%% 2.1) load swc files
edit load_tree

%% 2.2) save swc files
edit swc_tree

%% 3) Using CMTK & traceObj to reformat traces in matlab
% make object
obj = traceObj;
fly2nsybIVA_xform_perfile(obj, floatTrees, refIm, xDir, iDir, oDir, eSuffix, redo);
% help to_fly2nsybIVA_xform_perfile

% examples
edit neuron_dsx
edit neuron_G6S_GC
edit neuron_NP2631FruFLP

%% 4) Using traceObj to vizualize and perform basic operations with traces in matlab
%% 4.1) make object and load trees

obj = traceObj; % make object
loadtrees(obj, [], 'P1_nsybIVAi') % load trees

%% 4.2) extract and edit XYZ matrix

xyzmatrix(obj, []) % extract XYZ matrix

smoothtree(obj, 5) % smooth XYZ matrix

pruneNaNs(obj) % prune trees with consecutive Nans

%% 4.3) operation with trees (see plots in 4.4)

f2use = 1; norm_step = 5; imatrix = 1; max_dist = 40; verbose = 0; i_buffer = 5; neigh_dist = 3;

% (these steps now use the to_* functions, which do not depend on obj;
%  pick the XYZ field with to_pick_xyz_object(obj, idx, imatrix))
% i) get intersecting planes to selected tree
iObj_ref = to_pick_xyz_object(obj, f2use, imatrix);
iObj_float = to_pick_xyz_object(obj, [], imatrix);
[ref_point, ref_normal] = to_get_point_normal(iObj_ref, norm_step);

% to_get_point_normal(iObject, norm_step)
% ii) get coordinates and idx of points where planes intersect each treeS
[i_idx, i_coord] = to_get_planeintersects(iObj_float, ref_point, ref_normal, max_dist, verbose);

% to_get_planeintersects(iObject, ref_point, ref_normal, max_dist, verbose)
% compute avg norm and recalculate intersections
[sel_idx, sel_coord, sel_point, sel_norm] = to_get_avg_norm(iObj_float, ref_point, ref_normal, i_idx, i_coord, max_dist, i_buffer, neigh_dist);

% to_get_avg_norm(float_trees, ref_point, ref_normal, intersect_idx, intersect_coord, max_dist, i_buffer, neigh_dist)
% compute jitter
[axis_jitter, distance_jitter] = to_get_xyz_jitter(sel_coord);

% to_get_xyz_jitter(i_coord, i_centroid)
% batch 
[axis_jitter, distance_jitter] = ...
                get_jitter_from_treegroup(obj, f2use, norm_step, max_dist, neigh_dist, i_buffer, imatrix, verbose);
            
%% 4.4) plotting trees
% plot group of trees
plottree_2D(obj, [], 'XY')
plottree_2D(obj, [], 'YZ')
plottree_2D(obj, [], 'XZ')
plottree_3D(obj, []); view([0 90])
plottree_3D(obj, [], 1)

% ** i) plot indeces of intersection points
figure(); imagesc(i_idx)

% ** ii) plot tree and normal planes 2D
xy_range = [240 300]; i_delta = 10; f2use = 1; imatrix = 1;
plotorthoplanes_2D(obj, f2use, 'XY', ref_point, ref_normal, xy_range, i_delta);

% plotorthoplanes_2D(obj, f2use, axis2proj, i_point, i_normal, x_range, i_delta)
% ** iii) plot tree and normal planes 3D
i_delta = 10; i_sel = [37 50 100]; r_ = 20;
xy_range = [-r_ r_ -r_ r_];
[figH, axH, pxH] = plotorthoplanes_3D(obj, [], ref_point, ref_normal, xy_range, i_delta, imatrix, i_sel, [0 0 0], [0 1 0]);

% plotorthoplanes_3D(obj, f2use, i_point, i_normal, xy_range, i_delta, imatrix, i_sel, p_colorvect, t_colorvect)
% ** iv) plot tree points closest to intersection
i_xyz = to_plot_point_3D(obj, i_idx(i_sel, :), 1);

% ** v) plot tree-plane intersections
i_coord2plot = squeeze(i_coord(i_sel, :, :))';
scatter3(i_coord2plot(:, 1), i_coord2plot(:, 2), i_coord2plot(:, 3), '.', 'MarkerEdgeColor', [1 0 1])

% ** vi) plot reference tree-point
scatter3(obj.xyz_s{1}(i_sel, 1), obj.xyz_s{1}(i_sel, 2), obj.xyz_s{1}(i_sel, 3), 'o', 'MarkerEdgeColor', [1 1 0])

% ** vii) plot jitter
figHandle = figure('Name', 'jitter'); axHandle = subplot(1, 1, 1);
hbins = -30:4:30; colorvect = [1 0 0; 0 0 1; 0 0 0];
for j = 1:3 % XYZ axis
   vect = axis_jitter(:, j, :); 
   vect = vect(:); vect(isnan(vect)) = [];
   [y, ~] = hist(vect, hbins);
   plot(hbins, (y/sum(y))*100, 'Color', colorvect(j, :), ...
       'Linewidth', 2, 'Parent', axHandle); hold(axHandle, 'on')
end
axHandle.XTickMode = 'manual';
axHandle.XTick = [-25:5:25];
xlim(axHandle, [-25 25]);
figEdit(axHandle, figHandle, 'w')
set(axHandle, 'FontSize', 15);
