classdef traceObj < handle
    % object to perform basic geometric operations on traces (trees)
   
   properties
        iDir;         % input file related (spattern or ssuffix)
        traces;       % trees to work with
        xyz;          % matrix of XYZ axis
        xyz_s;        % smooth XYZ axis
        xyz_pix;      % XYZ matrix but in pixel coordinates
        atlas;        % info of atlases (required for reformatting/mirroring)
   end
   
   methods (Access = 'public')
       
        % **************************** build object ****************************
        
        function obj = traceObj
            % traceObj: initialize traceObj (creates: iDir, traces, atlas)
            
            % initialize
            obj.iDir = pwd; 
            obj.traces = struct([]);
            
            if exist([strrep(which('traceObj'), ...
                    'traceObj.m', ''), ...
                    'atlas_meta.mat'], 'file')
                
                load([strrep(which('traceObj'), ...
                    'traceObj.m', ''), ...
                    'atlas_meta.mat'])
                obj.atlas = atlas;
                
            end
            
        end
        
        % **************************** build objects ****************************
        
        function loadtrees(obj, iDir, iSuffix, FileName, catgate, verbose)
            % loadtrees(obj, iDir, iSuffix) load trees from swc files
            
            if exist('iDir', 'var') && ~isempty(iDir)
                obj.iDir = iDir;
            end
            
            if ~exist('iSuffix', 'var')
                
                iSuffix = [];
                
            else
                
                if ~isempty(iSuffix) && ~contains(iSuffix(1), '*')
                    iSuffix = ['*', iSuffix]; 
                end
                
            end
            
            if ~exist('FileName', 'var')
                FileName = [];
            end
            
            if ~exist('catgate', 'var') || isempty(catgate)
                catgate = 0;
            end
                        
            if ~exist('verbose', 'var') || isempty(catgate)
                verbose = 0;
            end
            
            f2load = rdir([obj.iDir, filesep, iSuffix, '*.swc']);
            f2load = str2match(FileName, f2load);
            f2load = {f2load.name}';
            
            if catgate
                idx = numel(obj.traces) + (1:numel(f2load));
            else
                obj.traces = struct([]);
                idx = 1:numel(f2load);
            end
            
            k = 1;
            for i = idx
                
                try
                    
                    % flag to not repair
                    option = 'no';
                    pretraces = load_tree(f2load{idx == i}, option);
                    
                    % if tree is composed of multiple                    
                    if numel(pretraces) > 1
                        keyboard
                    end
                    
                    % only perform this repairs
                    % eliminate trifurcations by adding short segments:
                    pretraces = elimt_tree (pretraces);        
                    % eliminate 0-length compartments:
                    pretraces = elim0_tree (pretraces); 
                    
                    if i == 1 && ~catgate
                        obj.traces = struct(pretraces);
                    end
                    obj.traces(k, 1) = pretraces;
                    k = k + 1;
                    
                    % plot loaded tree
                    if verbose
                        to_plottree_3d(getObj_xyz(pretraces), ...
                            idpar_tree(pretraces), [], 1, [0 1 1])
                    end
                    
                catch
                    
                    fprintf(['skeleton ', f2load{idx == i}, ' failed'])
                    
                end
                
            end
            
        end
        
        function collect_atlas_info(obj, atlasname, aDir, deviceId)
            % collect_atlas_info: code to collect atlas information and save it in a mat file
            % Note:
            %   see to_collect_atlas_info

            if ~exist('deviceId', 'var') || ...
                    isempty(deviceId)
                deviceId = 1;
            end
            
            if ~exist('aDir', 'var') || ...
                    isempty(aDir)
                aDir = customdirs_deafult(deviceId);
            end
            
            if ~exist('atlasname', 'var') || ...
                    isempty(atlasname)
                atlasname = {'JFRC2', 'VNCIS1', ...
                    'IS2', 'FCWB', 'IBNWB', 'nsybIVAi'};
            end

            to_collect_atlas_info(obj, atlasname, aDir, deviceId)
            
        end
        
        % **************************** batch operations ****************************
        
        function [axis_jitter, distance_jitter, avg_point, agv_norm] = ...
                get_jitter_from_treegroup(obj, reftree, ...
                norm_step, max_dist, neigh_dist, i_buffer, ...
                imatrix, verbose, floattrees_idx)
            % [axis_jitter, distance_jitter] = ...
            %    get_jitter_from_treegroup(obj, reftree, ...
            %    norm_step, max_dist, neigh_dist, i_buffer, ...
            %   imatrix, verbose, floattrees_idx)
            % function performs batch jobs (many functions in one) to get XYZ jitter
            
            if ~exist('imatrix', 'var') || isempty(imatrix); imatrix = 0; end
            if ~exist('floattrees_idx', 'var') || isempty(floattrees_idx)
                floattrees_idx = 1:numel(obj.xyz);
            end
            if ~exist('reftree', 'var') || isempty(reftree); reftree = 1; end
            
            % define floating trees
            iObj_float = to_pick_xyz_object(obj, floattrees_idx, imatrix);

            % define reference tree
            iObj_ref = to_pick_xyz_object(obj, reftree, imatrix);
                 
            % 1) get reference points and normals
            iTic = tic;
            
            fprintf('Calculate points and normals\n')
            
            tic
            [ref_point, ref_normal] = to_get_point_normal(iObj_ref, norm_step);
            % help to_get_point_normal
            toc
                        
            % 2) get coordinates and idx of points where planes intersect each tree
            fprintf('Get intersections\n')
            
            tic
            [i_idx, i_coord] = ...
                to_get_planeintersects(iObj_float, ref_point, ref_normal, max_dist, verbose);
            % help to_get_planeintersects
            toc
            
            % 3) compute avg norm and recalculate intersections
            fprintf('Get avg normal and get intersections\n')
            
            tic;
            [~, sel_coord, avg_point, agv_norm] = to_get_avg_norm(iObj_float, ...
                ref_point, ref_normal, i_idx, i_coord, max_dist, i_buffer, neigh_dist);
            toc
            
            % 4) compute jitter (remove points where there is no overlap with all traces)
            fprintf('Get jitter\n')
            points2del = ~sum(sum(isnan(sel_coord), 3), 2);
            fprintf(['Deleting ', num2str(sum(~points2del)), ...
                ' points (have nans) out of ', num2str(numel(points2del)), '\n'])
            
            [axis_jitter, distance_jitter] = ...
                to_get_xyz_jitter(sel_coord(points2del, :, :), avg_point(points2del, :));
            
            toc(iTic)
            
        end
        
        % **************************** reformatting files using CMTK tools ****************************
        
        function fly2nsybIVA_xform_perfile(obj, ...
                floatTrees, refIm, xDir, iDir, ...
                oDir, esuffix, redo)
            % function that performs reformatting of trees from fly to
            % nsybIVAi reference coordinates using CMTK  and saves new SWC files
            % Note:
            %   see to_fly2nsybIVA_xform_perfile
            
            if ~exist('esuffix', 'var') || ...
                    isempty(esuffix)
                esuffix = '';
            end
            
            if ~exist('xDir', 'var') || ...
                    isempty(xDir)
                xDir = '.';
            end
            
            if ~exist('iDir', 'var') || ...
                    isempty(iDir)
                iDir = '.';
            end
            
            if ~exist('oDir', 'var') || ...
                    isempty(oDir)
                oDir = '.';
            end
            
            if ~exist('redo', 'var') || ...
                    isempty(redo)
                redo = 0;
            end

            to_fly2nsybIVA_xform_perfile(...
                floatTrees, refIm, xDir, ...
                iDir, oDir, esuffix, redo)
            
        end

        function atlas2tlas_xform_perfile(obj, ...
                floatTrees, refi, refo, xDir, ...
                iDir, oDir, esuffix, redo)
            % function that performs reformatting of trees from atlas to
            % atlas coordinates using CMTK and saves new SWC files
            % Note:
            % see to_atlas2atlas_xform_perfile
            
            if ~exist('esuffix', 'var') || isempty(esuffix)
                esuffix = '';
            end
            
            if ~exist('xDir', 'var') || isempty(xDir)
                xDir = '.';
            end
            
            if ~exist('iDir', 'var') || isempty(iDir)
                iDir = '.';
            end
            
            if ~exist('oDir', 'var') || isempty(oDir)
                oDir = '.';
            end
            
            if ~exist('redo', 'var') || isempty(redo)
                redo = 0;
            end

            to_atlas2atlas_xform_perfile(floatTrees, ...
                refi, refo, xDir, iDir, oDir, esuffix, redo)
            
        end
        
        function atlas2tlas_xform(obj, refi, refo, xDir, esuffix)
            % function that performs reformatting of trees from atlas to atlas coordinates using CMTK
            % Note:
            % see to_atlas2atlas_xform
            
            if ~exist('esuffix', 'var') || ...
                    isempty(esuffix)
                esuffix = '';
            end
            
            if ~exist('xDir', 'var') || ...
                    isempty(xDir)
                xDir = '.';
            end

            for i = 1:numel(obj.trace)
                [obj.traces(i), ~] = ...
                    to_atlas2atlas_xform(obj.traces(i), refi, refo, xDir, esuffix);
            end
            
        end        
        
        % **************************** internal tools ****************************
        
        function xyzmatrix(obj, f2use)
            % xyzmatrix: get XYZ coordinates from trees (stored as a cell)
            % Usage:
            %   xyzmatrix(obj, f2use)
            %
            % Args:
            %   obj: traceObj obj
            %   f2use: files to use
            %
            % Note:
            % see getObj_xyz
            
            if ~exist('f2use', 'var') || isempty(f2use)
                f2use = 1:numel(obj.traces);
            end
            
            obj.xyz = getObj_xyz(obj.traces(f2use));

        end
        
        function pruneNaNs(obj, imatrix)
            % pruneNaNs: fix single NaN points, deletes trees with more than 1 consecutive NaNs
            % Usage:
            %   pruneNaNs(obj, resolution_)
            %
            % Args:
            %   obj: traceObj obj
            %   imatrix: field to use
            %
            % Note:
            % see to_get_trace_with_nans to_pick_xyz_object getObj_xyz
            
            if ~exist('imatrix', 'var'); imatrix = 0; end
            
            % find traces with nan and fix the ones with single cases
            % it updates obj.xyz
            [trace2keep_idx, obj.xyz] = to_get_trace_with_nans(obj, imatrix);
           
            fprintf(['Trees deleted ', num2str(sum(~trace2keep_idx)), ...
                'out of ', num2str(length(trace2keep_idx))])
            
            display([{obj.traces.name}', num2cell(~trace2keep_idx)])
            
            obj.traces = obj.traces(trace2keep_idx);
            obj.xyz = obj.xyz(trace2keep_idx);
            
        end
        
        function path2root(obj, f2use)
            % path2root: find the longest path to root and updates trees
            %
            % Usage:
            %   path2root(obj, resolution_)
            %
            % Args:
            %   obj: traceObj obj
            %   resolution_: new resolution in um
            %
            % Note:
            % see resample_tree
            
            if ~exist('f2use', 'var') || ...
                    isempty(f2use)
                f2use = 1:numel(obj.traces);
            end
            
            for i = f2use
                
                paths2root = ipar_tree(obj.traces(i));
                path_idx = sum(paths2root > 0, 2) == ...
                    max(sum(paths2root > 0, 2));
                nodes = sort(paths2root(path_idx, :));
                nodes(nodes == 0) = [];
                
                [~, obj.traces(i)] = sub_tree(obj.traces(i), nodes);
                
                clear nodes path_idx paths2root
                
            end            
            
        end
        
        function resampletree(obj, resolution_)
            % resampletree: resample trees along spatial axis
            % Usage:
            %   resampletree(obj, resolution_)
            %
            % Args:
            %   obj: traceObj obj
            %   resolution_: new resolution in um
            %
            % Note:
            % see resample_tree
            
            for i = 1:numel(obj.traces)
                
                obj.traces(i) = resample_tree(obj.traces(i), resolution_);
                fprintf(['finished skeleton ', obj.traces(i).name, '\n'])
                
            end
            
        end
        
        function iObject = mirror_brain(obj, iObject, axis2flip, atlasname)
            % mirror_brain: mirror iObject (xyz coordinates within it)
            % Note:
            % see to_mirror_axis
            atlas_idx = contains({obj.atlas.name}, atlasname);
            boundingbox = obj.atlas(atlas_idx).boundingbox;
            
            iObject = to_mirror_axis(iObject, axis2flip, boundingbox);
            % help to_mirror
        end
        
        function smoothtree(obj, span, imatrix)
            % smoothtree: smooth trees over space
            % Usage:
            %   smoothtree(obj, span)
            %
            % Args:
            %   obj: traceObj obj
            %   span: span defines the window
            %
            % Note:
            % see to_pick_xyz_object getObj_xyz

            if ~exist('imatrix', 'var'); imatrix = 0; end

            iObject = to_pick_xyz_object(obj, [], imatrix);
            
            if ~iscell(getObj_xyz(iObject))
                xyz_int{1} = getObj_xyz(iObject);
            else
                xyz_int = getObj_xyz(iObject);
            end

            for i = 1:numel(obj.xyz)
                
                obj.xyz_s{i, 1} = [smooth(xyz_int{i}(:, 1), span), ...
                    smooth(xyz_int{i}(:, 2), span), ...
                    smooth(xyz_int{i}(:, 3), span)];
                
            end
            
        end
        
        function xyz_um2pix(obj, imatrix, iscale)
            % xyz_um2pix: convert um to pixel values
            % Usage:
            %   xyz_um2pix(obj, imatrix, iscale)
            %
            % Args:
            %   obj: traceObj obj
            %   imatrix: field to use
            %   iscale: indeces of parents to xyz
            %
            % Note: 
            % see to_pick_xyz_object getObj_xyz
            
            if ~exist('imatrix', 'var'); imatrix = 0; end

            iObject = to_pick_xyz_object(obj, [], imatrix);
            
            if ~iscell(getObj_xyz(iObject))
                xyz_int{1} = getObj_xyz(iObject);
            else
                xyz_int = getObj_xyz(iObject);
            end

            for i = 1:numel(obj.xyz)
            	
                obj.xyz_pix{i, 1} = bsxfun(@rdivide, xyz_int{i, 1}, iscale);
                
            end
            
        end
        
        % **************************** plotting ****************************
        
        function [figH, axH] = plottree_2D(obj, f2use, axis2proj, ...
                axH, imatrix, verbose, color2use, tstyle, psoma)
            % plottree_2D: plot traces in 2D format
            % Usage:
            %   plottree_2D(obj, f2use, axis2proj, ...
            %       axH, imatrix, verbose, color2use, tstyle, psoma)
            % Note: 
            % see to_plottree_2d
                      
            if ~exist('verbose', 'var'); verbose = 0; end
            if isempty(axis2proj); axis2proj = 'XY'; end
            if ~exist('axH', 'var') || isempty(axH)
                figH = figure(); axH = subplot(1, 1, 1);
            end
            if isempty(f2use); f2use = 1:numel(obj.traces); end
            if ~exist('imatrix', 'var') || isempty(imatrix); imatrix = 0; end
            if ~exist('color2use', 'var') || isempty(color2use)
                color2use = [0 0 0]; 
                if numel(f2use) ~= 1
                    color2use = gradientgen(7, numel(f2use));
                end
            else
                if size(color2use, 1) == 1
                    color2use = repmat(color2use, [numel(f2use) 1]);
                end
            end
            if ~exist('tstyle', 'var') || isempty(tstyle); tstyle = '-'; end
            if ~exist('psoma', 'var') || isempty(psoma); psoma = 1; end

            % define xzy values to use
            iObject = to_pick_xyz_object(obj, [], imatrix);
            
            % plot projections
            for i = f2use
                
                idpar = idpar_tree(obj.traces(i));
                
                to_plottree_2d(getObj_xyz(iObject(i)), idpar, ...
                    axis2proj, axH, verbose, color2use(i, :), tstyle, ...
                    psoma)

            end
            
        end
        
        function [figH, axH, linH, somaH] = plottree_3D(obj, f2use, axH, ...
                imatrix, verbose, color2use, tstyle, psoma, tagflag)
            % plottree_3D: plot traces in 3D format
            % Usage:
            %   plottree_3D(obj, f2use, axH, ...
            %      imatrix, verbose, color2use, tstyle, psoma, tagflag)
            %
            % Note:
            % see to_plottree_3d

            if ~exist('verbose', 'var')
                verbose = 0;
            end
            
            if ~exist('axH', 'var') || isempty(axH)
                figH = figure();
                axH = subplot(1, 1, 1);
            end
            
            if isempty(f2use)
                f2use = 1:numel(obj.traces);
            end
            
            if ~exist('imatrix', 'var') || ...
                    isempty(imatrix)
                imatrix = 0;
            end
            
            if ~exist('color2use', 'var') || ...
                    isempty(color2use)
                
                color2use = [0 0 0]; 
                if numel(f2use) ~= 1
                    color2use = gradientgen(7, numel(f2use));
                end
                
            else
                
                if size(color2use, 1) == 1
                    color2use = repmat(color2use, ...
                        [numel(f2use) 1]);
                end
                
            end
            
            if ~exist('tstyle', 'var') || ...
                    isempty(tstyle)
                tstyle = '-';
            end
            
            if ~exist('psoma', 'var') || ...
                    isempty(psoma)
                psoma = 1;
            end

            if ~exist('tagflag', 'var')
                tagflag = 0;
            end
            
            % define xzy values to use
            iObject = to_pick_xyz_object(obj, [], imatrix);
            
            % check if there exist tags
            if isfield(obj.traces(1), 'tagnodes') && tagflag
                tags = {obj.traces.tagnodes};
            else
                tags = cell(numel(obj.traces), 1);
            end
            
            % plot projections
            k = 1;
            
            for i = f2use
                
                idpar = idpar_tree(obj.traces(i));
                
                if isstruct(iObject) || iscell(iObject)
                    
                    [linH(i), somaH(i)] = ...
                        to_plottree_3d(getObj_xyz(iObject(i)), ...
                        idpar, axH, verbose, color2use(k, :), ...
                        tstyle, psoma, tags{i});
                    
                else
                    
                    [linH(i), somaH(i)] = ...
                        to_plottree_3d(getObj_xyz(iObject), ...
                        idpar, axH, verbose, color2use(k, :), ...
                        tstyle, psoma, tags{i});
                   
                end
                
                k = k + 1;
                
            end
            
            if ~exist('figH', 'var')
                figH = gcf;
            end
            
        end
        
        function plotorthoplanes_2D(obj, f2use, axis2proj, i_point, i_normal, x_range, i_delta)
            % plotorthoplanes_2D(obj, f2use, axis2proj, i_point, i_normal, x_range, i_delta)
            % plane equation: (p - i_point).i_normal = 0 || ax + by + cz + d = 0
            if isempty(i_delta); i_delta = 1; end
            figH = figure(); axH = subplot(1, 1, 1);
            a = i_normal(:, 1); b = i_normal(:, 2); c = i_normal(:, 3);
            colorvect = gradientgen(7, numel(a));
            % plot projections
            for i = 1:i_delta:numel(a)
                d = -i_point(i, :)*i_normal(i, :)';
                switch axis2proj
                    case 'XY'
                        y = @(x) (-c(i)*i_point(i, 3)-d-a(i)*x)/b(i);
                    case 'YZ'
                        y = @(x) (-a(i)*i_point(i, 1)-d-b(i)*x)/c(i);
                    case 'XZ'
                        y = @(x) (-b(i)*i_point(i, 2)-d-a(i)*x)/c(i);
                end
                fplot(y, x_range, 'Color', colorvect(i, :))
                if i == 1; hold on; end
            end
            plottree_2D(obj, f2use, axis2proj, axH)
            % edit labels
            switch axis2proj
                case 'XY'
                    xlabel(axH, 'X'); ylabel(axH, 'Y');
                case 'YZ'
                    xlabel(axH, 'Y'); ylabel(axH, 'Z');
                case 'XZ'
                    xlabel(axH, 'X'); ylabel(axH, 'Z');
            end
        end
        
        function [figH, axH, pxH] = plotorthoplanes_3D(obj, f2use, i_point, i_normal, xy_range, i_delta, imatrix, i_sel, p_colorvect, t_colorvect)
            % plotorthoplanes_3D(obj, f2use, i_point, i_normal, xy_range, i_delta, imatrix, i_sel, p_colorvect, t_colorvect)
            % plane equation: (p - i_point).i_normal = 0 || ax + by + cz + d = 0
            
            siz = size(i_normal);
            if isempty(i_delta); i_delta = 1; end
            if ~exist('imatrix', 'var') || isempty(imatrix); imatrix = 0; end
            if ~exist('p_colorvect', 'var') || isempty(p_colorvect)
                p_colorvect = gradientgen(7, siz(1));
            else
                if size(p_colorvect, 1) == 1
                    p_colorvect = repmat(p_colorvect, [siz(1) 1]);
                end
            end
            if ~exist('t_colorvect', 'var'); t_colorvect = []; end
            if ~exist('i_sel', 'var') || isempty(i_sel)
                i_sel = 1:i_delta:siz(1);
            end
            
            figH = figure();
            axH = subplot(1, 1, 1);
            
            ii = 1;
            a = i_normal(:, 1);
            b = i_normal(:, 2);
            c = i_normal(:, 3);

            for i = i_sel
                d = -i_point(i, :)*i_normal(i, :)';
                z = @(x, y) (-a(i)*x - b(i)*y - d)/c(i);
                i_xy_range = xy_range + i_point(i, [1 1 2 2]);
                pxH(ii) = fsurf(z, i_xy_range, 'FaceAlpha', 0.3, ...
                    'EdgeColor', p_colorvect(i, :), 'FaceColor', p_colorvect(i, :));
                hold(axH, 'on'); ii = ii + 1;
            end
            
            plottree_3D(obj, f2use, axH, imatrix, 0, t_colorvect)
            
        end
        
   end
   
end
