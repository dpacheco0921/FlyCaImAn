function [linH, somaH] = to_plottree_3d(...
    iObject, idpar, axH, verbose, ...
    colormatrix, tstyle, psoma, nodetags)
% to_plottree_3d: plot traces in 3D format
%
% Usage:
%   to_plottree_3d(iObject, idpar, axH, verbose, ...
%      colormatrix, tstyle, psoma, nodetags)
%
% Args:
%   iObject: tree or xyz matrix
%   idpar: indeces of parents to xyz (see idpar_tree)
%   axH: axes handle
%   verbose: verbose
%   colormatrix: matrix with colors to use
%   tstyle: line style
%   psoma: plot soma
%   nodetags: nodes to tag (plot in different color - black)
%
% Todo:
% add flag to select some nodes
% make soma color independent of trace color

% default params

if ~exist('verbose', 'var') || ...
        isempty(verbose)
    verbose = 0;
end

if ~exist('colormatrix', 'var') || ...
        isempty(colormatrix)
    colormatrix = [0 0 0];
end

if ~exist('axH', 'var') || ...
        isempty(axH)
    figH = figure();
    axH = subplot(1, 1, 1);
end

if ~exist('tstyle', 'var') || ...
    isempty(tstyle)
    tstyle = '-';
end

if ~exist('psoma', 'var') || ...
        isempty(psoma)
    psoma = 1;
end

if ~exist('nodetags', 'var') || ...
        isempty(nodetags)
    nodetags = [];
end

if ~exist('idpar', 'var') || ...
        isempty(idpar)
    if isstruct(iObject) && isfield(iObject, 'dA')
        idpar = idpar_tree(iObject);
    else
        idpar = [];
    end
end

xyz = getObj_xyz(iObject);
tree_siz = size(xyz);

if verbose
    fprintf(['Number of NaNs ', ...
        num2str(sum(isnan(xyz(:)))), '\n']);
end

% make vector with all segments concatenated
if ~isempty(idpar)
    
    k = 1;

    for j = 1:tree_siz(1)

        t_xyz(:, k:k+1) = [xyz(idpar(j), :)', xyz(j, :)'];
        t_xyz(:, end + 1) = nan;
        k = size(t_xyz, 2) + 1;

    end

    linH = plot3(t_xyz(1, :), t_xyz(2, :), t_xyz(3, :), ...
        tstyle, 'Color', colormatrix, 'Linewidth', 3, ...
        'Parent', axH);
    
    if ~isempty(nodetags) 
        
        hold(axH, 'on');
        
        nodetags = [nodetags*3 - 2, nodetags*3 - 1, nodetags*3];
        nodetags = reshape(nodetags', [1 numel(nodetags)]);
        
        plot3(t_xyz(1, nodetags), t_xyz(2, nodetags), t_xyz(3, nodetags), ...
            '.', 'Color', [0 0 0], 'Linewidth', 3, ...
            'Parent', axH)
    
    end
    
else
    
    linH = plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3), ...
        tstyle, 'Color', colormatrix, 'Linewidth', 3, ...
        'Parent', axH);

    if ~isempty(nodetags)
        
        hold(axH, 'on');
        
        plot3(xyz(nodetags, 1), xyz(nodetags, 2), xyz(nodetags, 3), ...
        '.', 'Color', [0 0 0], 'Linewidth', 3, ...
        'Parent', axH)
    
    end
    
end

hold(axH, 'on');
    
if psoma
    
    % assumes that first point is soma
    somaH = plot3(xyz(1, 1), xyz(1, 2), xyz(1, 3), ....
        'Color', colormatrix, 'Marker', 'o', ...
        'MarkerSize', 10, 'MarkerEdgeColor', [1 1 1], ...
        'MarkerFaceColor',  colormatrix, 'Parent', axH);
    
else
    
    somaH = 0;
    
end

% edit labels
xlabel(axH, 'X');
ylabel(axH, 'Y');
zlabel(axH, 'Z');

end
