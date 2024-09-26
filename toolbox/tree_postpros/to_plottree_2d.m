function to_plottree_2d(...
    iObject, idpar, axis2proj, axH, ...
    verbose, colormatrix, tstyle, psoma)
% to_plottree_2d: plot traces in 2D format
%
% Usage:
%   to_plottree_2d(...
%       iObject, idpar, axis2proj, axH, ...
%       verbose, colormatrix, tstyle, psoma)
%
% Args:
%   iObject: tree or xyz matrix
%   idpar: indeces of parents to xyz
%   axis2proj: which axis to use
%   axH: axes handle
%   verbose: verbose
%   colormatrix: matrix with colors to use
%   tstyle: line style
%   psoma: plot soma

% deafult params

if ~exist('idpar', 'var') || ...
        isempty(idpar)
    idpar = [];
end

if ~exist('axis2proj', 'var') || ...
        isempty(axis2proj)
    axis2proj = 'XY';
end

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


xyz = getObj_xyz(iObject);

tree_siz = size(xyz);

if verbose
    fprintf(['Number of NaNs ', ...
        num2str(sum(isnan(xyz(:)))), '\n']);
end

% define axis
switch axis2proj
    
    case 'XY'
        
        c_idx = [1 2];
    case 'YZ'
        
        c_idx = [2 3];
    case 'XZ'
        
        c_idx = [1 3];
        
end

if ~isempty(idpar)
    
    k = 1;

    for j = 1:tree_siz(1)

        t_xyz(:, k:k+1) = [xyz(idpar(j), :)', xyz(j, :)'];
        t_xyz(:, end + 1) = nan;
        k = size(t_xyz, 2) + 1;

    end

    plot(t_xyz(:, c_idx(1)), t_xyz(:, c_idx(2)), ...
        tstyle, 'Color', colormatrix, 'Linewidth', 3, ...
        'Parent', axH);
        
else
    
    plot(xyz(:, c_idx(1)), xyz(:, c_idx(2)), ...
    	tstyle, 'Color', colormatrix, 'Linewidth', 3, ...
        'Parent', axH);
    
    hold(axH, 'on')
    
end

if psoma % assumes that first point is soma
    plot(xyz(1, 1), xyz(1, 2), ....
        'Color', colormatrix, 'Marker', 'o', ...
        'MarkerSize', 10, 'MarkerEdgeColor', [1 1 1], ...
        'MarkerFaceColor',  colormatrix, 'Parent', axH)
end

% edit labels
switch axis2proj
    
    case 'XY'
        
        xlabel(axH, 'X');
        ylabel(axH, 'Y');
        
    case 'YZ'
        
        xlabel(axH, 'Y');
        ylabel(axH, 'Z');
        
    case 'XZ'
        
        xlabel(axH, 'X'); 
        ylabel(axH, 'Z');
        
end

end
