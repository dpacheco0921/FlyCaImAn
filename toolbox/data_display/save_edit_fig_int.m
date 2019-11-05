function save_edit_fig_int(axhandle, fighandle, ...
    oDir, figname, figformat, fitsize, ...
    axcolor, figcolor, xyzcolor, tickgate, ...
    resolution_, fontsiz, fontsizrel)
% save_edit_fig_int: deafult editing and saving figures
%
% Usage:
%   save_edit_fig_intsave_edit_fig_int(axhandle, fighandle, ...
%   oDir, figname, figformat, fitsize, ...
%   axcolor, figcolor, xyzcolor, tickgate, ...
%   resolution_, fontsiz, fontsizrel)
%
% Args:
%   axhandle: axes(s) handle
%   fighandle: figure(s) handle
%   oDir: target directory
%   figname: figure name
%   figformat: which formats to save (R = [1x12])
%       (1, fig, 2, tif)
%       (3, eps (build-in), 4, eps-ogl (build-in))
%       (5, expfig_eps_nocrop, 6 expfig_eps_crop)
%       (7, expfig_eps-ogl_nocrop, 8 expfig_eps-ogl_crop)
%       (9, png (build-in), 10, expfig_png)
%       (11, svg (build-in), 12 svg)
%   fitsize: fit axes size to whole figure
%       (default, 1)
%   axcolor: color of axes background
%       (default, 'none')
%   figcolor: color of figure background
%       (default, 'none')
%   xyzcolor: color of axis
%       (default, 'k')
%   tickgate: gate to show ticks (on, off)
%       (default, 'on')
%   resolution_: resolution, in dots per pixel
%       (default, -r300)
%   fontsiz: font size per axes
%       (default, 15)
%   fontsizrel: fontsize multiplier
%       (default, 1.5)
%
% Notes: figEdit savefig_int

if ~exist('fitsize', 'var') || isempty(fitsize)
    fitsize = 1;
end

if ~exist('axcolor', 'var') || isempty(axcolor)
    axcolor = 'none';
end

if ~exist('figcolor', 'var') || isempty(figcolor)
    figcolor = 'none';
end

if ~exist('xyzcolor', 'var') || isempty(xyzcolor)
    xyzcolor = 'k';
end

if ~exist('tickgate', 'var') || isempty(tickgate)
    tickgate = 'on';
end

if ~exist('resolution_', 'var') || isempty(resolution_)
    resolution_ = '-r300';
end

if ~exist('fontsiz', 'var') || isempty(fontsiz)
    fontsiz = 15;
end

if ~exist('fontsizrel', 'var') || isempty(fontsizrel)
    fontsizrel = 1.5;
end

if fitsize
    % make axes the size of the figure
    set(axhandle, 'Position', [0 0 1 1])
end

figEdit(axhandle, fighandle, axcolor, ...
    figcolor, xyzcolor, tickgate, ...
    fontsiz, fontsizrel)

if ~isempty(oDir)
    savefig_int(fighandle, oDir, ...
        figname, figformat, resolution_)
else
    fprintf('Error - no oDir provided')
end

end
