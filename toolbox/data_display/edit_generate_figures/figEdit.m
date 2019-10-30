function figEdit(axhandle, fighandle, axcolor, ...
    figcolor, xyzcolor, tickgate, fontsiz, fontsizrel)
% figEdit: function to edit basic features of figures for default vizualization
%   background
%
% Usage:
%   figEdit(axhandle, fighandle, axcolor, ...
%       figcolor, xyzcolor, tickgate, fontsiz, fontsizrel)
%
% Args:
%   axhandle: axes(s) handle
%   fighandle: figure(s) handle
%   axcolor: color of axes background
%       (default, 'w')
%   figcolor: color of figure background 
%       (default, 'w')
%   xyzcolor: color of axis
%       (default, 'k')
%   tickgate: gate to show ticks (on, off)
%       (default, 'on')
%   fontsiz: font size per axes
%       (default, 15)
%   fontsizrel: fontsize multiplier
%       (default, 1.5)
%
% Notes:
%
% ToDo: edit also size features
%   figHandle.Units = 'centimeters';
%   figHandle.PaperPositionMode = 'manual';
%   figHandle.PaperSize = figHandle.PaperSize/10;
%   figHandle.PaperPosition = figHandle.PaperPosition/10;

if ~exist('axcolor', 'var') || isempty(axcolor); axcolor = 'w'; end
if ~exist('figcolor', 'var') || isempty(figcolor); figcolor = 'w'; end
if ~exist('xyzcolor', 'var') || isempty(xyzcolor); xyzcolor = 'k'; end
if ~exist('tickgate', 'var') || isempty(tickgate); tickgate = 'on'; end
if ~exist('fontsiz', 'var') || isempty(fontsiz); fontsiz = 15; end
if ~exist('fontsizrel', 'var') || isempty(fontsizrel); fontsizrel = 1.5; end

for i = 1:numel(axhandle)
    try
        
        set(axhandle(i), 'Box', 'off');
        
        set(axhandle(i), 'color', axcolor);
        
        if ~iscell(xyzcolor)
            set(axhandle, 'xcolor', xyzcolor); 
            set(axhandle, 'ycolor', xyzcolor);
            set(axhandle, 'zcolor', xyzcolor);
        else
            set(axhandle, 'xcolor', xyzcolor{1}); 
            set(axhandle, 'ycolor', xyzcolor{2});
            set(axhandle, 'zcolor', xyzcolor{3});            
        end
        
        if contains(tickgate, 'off')
            set(axhandle, 'xtick', []); 
            set(axhandle, 'ytick', []);
            set(axhandle, 'ztick', []);
        end
        
        set(axhandle(i), 'FontSize', fontsiz, ...
            'LabelFontSizeMultiplier', fontsizrel);
        set(axhandle(i), 'FontName', 'Arial');
        
    end
end

for i = 1:numel(fighandle)
    set(fighandle(i), 'color', figcolor);
end

end