function figResetYlim(axH, figH, irange)
% figResetYlim: change ylim
%
% Usage:
%   figResetYlim(AxH, figH, irange)
% 
% Args:
%   AxH: axes(s) handle
%   figH: figure(s) handle
%   irange: Xlim range per AxH

for i = 1:numel(axH)
    
    set(axH(i), 'Box', 'off');
    set(axH(i), 'color', 'w');
    
    if size(irange, 1) > 1
        ylim(axH(i), [irange(i, 1) irange(i, end)])
    else
        ylim(axH(i), [irange(1) irange(end)])
    end
    
end

set(figH, 'color', 'w');

end
