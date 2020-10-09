function figclosekey(figH, key)
% figclosekey: Add a hotkey for closing the figure.
%
% Usage:
%   figclosekey(figH, key)
% 
% Args:
%   figH: figure handle (default: gcf)
%   key: hotkey to use (default: 'q')
% 
% See also: event, addlistener

if nargin < 2 || isempty(key)
    key = 'q';
end

if nargin == 1 && ischar(figH)
    [figH, key] = swap(figH, key);
end

if nargin < 1 || isempty(figH)
    figH = gcf();
end

set(figH, 'KeyPressFcn', ...
    @(figH, evt)KeyPressFcn_cb(figH, evt, key));

end

function KeyPressFcn_cb(figH, evt, key)

if strcmp(evt.Key, key)
    delete(figH)
end

end
