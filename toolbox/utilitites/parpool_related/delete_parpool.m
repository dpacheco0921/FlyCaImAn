function delete_parpool(ppobj)
% delete_parpool: shutting down parpool
%
% Usage:
%   delete_parpool
%   delete_parpool(ppobj)
% 
% Args:
%   ppobj: parpool object returned by setup_parpool (optional; default, the
%       current pool). Nothing is returned.
% 
% See also: parcluster parpool

if ~exist('ppobj', 'var') || isempty(ppobj)
    ppobj = [];
end

% with no pool object given, shut down whatever pool is running (if any)
if isempty(ppobj)
    delete(gcp('nocreate'));
else
    delete(ppobj);
end

end
