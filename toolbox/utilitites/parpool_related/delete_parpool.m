function delete_parpool(ppobj)
% delete_parpool: shutting down parpool
%
% Usage:
%   pphandle = setup_parpool(deviceID, corenum)
% 
% Args:
%   ppobj: parpool obj
% 
% See also: parcluster parpool

if ~exist('ppobj', 'var') || isempty(ppobj)
    ppobj = [];
end

if ~exist('ppobj', 'var') || isempty(ppobj)
    delete(gcp('nocreate'));
else
    delete(ppobj);
end

end
