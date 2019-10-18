function [ppobj, ppprofile] = setup_parpool(deviceID, corenum)
% setup_parpool Returns the index at which the max is found.
% 
% Usage:
%   pphandle = setup_parpool(deviceID, corenum)
% Args:
%   deviceID: 'int', internal, 'spock' or 
%       'della' or other names it assumes they are clusters
%   corenum: number of workers
% 
% Returns:
%   ppobj: parpool obj
%   ppprofile: parpool profile
% 
% See also: parcluster parpool

if ~exist('deviceID', 'var') || isempty(deviceID); deviceID = []; end
if ~exist('corenum', 'var') || isempty(corenum); corenum = 1; end

ppobj = []; ppprofile = [];

if ispc || ismac || contains(deviceID, 'int')
    
    if isempty(gcp('nocreate'))
        ppprofile = parcluster('local');
        ppobj = parpool(ppprofile, corenum);
    end
    
else
    
    ppprofile = parcluster('local');
    ppprofile.JobStorageLocation = strcat('/tmp/', getenv('USER'), ...
        '-', getenv('SLURM_JOB_ID'));
    ppobj = parpool(ppprofile, corenum);
    
end

end
