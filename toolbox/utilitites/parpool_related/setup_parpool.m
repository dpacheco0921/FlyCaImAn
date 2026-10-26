function [ppobj, ppprofile] = setup_parpool(deviceID, corenum)
% setup_parpool: start a local parallel pool (on a cluster node, use a temp job folder).
% 
% Usage:
%   [ppobj, ppprofile] = setup_parpool(deviceID, corenum)
%
% Args:
%   deviceID: 'int', internal, 'spock' or 
%       'della' or other names it assumes they are clusters
%   corenum: number of workers (default, 1)
% 
% Returns:
%   ppobj: parpool obj
%   ppprofile: parpool profile (both are [] if a pool was already running)
% 
% See also: parcluster parpool

if ~exist('deviceID', 'var') || isempty(deviceID)
    deviceID = ''; % char, so contains(deviceID, ...) below does not error
end

if ~exist('corenum', 'var') || isempty(corenum)
    corenum = 1;
end

ppobj = [];
ppprofile = [];

if ispc || ismac || contains(deviceID, 'int')
    
    if isempty(gcp('nocreate'))
        ppprofile = parcluster('local');
        ppobj = parpool(ppprofile, corenum);
    end
    
else
    
    if isempty(gcp('nocreate'))
        ppprofile = parcluster('local');

        % use node-local temp folder (unique per user and SLURM job) because
        % the default job storage location is not safe on shared cluster disks
        ppprofile.JobStorageLocation = ...
            strcat('/tmp/', getenv('USER'), ...
            '-', getenv('SLURM_JOB_ID'));

        ppobj = parpool(ppprofile, corenum);
    end
    
end

end
