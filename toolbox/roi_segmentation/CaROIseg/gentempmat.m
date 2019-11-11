function dataObj = gentempmat(filename, isuffix)
% gentempmat: memmap file which contains
%   contains fields Y, Yr, sizY, nY.
%
% Usage:
%   dataObj = gentempmat(filename, isuffix)
%
% Args:
%   filename: filename
%   isuffix: suffix that defines the file 
%       [filename, suffix] where all the fields exist.
%
% Notes:
%   Need to add a function to reformat files in case the formating is not
%       ready

% default params
if ~exist('isuffix', 'var') || ...
        isempty(isuffix)
    isuffix = '_prosdata.mat';
end

% remove suffix (.mat) from filename
filename = strrep(filename, '.mat', ''); % remove .mat

% add suffix (.mat) to isuffix and osuffix
if ~contains(isuffix, '.mat')
    isuffix = [isuffix, '.mat'];
end

% check if file with Data or Y variable exist
if ~exist([filename, isuffix], 'file')
    % generate the memmap mat file
    fprintf('Error data file needs to be formatted')
else
    % memmap *.mat file
    dataObj = matfile([filename, isuffix], 'Writable', true);
end

end
