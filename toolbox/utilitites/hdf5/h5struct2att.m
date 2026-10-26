function h5struct2att(filepath, location, S)
%H5STRUCT2ATT Writes attributes to an HDF5 file from a scalar structure.
% Usage:
%   h5struct2att(filepath, location, S)
%
% Args:
%   filepath: path to the HDF5 file
%   location: group or dataset path that receives the attributes
%   S: scalar structure; each field becomes one attribute (logicals are
%       saved as uint8, cellstr as one string joined by newlines)
%   (if called as h5struct2att(filepath, S), the variable name of S is
%   used as location)
%
% See also: h5att2struct

if nargin == 2 && isstruct(location); S = location; location = inputname(2); end
if isempty(location); location = ''; end
if ~startsWith(location,'/'); location = ['/' location]; end

names = fieldnames(S);
for i = 1:numel(names)
    if islogical(S.(names{i})); S.(names{i}) = uint8(S.(names{i})); end
    if iscellstr(S.(names{i})); S.(names{i}) = strjoin(S.(names{i}),'\n'); end
    h5writeatt(filepath, location, names{i}, S.(names{i}))
end


end

