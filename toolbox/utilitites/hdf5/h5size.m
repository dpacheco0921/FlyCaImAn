function [sz, maxSize] = h5size(filepath, dataset, dim)
%H5SIZE Returns the size of the specified dataset.
% Usage:
%   [size, maxSize] = h5size(filepath, dataset)
%   [size, maxSize] = h5size(filepath, dataset, dim)
%
% Args:
%   filepath: path to the HDF5 file
%   dataset: path to the dataset within the file
%   dim: (optional) dimension(s) to return instead of all
%
% Returns:
%   sz: size of the dataset (number of elements per dimension)
%   maxSize: maximum size of the dataset (only if two outputs requested)
%
% Notes:
%   With a single output it calls h5file (external HDF5 toolbox, not in
%   this repo); with two outputs it only needs h5info.

if nargin < 3; dim = []; end

if nargout > 1
    info = h5info(filepath, dataset);

    sz = info.Dataspace.Size;
    maxSize = info.Dataspace.MaxSize;

    if ~isempty(dim)
        sz = sz(dim);
        maxSize = maxSize(dim);
    end
else
    sz = size(h5file(filepath, dataset));
    if ~isempty(dim)
        sz = sz(dim);
    end
end
end

