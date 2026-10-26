function [Y, cm] = labeledim_to_imperlabel(Im)
% labeledim_to_imperlabel: it converts a labeled matrix (each label is an integer)
%   and it expands an aditional dimension where each point refers to each
%   label
%
% Usage:
%   [Y, cm] = labeledim_to_imperlabel(Im)
%
% Args:
%   Im: input labeled image (2D or 3D); 0 and NaN are background
%
% Returns:
%   Y: sparse matrix [numel(Im), n_labels], column k is the binary
%       (flattened) mask of the k-th smallest label
%   cm: center of mass per label [n_labels, 3] (pixels; from com)

% get original size
d1 = size(Im, 1);
d2 = size(Im, 2);
d3 = size(Im, 3);

% flatten image
Im = Im(:);

% get all labels (discard nan and 0)
labels = unique(Im);
labels(isnan(labels)) = [];
labels(labels == 0) = [];

% make sparse matrix
k = 1;
Y = sparse(numel(Im), numel(labels));
cm = zeros(numel(labels), 3);

for i = labels'
    Y(Im == i, k) = 1;
    k = k + 1;
end

cm = com(Y, d1, d2, d3);

end
