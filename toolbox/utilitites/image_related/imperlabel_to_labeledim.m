function Y = imperlabel_to_labeledim(iY, siz_iY)
% imperlabel_to_labeledim: it does the reverse of labeledim_to_imperlabel
%
% Usage:
%   Y = imperlabel_to_labeledim(iY, siz_iY)
%
% Args:
%   iY: input 2D matrix [n_pixels, n_labels], nonzero = pixel in label
%   siz_iY: original image dimensions (2D or 3D)
%
% Returns:
%   Y: labeled image of size siz_iY; pixel value = column index of iY
%       (if labels overlap, the highest column index wins)

% each row of iY must be one pixel of the original image
assert(size(iY, 1) == prod(siz_iY), ...
    'imperlabel_to_labeledim: rows of iY do not match prod(siz_iY)')

Y = zeros(siz_iY);

for i = 1:size(iY, 2)
    Y(iY(:, i) ~= 0) = i;
end

end
