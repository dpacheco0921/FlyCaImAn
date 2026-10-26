function Y = dilate_multiple_matrices(iY, siz_iY, se)
% dilate_multiple_matrices: it dilates a flattened 2D or 3D images 
%   independently.
%
% Usage:
%   Y = dilate_multiple_matrices(iY, siz_iY, se)
%
% Args:
%   iY: flatten 2D or 3D image, where columns are different images
%       it assumes each column is binary
%   siz_iY: original 2D or 3D dimensions
%   se: structuring element object
%       or array of structuring element objects (see imdilate)
%
% Returns:
%   Y: sparse binary matrix, same size as iY, with each column dilated
%
% Notes:
%   see imdilate

assert(size(iY, 1) == prod(siz_iY), ...
    'dilate_multiple_matrices: rows of iY do not match prod(siz_iY)')

Y = sparse(size(iY, 1), size(iY, 2));

for im_i = 1:size(iY, 2)
    
    temp_im = reshape(full(iY(:, im_i)), siz_iY);
    temp_im = reshape(imdilate(temp_im, se), [prod(siz_iY), 1]);
    Y(temp_im ~= 0, im_i) = 1;
    
end

end
