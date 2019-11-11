function patches = construct_patches_seg(sizY, seg_idx)
% construct_patches_seg: constructs the coordinates of the 
% different patches for batch processing using segment indeces
%
% Usage:
%   patches = construct_patches_seg(sizY, seg_idx)
%
% Args:
%   sizY: size of input data (excluding time dimension)        
%   seg_idx: indeces or each segment (vector with integer values)

% dimension of dataset (3d)
dimY = length(sizY);

seg_n = unique(seg_idx, 'stable')';

k = 1;
for i = seg_n
    % start XY
    X1(k) = 1; Y1(k) = 1;
    % end XY
    X2(k) = sizY(1); Y2(k) = sizY(2);
    % Planes to use per segment
    Z1(k) = find(seg_idx == i, 1, 'first');
    Z2(k) = find(seg_idx == i, 1, 'last');
    
    k = k + 1;
end

patches = mat2cell([X1(:), X2(:), ...
    Y1(:), Y2(:), Z1(:), Z2(:)], ...
    ones(numel(X1), 1), 2*dimY);

end