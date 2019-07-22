function iIm = pruneIm(iIm, mask)
% pruneIm: uses a 2D mask to prune volumes (deletes edges: rows and columns)
%
% Usage:
%   newvar_sorted = st_add_newvarperroi(obj, fly_idx, roi_idx, input_var)
%
% Args:
%   iIm: input volume
%   mask: 2D/3D image with 0 and 1, zero are pixels to keep
%
% Notes
% 3D images are max projected

if ~isempty(iIm)
    
    dDim = size(mask);
    
    if size(mask, 3)
        mask = max(mask, [], 3) > 0;
    end
    
    X = sum(mask, 2) == dDim(2);
    Y = sum(mask, 1) == dDim(1);
    
    if length(size(iIm)) == 3
        iIm = iIm(X~= 1, Y~= 1, :);
    elseif length(size(iIm)) == 4
        iIm = iIm(X~= 1, Y~= 1, :, :);
    end
    
else
    
    iIm = [];
    
end

end