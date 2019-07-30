function rho = neighcorr_2Dof3D(Y)
% neighcorr_2Dof3D: get 2D local correlation of 3D images.
%
% Usage:
%   rho = neighcorr_2Dof3D(Y)
%
% Args:
%   Y: M x N x O x T movie or memmap object (looks for Y.Y)
%
% Returns:
%   rho: M x N x O matrix, cross-correlation with adjacent pixel

memmaped = isobject(Y);

if memmaped
    
    dDim = Y.sizY;
    fprintf(['Running localcorr2Dof3D for big volumes (', ...
        num2str(dDim(3)), ') planes ... '])
    
    for z = 1:dDim(3)
        fprintf('*')
        rho(:, :, z) = neighcorr_2D(squeeze(Y.Y(:, :, z, :)));
    end
    
    fprintf('\n')
    
else
    
    for z = 1:size(Y, 3)
        rho(:, :, z) = neighcorr_2D(squeeze(Y(:, :, z, :)));
    end
    
end

end
