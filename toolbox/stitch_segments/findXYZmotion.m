function coridx_per_xyz = findXYZmotion(refIm, floatIm, maxshift)
% findXYZmotion: gets the correlation of refIm to floatIm plus some shift
%   in x and y (from -maxshift to +maxshift)
%
% Usage:
%   coridx_per_xyz = findXYZmotion(refIm, floatIm, maxshift)
%
% Args:
%   refIm: reference image [height width]
%   floatIm: floating image(s) [height width planes], same height and width
%   maxshift: maximun pixel shift in X/Y (pixels)
%
% Returns:
%   coridx_per_xyz: correlation coefficients (unitless) of refIm with each
%       plane of floatIm, [2*maxshift+1 (rows shift), 2*maxshift+1 (columns
%       shift), planes]. NaN pixels give NaN correlations.

xrange = -maxshift:1:maxshift;
yrange = -maxshift:1:maxshift;
coridx_per_xyz = zeros(numel(xrange), numel(yrange), size(floatIm, 3));

for z = 1:size(floatIm, 3)
    
    for x = xrange
        
        xi = find(xrange == x);
        
        for y = yrange
            
            yi = find(yrange == y);
            coridx_per_xyz(xi, yi, z) = ...
                corr2(refIm((1+maxshift):(end-maxshift-1), ...
                (1+maxshift):(end-maxshift-1)), ...
                floatIm((1+maxshift+x):(end-maxshift-1+x), ...
                (1+maxshift+y):(end-maxshift-1+y), z));
            
        end
        
    end
    
end

end
