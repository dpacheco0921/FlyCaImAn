function data = imblur(data, sig, siz, nDimBlur)
% imblur: Gaussian blur for high dimensional data
%
% Usage:
%   data = imblur(data, sig, siz, nDimBlur)
%
% Args:
%   data: original data
%   sig: std of gaussian kernel
%   siz: size of kernel
%   nDimBlur: number of dims to blur(optional)
%
% Returns:
%   data: result after the Gaussian blur
%
% Notes
% Author: Yuanjun Gao

if ~exist('nDimBlur', 'var')
    nDimBlur = ndims(data) - 1;
end

if length(sig) == 1 % isotropic case
    sig = sig * ones(1, nDimBlur);
end
assert(nDimBlur == length(sig));

if length(siz) == 1
    siz = siz * ones(1, nDimBlur);
end
assert(nDimBlur == length(siz));

for i = 1:nDimBlur
    
    if sig(i)
        x = -floor(siz(i) / 2):floor(siz(i) / 2);
        H = exp(-(x.^2/ (2 * sig(i)^2)));
        H = H' / sum(H(:));
        
        if nDimBlur > 1
            indH = 1:nDimBlur;
            indH(i) = 1;
            indH(1) = i;
            H = permute(H, indH);
        end
        
        data = imfilter(data, H, 'same', 0);
        
    end
    
end
