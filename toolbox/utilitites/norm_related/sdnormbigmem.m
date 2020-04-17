function Y = sdnormbigmem(Y, chunk_siz)
% sdnormbigmem: sd-norm data to the last dimention
%   1DxT or 2DxT or 3DxT matrices
%
% Usage:
%   Y = sdnormbigmem(Y, chunk_siz)
%
% Args:
%   Y: 1DxT or 2DxT or 3DxT matrix
%   chunk_siz: chunk size to run at a time in dimension 1
%
% Notes
% equation used: 
%   S = ((1/N)sum((Ai??)^2))^0.5 (equivalent to std(A, 1))

if ~exist('chunk_siz', 'var') || isempty(chunk_siz)
    chunk_siz = 5e3;
end

siz = size(Y);
nDim = length(siz);

if nDim > 2
    chunk_siz = round(chunk_siz/prod(siz(2:nDim-1)));
end

if prod(siz(1:end-1)) > chunk_siz
    
    for i = 1:chunk_siz:siz(1)
        
        lidx = (i:min(i + chunk_siz - 1, siz(1)));
        Y(lidx, :, :, :) = bsxfun(@times, Y(lidx, :, :, :), ...
            1./sqrt(nanmean(Y(lidx, :, :, :).^2, length(siz))));
        clear lidx
        
    end
    
else
    
    Y = bsxfun(@times, Y, 1./sqrt(nanmean(Y.^2, length(siz))));
    
end

end
