function Y = sqrtbigmem(Y, chunk_siz)
% sqrtbigmem: calculate square root per raw for big 2D matrices
%
% Usage:
%   Y = sqrtbigmem(Y, chunk_siz)
%
% Args:
%   Y: 2DxT matrix
%   chunk_siz: chunk size to run at a time in dimension 1

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
        Y(lidx, :) = sqrt(Y(lidx, :));
        clear lidx
        
    end
    
else
    Y = sqrt(Y);
end

end
