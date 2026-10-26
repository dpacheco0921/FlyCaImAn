function Y = sqrtbigmem(Y, chunk_siz)
% sqrtbigmem: calculate square root per row for big 2D matrices
%
% Usage:
%   Y = sqrtbigmem(Y, chunk_siz)
%
% Args:
%   Y: 2DxT matrix
%   chunk_siz: chunk size to run at a time in dimension 1 (default, 5e3)
%
% Returns:
%   Y: element-wise square root of Y (same size)

if ~exist('chunk_siz', 'var') || isempty(chunk_siz)
    chunk_siz = 5e3;
end

siz = size(Y);
nDim = length(siz);

if nDim > 2
    % rows per chunk shrinks as the middle dimensions grow; keep at least 1 row so the
    % loop below always runs (a 0 would silently leave Y unprocessed)
    chunk_siz = max(1, round(chunk_siz/prod(siz(2:nDim-1))));
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
