function X = dotbigmem(A, B, dim, chunksize)
% dotbigmem: Computes dot product as in dot.m 
%   of matrices, where each row or column are 
%   used from matrices A and B. it splits A & B if the matrix is too big
%
% Usage:
%   dotbigmem(A, B, dim, chunksize)
%
% Args:
%   A, B: input matrices
%   dim: dimension over which to compute the dot product
%   chunksize: size of chunks
%
% See also: dot

if ~exist('chunksize', 'var') || ...
        isempty(chunksize)
    chunksize = 5e3;
end

dDim = size(A);

if dDim(dim) > chunksize
    
    for i = 1:chunksize:dDim(dim)
        
        lidx = (i:min(i + chunksize - 1, dDim(dim)));
        
        if dim == 1
            X(1, lidx) = dot(A(:, lidx), B(:, lidx), dim);
        else
            X(lidx, 1) = dot(A(lidx, :), B(lidx, :), dim);
        end
        
        clear lidx
        
    end
    
else
    
    X = dot(A, B, dim);
    
end

end
