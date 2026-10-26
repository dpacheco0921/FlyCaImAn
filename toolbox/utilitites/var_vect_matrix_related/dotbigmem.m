function X = dotbigmem(A, B, dim, chunksize)
% dotbigmem: Computes dot product as in dot.m 
%   of matrices, where each row or column are 
%   used from matrices A and B. it splits A & B if the matrix is too big
%
% Usage:
%   dotbigmem(A, B, dim, chunksize)
%
% Args:
%   A, B: input 2D matrices of the same size
%   dim: dimension over which to compute the dot product (1 or 2)
%   chunksize: number of columns (dim = 1) or rows (dim = 2) processed at a time
%       (default, 5e3)
%
% Returns:
%   X: dot products; 1 x size(A,2) if dim = 1, size(A,1) x 1 if dim = 2
%
% See also: dot

if ~exist('chunksize', 'var') || ...
        isempty(chunksize)
    chunksize = 5e3;
end

dDim = size(A);

% the dot product collapses dimension "dim", so the output has one value per
% element of the OTHER dimension; chunk over that other dimension
other_dim = 3 - dim;

if dim ~= 1 && dim ~= 2
    error('dotbigmem: dim must be 1 or 2');
end

if dDim(other_dim) > chunksize
    
    for i = 1:chunksize:dDim(other_dim)
        
        lidx = (i:min(i + chunksize - 1, dDim(other_dim)));
        
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
