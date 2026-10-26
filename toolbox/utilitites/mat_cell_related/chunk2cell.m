function [C, G] = chunk2cell(X, chunkSize, dim)
% chunk2cell: Split an array into chunks of a specific size.
%
% Usage:
%   [C,G] = chunk2cell(X, chunkSize, dim)
%
% Args:
%   X: any array
%   chunkSize: max number of elements per chunk
%   dim: dimension along which to split (default: 1; 2 for horizontal vectors)
% 
% Returns:
%   C: cell array of chunks (numel = ceil(size(X,dim)/chunkSize)); the last
%       chunk is shorter when the length is not a multiple of chunkSize
%   G: a grouping variable of the same length as size(X,dim) indicating the
%       chunk index each number corresponds to
% 
% See also: grp2cell

if nargin < 3; dim = []; end
if isempty(dim) && isvector(X)
    % vectors are split along their long dimension (same as argmax, without depending on a function outside this repo)
    [~, dim] = max(size(X));
end
if isempty(dim); dim = 1; end

N = size(X,dim);
G = ceil(vert(1:N) ./ chunkSize);

if dim == 1 && ismatrix(X)
    
    numChunks = ceil(N / chunkSize);
    D = ones(numChunks,1) .* chunkSize;
    D(end) = D(end) + (N - numChunks * chunkSize);
    C = mat2cell(X, D);
    
else
    
    C = grp2cell(X, G, dim);
    
end

end
