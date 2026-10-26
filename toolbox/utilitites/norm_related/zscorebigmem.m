function Y = zscorebigmem(Y, chunk_siz)
% zscorebigmem: z-score X to the last dimention 
%   1DxT or 2DxT or 3DxT matrices
%
% Usage:
%   Y = zscorebigmem(Y, chunk_siz)
%
% Args:
%   Y: 1DxT or 2DxT or 3DxT matrix
%   chunk_siz: chunk size to run at a time in dimension 1 (default, 5e3)
%
% Returns:
%   Y: double array, z-scored along the last dimension (a vertical vector input
%       is returned as a horizontal vector)
%
% Notes:
%   the std normalization is equivalent to std(A, 1) (divides by N, not N-1)

if ~exist('chunk_siz', 'var') || isempty(chunk_siz)
    chunk_siz = 5e3;
end

% flip dimension if a vertical vector is provided
if length(size(Y)) <= 2 && ...
    size(Y, 2) == 1
    Y = Y';
end

Y = double(Y);
Y = centerbigmem(Y, chunk_siz);
Y = sdnormbigmem(Y, chunk_siz);

end
