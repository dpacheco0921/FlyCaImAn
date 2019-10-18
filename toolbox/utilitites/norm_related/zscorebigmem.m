function Y = zscorebigmem(Y, chunk_siz)
% zscorebigmem: z-score X to the last dimention 
%   1DxT or 2DxT or 3DxT matrices
%
% Usage:
%   Y = zscorebigmem(Y, chunk_siz)
%
% Args:
%   Y: 1DxT or 2DxT or 3DxT matrix
%   chunk_siz: chunk size to run at a time in dimension 1

if ~exist('chunk_siz', 'var') || isempty(chunk_siz)
    chunk_siz = 5e3;
end

Y = double(Y);
Y = centerbigmem(Y, chunk_siz);
Y = sdnormbigmem(Y, chunk_siz);

end