function [chunks, nchunks, chunk_idx] = ...
    ppool_makechunks(chunksiz, corenum, ...
    vect_length, vect_init)
% ppool_makechunks: split a vector into chunks that 
%   could be then used by parpool to run data chunk by chunk
%
% Usage:
%   [chunks, nchunks, chunk_idx] = ...
%       ppool_makechunks(chunksiz, corenum, ...
%       vect_length, vect_init)
%
% Args:
%   chunksiz: size of chunks to split the input vector into
%   corenum: number of cores available
%   vect_length: length of vector to split
%   vect_init: first index of this vector
% 
% Returns
%   chunks: first vector index per chunk
%   nchunks: number of chunks
%   chunk_idx: chunks indeces grouped in corenum batches

if ~exist('vect_init', 'var')
    vect_init = 1;
end

chunks = vect_init:chunksiz:vect_length;
nchunks = numel(chunks);

% indices of the patches for each batch
chunk_idx = arrayfun(@(i) chunks(i:min((i+corenum-1), nchunks)), ...
    1:corenum:nchunks, 'UniformOutput', false);

end
