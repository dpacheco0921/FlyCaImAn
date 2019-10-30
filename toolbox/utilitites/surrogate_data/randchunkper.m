function [Yo, idx_out] = randchunkper(...
    Yi, n_points, chunk_min_size, repn, stim)
% randchunkper: Generate a Random permutation of chunks of a 
%   vector Yi (row vector), different from the original arragement
%
% Usage:
%   [Yo, idx_out] = randchunkper(...
%       Yi, n_points, chunk_min_size, repn, stim)
%
% Args:
%   Yi: row vector to split into chunks
%   n_points: number of split points
%   chunk_min_size: min size of chunk
%   repn: number of permutations
%   stim: vector with stimuli information
%
% Output:
%   Yo: Yi sorted with idx_out
%   idx_out: matrix with all chunk-permuted indeces
% Notes:
% Also, during the permutation of the chunks it excludes consecutive chunks

if ~exist('stim', 'var')
    stim = [];
end

rng('shuffle'); % seems to be an issue with matlab 2018

idx_vect(1, :) = 1:numel(Yi);
Yo = zeros(repn, numel(Yi));
idx_out = zeros(repn, numel(Yi));

for i = 1:repn
    
    i_gate = 0;
    
    while i_gate == 0
        
        % get chunks init and end indeces
        chunk_i_e_idx = get_chunks_i_e(idx_vect, n_points, chunk_min_size);
        
        % split vectos into chunks
        chunks_cell = [];
        chunks_cell = split_vector_chunks(idx_vect, n_points, chunk_i_e_idx);
        
        % permute order, excluding cases of consecutive chunks
        dgate = 1;
        while dgate > 0
            neworder = randperm(n_points + 1); 
            dgate = numel(find(unique(diff(neworder)) == 1));
        end
        
        % permute chunks and stitch
        idx_out(i, :) = cell2mat(chunks_cell(neworder));
        
        % exclude cases where it matches the initial index order
        if ~ismember(idx_out(i, :), idx_vect, 'rows'); i_gate = 1; end
        
        % exclude cases where an input stimuli vector is the same under
        % permutation
        
        if ~isempty(stim)
            if ~ismember(stim(idx_out(i, :)), stim, 'rows')
                i_gate = 1; 
            end 
        end
        
    end
    
    Yo(i, :) = Yi(:, idx_out(i, :));
    
end

% check if there are any repeats
% arrayfun(@(k) max(corr(Yo(1:size(Yo, 1) ~= k, :)', Yo(k,:)')), 1:size(Yo, 1), 'Uniform', 1);

end

function chunk_i_e_idx = get_chunks_i_e(...
    idx_vect, n_points, chunk_min_size)

chunk_i_e_idx = zeros(1, n_points);
chunk_size = [];

for i = 1:n_points
    
    % algorithm # 1: the max size is the whole minus the min size
    % times the chunks
    % chunk_size = [min_size, numel(Idx_temp) - (k - ii + 1)*min_size - sum(chunk_i)];
    % algorithm # 2: the max size is always what is left minus a min size chunk divided by
    % the number of chunks left
    chunk_size = [chunk_min_size, ...
        (numel(idx_vect) - sum(chunk_i_e_idx)...
        - chunk_min_size)/(n_points - i + 1)];
    chunk_i_e_idx(1, i) = ...
        round(rand(1)*diff(chunk_size)) + chunk_size(1);

end

chunk_i_e_idx = [0 cumsum(chunk_i_e_idx) numel(idx_vect)];

end

function chunks_cell = split_vector_chunks(...
    idx_vect, n_points, chunk_i_e_idx)

chunks_cell = cell(1, n_points + 1);

% split idx

for i = 1:n_points + 1
    chunks_cell{i} = idx_vect(...
        chunk_i_e_idx(i) + 1:chunk_i_e_idx(i + 1)); 
end

end
