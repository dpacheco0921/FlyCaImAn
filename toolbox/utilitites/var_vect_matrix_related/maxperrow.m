function idx = maxperrow(imatrix)
% maxperrow: get the column index of the maximum value in each row
%   (if several columns tie, the first one is returned)
% 
% Usage:
%   idx = maxperrow(imatrix)
%
% Args:
%   imatrix: input matrix, N rows x M columns
%
% Returns:
%   idx: N x 1 column vector, column index of the max of each row
%       (NaNs are ignored; an all-NaN row returns 1)

% max along dimension 2 returns one index per row even when values tie
%   (the previous find(row == max(row)) returned several indices on ties
%   and nothing on all-NaN rows, which made the assignment fail)
[~, idx] = max(imatrix, [], 2);

end
