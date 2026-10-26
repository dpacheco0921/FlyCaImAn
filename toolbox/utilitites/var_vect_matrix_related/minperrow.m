function idx = minperrow(imatrix)
% minperrow: get the column index of the minimum value in each row
%   (if several columns tie, the first one is returned)
% 
% Usage:
%   idx = minperrow(imatrix)
%
% Args:
%   imatrix: input matrix, N rows x M columns
%
% Returns:
%   idx: N x 1 column vector, column index of the min of each row
%       (NaNs are ignored; an all-NaN row returns 1)

% min along dimension 2 returns one index per row even when values tie
[~, idx] = min(imatrix, [], 2);

end
