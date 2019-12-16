function idx = minperrow(imatrix)
% minperrow: get the indeces of the maximun 
%   value across columns per row
% 
% Usage:
%   idx = maxperrow(imatrix)
%
% Args:
%   Y: input matrix

idx = zeros(size(imatrix, 1), 1);

for i = 1:size(imatrix, 1)
    
    idx(i, 1) = find(imatrix(i, :) == min(imatrix(i, :)));
    
end 

end