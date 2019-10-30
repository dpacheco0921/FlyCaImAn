function empties = areempty(cellarr)
% empties: returns a logical array of the 
%   size of the cell array indicating whether each cell is empty.
%
% Usage:
%   empties = areempty(cellarr)
%
% See also: isempty

empties = cellfun('isempty', cellarr);

end