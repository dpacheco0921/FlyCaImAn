function empties = areempty(cellarr)
% areempty: returns a logical array of the
%   size of the cell array indicating whether each cell is empty.
%
% Args:
%   cellarr: cell array
%
% Returns:
%   empties: logical array, same size as cellarr, true where the cell is empty
%
% Usage:
%   empties = areempty(cellarr)
%
% See also: isempty

empties = cellfun('isempty', cellarr);

end