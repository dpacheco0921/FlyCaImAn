function [exist_index] = find_match_file(filename, suffix)
% find_match_file: find files that match the file name + suffix
%
% Usage:
%   find_match_file(filename, suffix)
%
% Args:
%   filename: cell array of filenames
%   suffix: string to add at the end of filename

exist_index = false(numel(filename), 1);

for i = 1:numel(filename)
    exist_index(i) = exist([filename{i}, suffix], 'file');
end

end