function tokens = strsplit2(string, delimiter)
% strsplit2: splits input string into cells using delimiter
%
% Usage:
%   tokens = strsplit2(string, delimiter)
%
% Args:
%   string: input a string
%   delimiter:can be a single character or a string of characters, each
%       single one of which is used to split STRING
% 
% Returns:
%   tokens: splited strings

tokens = regexp(string, ['\s*\' delimiter '\s*'], 'split');

end
