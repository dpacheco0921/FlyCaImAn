function i_file = str2rm(i_str, i_file)
% str2rm: function to find group of strings with
%   non-overlapping names to i_str
%
% Usage:
%   i_file = str2rm(i_str, i_file)
%
% Args:
%   i_str: is a string or a cell of strings (pattern to remove)
%   i_file: is a string or a cell of strings (input files)
% 
% Returns:
%   i_file: subset of i_file that matches i_str

if ~isempty(i_str)
    
    if ~iscell(i_str)
        i_str = {i_str};
    end
    
    Strnum = numel(i_str);
    
    if isstruct(i_file)
        
        for i = 1:Strnum
            file2clear = strfind({i_file.name}, i_str{i});
            i_file(~cellfun(@isempty, file2clear)) = [];
            clear file2clear
        end
        
    else
        
        for i = 1:Strnum
            file2clear = strfind(i_file, i_str{i});
            i_file(~cellfun(@isempty, file2clear)) = [];
            clear file2clear
        end
        
    end
    
else
    
    %fprintf('no input str \n')
    
end

end
