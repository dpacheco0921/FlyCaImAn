function i_str = addsuffix(i_str, i_suffix)
% addsuffix: function to add suffix at the end of a string
%
% Usage:
%   i_str = addsuffix(i_str, i_suffix)
%
% Args:
%   i_str: is a string or a cell composed by strings
%   i_suffix: suffix string

if ~isempty(i_str)
    
    if iscell(i_str)
        
        for i = 1:numel(i_str)
            i_str{i} = [i_str{i}, i_suffix];
        end
        
    else
        i_str = [i_str, i_suffix];
    end
    
end

end