function [i_file, sel_idx] = str2match(i_str, i_file, mtype)
% str2match: function to find group of strings with overlapping names
%
% Usage:
%   [i_file, sel_idx] = str2match(i_str, i_file, mtype)
%
% Args:
%   i_str: is a string or a cell of strings (pattern to look for)
%   i_file: is a string or a cell of strings (input files)
%   mtype: ignore case
%       (1: ignore lower/cap case)
%       (0: otherwise, default)
% 
% Returns:
%   i_file: subset of i_file that matches i_str
%   sel_idx: indeces of selected i_file

if ~exist('mtype', 'var') || isempty(mtype); mtype = 0; end

if ~isempty(i_str)
    
    if ~iscell(i_str); i_str = {i_str}; end
    
    Strnum = numel(i_str);
    
    for Str_idx = 1:Strnum
        
        if iscell(i_file)
            
            if mtype
                BinM = strcmpi(i_file, i_str{Str_idx});
            else
                BinM = contains(i_file, i_str{Str_idx});
            end
            
        else
            
            if mtype
                BinM = strcmpi(i_file, i_str{Str_idx});
            else
                BinM = contains({i_file.name}, i_str{Str_idx});
            end
            
        end
        
        if sum(BinM) == 0
            fprintf('No Match\n');
        else
            f2keep(Str_idx, :) = BinM;
        end
        
    end
    
    if exist('f2keep', 'var')
        
        sel_idx = sum(f2keep, 1) > 0;
        i_file(~sel_idx) = [];
        clear Str2Del
        
    else
        
        i_file(:) = [];
        sel_idx = true(numel(i_file), 1)';
        
    end
    
else
    
    sel_idx = true(numel(i_file), 1)';
    
end

end
