function write_table_in_textfile(filename, oDir, variable2write)
% write_table_in_textfile: writes a text file from table or cell of
%   variables
%
% Usage:
%   write_table_in_textfile(filename, oDir, variable2write)
%
% Args:
%   filename: name of file to save
%   oDir: output directory
%   variable2write: variable to write as text (table or cell)
%
% Notes:
%   if input variable is a cell, it assumes that the first row are variable
%       names

if istable(variable2write)
    
    % if table, extract variable names and variables
    inputVar = table2cell(variable2write);
    varNames = variable2write.Properties.VariableNames;
    
elseif iscell(variable2write)
    
    % if cell it assumes first row is variable names
    varNames = variable2write(1, :);
    inputVar = variable2write(2:end, :);
    
end

full_name = fullfile(oDir, [filename '.txt']);

% open/create log file
fid = fopen(full_name, 'w');

% set header for the new log file
varName_str = [];
for i = 1:numel(varNames)
    
    % convert numbers to strings
    if ~ischar(varNames{i})
       varNames{i} = num2str(varNames{i}); 
    end
    
    if i ~= numel(varNames)
        varName_str = ...
            [varName_str, varNames{i}, '\t'];
    else
        varName_str = ...
            [varName_str, varNames{i}, '\n'];
    end
    
end
fprintf(fid, varName_str);

% close log file
fclose(fid);

% open/create log file
fid = fopen(full_name, 'a+');

% write each row
for i = 1:size(inputVar, 1)
    
    row_str = [];
    
    for j = 1:size(inputVar, 2)
        
        if ischar(inputVar{i, j})
           str2write = inputVar{i, j};
        else
           str2write = num2str(inputVar{i, j});
        end
        
        if j ~= size(inputVar, 2)
            row_str = ...
                [row_str, str2write, '\t'];
        else
            row_str = ...
                [row_str, str2write, '\n'];
        end
        
    end
    
    fprintf(fid, row_str);
    
end

% close log file
fclose(fid);

end
