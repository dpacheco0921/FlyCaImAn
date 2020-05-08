function oVar = load_variable_across_files(Filename, ...
    filesuffix, variable_name, fieldstocompile, dir_depth)
% load_variable_across_files: function that load and compiles specific
%   field from set variables
%
% Usage:
%   oVar = load_variable_across_files(Filename, ...
%    filesuffix, variable_name, fieldstocompile, dir_depth)
%
% Args:
%   Filename: name of files to load
%   filesuffix: suffix of files to load
%   variable_name: name of variable to load
%   fieldstocompile: fields from variable to compile
%   dir_depth: depth of directory for search

if ~exist('dir_depth', 'var') || isempty(dir_depth)
    dir_depth = 0;
end

if ~exist('fieldstocompile', 'var') || isempty(fieldstocompile)
    fieldstocompile = [];
end

% get files to load
if dir_depth == 0
    f2run = rdir(['.', filesep, '*', ...
        filesuffix]);
elseif dir_depth == 1
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', filesuffix]);
else
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', filesep, '*', ...
        filesuffix]);
end

filename = str2match(Filename, f2run);
filename = {filename.name}';
[filename, iDir_1] = split_path(filename);

for i = 1:numel(filename)    
    
    fprintf('*')
    load([iDir_1{i}, filesep, filename{i}], variable_name)
    
    if ~isempty(fieldstocompile)
        
        for j = 1:numel(fieldstocompile)
            eval(['oVar(', num2str(i), ').', fieldstocompile{j}, ...
                ' = ', variable_name, '.', fieldstocompile{j}, ';'])
        end
        
    else
        
        eval(['oVar(', num2str(i), ') = ', variable_name, ';'])
        
    end
    
    eval(['clear ', variable_name])
    
end
fprintf('Done\n')

end



















