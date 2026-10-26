function S = load_sel_vars(filePath, var2load)
% load_sel_vars: load selected variables 'var2load'
%   from file 'filePath', and saves them in structure S
%
% Usage:
%   S = load_sel_vars(filePath, var2load)
% 
% Args:
%   filePath: file path
%   var2load: string or cell of strings; variables in the file whose name equals
%       one of them (ignoring upper/lower case) are loaded
%
% Returns:
%   S: structure with one field per loaded variable (empty struct if none matched)
%
% See also: str2match, loadvar

info = whos('-file', filePath);
allvars = {info.name};
allvars = str2match(var2load, allvars, 1); % exact name match, ignoring case (str2match mtype = 1)

% start from an empty struct so an empty result is returned if nothing matches
S = struct;

for i = 1:numel(allvars)
    S.(allvars{i}) = loadvar(filePath, allvars{i});
end

end
