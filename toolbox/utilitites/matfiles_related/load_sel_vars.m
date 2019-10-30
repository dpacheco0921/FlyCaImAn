function S = load_sel_vars(filePath, var2load)
% load_sel_vars: load selected variables 'var2load'
%   from file 'filePath', and saves them in structure S
%
% Usage:
%   S = load_sel_vars(filePath, var2load)
% 
% Args:
%   filePath: file path
%   var2load: variable to load
% 
% See also: 

info = whos('-file', filePath);
allvars = {info.name};
allvars = str2match(var2load, allvars, 1);

for i = 1:numel(allvars)
    S.(allvars{i}) = loadvar(filePath, allvars{i});
end

end
