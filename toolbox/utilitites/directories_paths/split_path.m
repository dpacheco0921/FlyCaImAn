function [iFile, iDir] = split_path(iFile)
% split_path: it splits path into file name and full directory
% 
% Usage:
%   iFile = split_path(InputFiles)
%
% Args:
%   iFile: input files to chop
% 
% Returns:
%   iFile: file name(s) in cells
%   iDir: file directory(s) in cells

idelimiter = filesep;

if iscell(iFile)
    
    for i = 1:numel(iFile)
        
        [iDir{i, 1}, ~, ~] = fileparts(iFile{i});
        iFile{i} = strrep(iFile{i}, [iDir{i, 1}, idelimiter], '');
        
    end
    
else
    
    [iDir, ~, ~] = fileparts(iFile);
    iFile = strrep(iFile, [iDir, idelimiter], '');
    
end

end
