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
        
        prestr = strsplit2(iFile{i}, idelimiter);
        iDir{i, 1} = strrep(iFile{i}, ...
            [idelimiter, prestr{end}], '');
        iFile{i, 1} = prestr{end};
        
    end
    
else
    
    prestr = strsplit2(iFile, idelimiter);
    iDir = strrep(iFile, [idelimiter, ...
        prestr{end}], '');
    iFile = prestr{end};
    
end

end
