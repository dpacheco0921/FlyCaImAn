function [animalname, animalnum, repnum, str_length] = ...
    rdir_namesplit(basename, suffixend, ...
    fisuffix, fi2reject, fi2match, nametype)
% rdir_namesplit: function to gather all files 
%   within directory and split the name into parts
%
% Usage:
%   [animalname, animalnum, repnum] = ...
%      rdir_namesplit(basename, suffixend, ...
%      fisuffix, fi2reject, fi2match, nametype)
%
% Args:
%   basename: file name pattern
%   suffixend: file format
%   fisuffix: extra pattern to attach to suffixend
%   fi2reject: pattern to reject files
%   fi2match: pattern to accept files
%   nametype: type of animalname to output
%
% Returns:
%   animalname: animal name
%   animalnum: animal number
%   repnum: repetition number
%   str_length: string length for first file
%
% Notes:
% Assumes files names have the following struture:
%   (year|month|day)_(animalnum)_(trialnum/repnum)

if ~exist('suffixend', 'var') || isempty(suffixend)
    suffixend = '.tif';
end

if ~exist('nametype', 'var') || isempty(nametype)
    nametype = 0;
end

if exist('fisuffix', 'var') && ~isempty(fisuffix)
    suffixend = [fisuffix, suffixend];
end

str_length = [];

f2run = rdir(['.', filesep, '*', suffixend]);
f2run = str2rm(fi2reject, f2run);

if exist('fi2match', 'var')  && ~isempty(fi2match)
    f2run = str2match(fi2match, f2run); 
end

f2run = str2match(basename, f2run);
f2run = {f2run.name};
% remove dir preffix
f2run = strrep(f2run, ['.', filesep], '');

% initialize output params
animalname = cell(1, numel(f2run));
animalnum = zeros(1, numel(f2run));
repnum = zeros(1, numel(f2run));

for FNum = 1:numel(f2run)
    
    TempS = strsplit2(f2run{FNum}, '_');
    
    if nametype == 0 
        
        % split string into: animalname , animalnum , reps
        animalname{1, FNum} = TempS{1};
        animalnum(1, FNum) = str2double(TempS{2});
        repnum(1, FNum) = str2double(TempS{3});
        
        str_length(FNum, :) = ...
            [length(TempS{1}), ...
            length(TempS{2}), ...
            length(TempS{3}), ...
            length(strrep(TempS{4}, suffixend, ''))];
        
    elseif nametype == 1
        
        % split string into: animalname + animalnum , reps
        animalname{1, FNum} = [TempS{1}, '_', TempS{2}];
        repnum(1, FNum) = str2double(TempS{3});
        
        str_length(FNum, :) = ...
            [length(animalname{1, FNum}), ...
            length(TempS{3}), ...
            length(strrep(TempS{4}, suffixend, ''))];
        
    elseif nametype == 2 
        
        % does not split: animalname + animalnum + reps
        animalname{1, FNum} = [TempS{1}, '_', TempS{2}, '_', TempS{3}];
        
        str_length(FNum, :) = ...
            [length(animalname{1, FNum}), ...
            length(strrep(TempS{4}, suffixend, ''))];
        
    end
    
end

end
