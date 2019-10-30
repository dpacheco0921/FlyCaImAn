function rmvar(filename, varargin)
% rmvar: Remove variable from MAT-File
%   RMVAR FILENAME VAR1 VAR2... removes the variables VAR1, VAR2... from
%   the Mat-File FILENAME. If FILENAME has no extension RMVAR looks for
%   FILENAME.mat
%
%   RMVAR('Filename','VAR1','VAR2' ...) does the same.
%
%   Example:
%      % Creates a file 'myfile.mat' containing 2 variables:
%      a='hello';
%      b=magic(3);
%      save myfile.mat a b
%      % Removes the variable 'a' and opens the result
%      clear
%      rmvar myfile a
%      load myfile
%      whos
%
%   F. Moisy
%   Revision: 1.00,  Date: 2008/03/31.
%
%   See also LOAD, SAVE.

% History:
% 2008/03/31: v1.00, first version.
% 2014/06/20: Leo Simon fix bug for non-existing variables

narginchk(2, inf);

% determine which variables exist in this file
WHOS = whos('-file', filename); 
removeThese = {};

for i = 1:numel(varargin)
    
    if sum(contains({WHOS.name}, varargin{i})) == 0
        disp([ varargin{i} ' isn''t saved in ' filename ]); 
    else 
        removeThese = [removeThese, varargin{i}]; 
    end
    
end

% remove overlapping variables
vars = rmfield(load(filename), removeThese); 

% save
save(filename, '-struct', 'vars'); 

end
