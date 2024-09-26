function to_fly2nsybIVA_xform_perfile(floatTrees, ...
    refIm, xDir, iDir, oDir, esuffix, redo)
% to_fly2nsybIVA_xform_perfile: function that performs reformatting of 
%   trees from fly to nsybIVAi reference coordinates using CMTK
%
% Usage:
%   to_fly2nsybIVA_xform_perfile(floatTrees, ...
%       refIm, xDir, iDir, oDir, esuffix, redo)
%
% Args:
%   floatTrees: name of trees (swc files)
%       naming format (year/month/day)_(flynumber)_(type:w|wm)_(channel:01|02))_(neuronname:P1,etc)).SWC)
%   xDir: directory of transformation
%   iDir: *.swc files directory
%   oDir: output directory
%   esuffix: extra suffix
%   redo: redo gate

if ~exist('esuffix', 'var') || isempty(esuffix); esuffix = ''; end
if ~exist('xDir', 'var') || isempty(xDir); xDir = '.'; end
if ~exist('iDir', 'var') || isempty(iDir); iDir = '.'; end
if ~exist('oDir', 'var') || isempty(oDir); oDir = '.'; end
if ~exist('redo', 'var') || isempty(redo); redo = 0; end

for i = 1:numel(floatTrees)
    
    if exist([iDir, filesep, strrep(floatTrees{i}, '.swc', ...
            ['_', refIm, esuffix, '.swc'])], 'file') ~= 2 || redo
        
        fprintf(['Running file: ', floatTrees{i}, '\n'])
        tree = load_tree([iDir, filesep, floatTrees{i}]);
        
        ibpars.xDir = xDir;
        ibpars.esuffix = ['*', esuffix];
        [~, ~, repoDirs] = customdirs_deafult;
        ibpars.transform_dir = repoDirs{3};

        floatTrees_ = strsplit2(floatTrees{i}, '_');
        floatTrees_ = [floatTrees_{1}, '_', floatTrees_{2},  '_', floatTrees_{3}, '_'];
        
        [tree, err_status] = cmtk_xform(tree, floatTrees_, refIm, ibpars);
        % edit cmtk_streamxform
        % edit cmtk_xform
        
        if err_status; keyboard; end
        
        if numel(tree.X) ~= 0
            
            tree.name = [tree.name, '_', refIm];
            swc_tree(tree, [oDir, filesep, strrep(floatTrees{i}, '.swc', ['_', refIm, esuffix, '.swc'])]);
            
        else
            fprintf('Empty tree - not saving\n')
        end
        
        clear tree floatIm ibpars
        
    else
        fprintf('Already corrected')
    end
    
    fprintf('\n')
    
end

end
