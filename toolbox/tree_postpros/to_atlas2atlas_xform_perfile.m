function to_atlas2atlas_xform_perfile(floatTrees, ...
    refi, refo, xDir, iDir, oDir, esuffix, redo)
% to_atlas2atlas_xform_perfile: function that performs 
%   reformatting of trees from atlas to atlas coordinates using CMTK
%
% Usage:
%   to_atlas2atlas_xform_perfile(floatTrees, ...
%       refi, refo, xDir, iDir, oDir, esuffix, redo)
%
% Args:
%   floatTrees: cell array with names of trees (swc files)
%       naming format (year/month/day)_(flynumber)_(type:w|wm)_(channel:01|02))_(neuronname:P1,etc)).SWC)
%       ending in '_<refi>.swc'
%   refi: input reference atlas
%   refo: output reference atlas
%   xDir: directory of transformation (default, '.')
%   iDir: *.swc files directory (default, '.')
%   oDir: output directory (default, '.')
%   esuffix: extra suffix (default, '')
%   redo: redo gate, 1 to redo trees already reformatted (default, 0)
%
% Returns:
%   nothing; saves '<name>_<refo><esuffix>.swc' in oDir for each tree.
%   Trees that fail are reported (with the error message) and skipped.

if ~exist('esuffix', 'var') || isempty(esuffix)
    esuffix = '';
end

if ~exist('xDir', 'var') || isempty(xDir)
    xDir = '.';
end

if ~exist('iDir', 'var') || isempty(iDir)
    iDir = '.';
end

if ~exist('oDir', 'var') || isempty(oDir)
    oDir = '.';
end

if ~exist('redo', 'var') || isempty(redo)
    redo = 0;
end

for i = 1:numel(floatTrees)
    
    % skip trees whose reformatted swc already exists in the OUTPUT folder
    %   (it was checked in iDir, but the result is saved in oDir)
    if exist([oDir, filesep, strrep(floatTrees{i}, ['_', refi, '.swc'], ...
            ['_', refo, esuffix, '.swc'])], 'file') ~= 2 || redo
        
        try 
            fprintf(['Running file: ', strrep(floatTrees{i}, '\', ' '), '\n'])

            % flag to not repair
            option = 'no';
            tree = load_tree([iDir, filesep, floatTrees{i}], option);

            if numel(tree) > 1

                % swc file with disconnected pieces: join them
                tree = to_stitch_trees(tree);

            end
            
            % only perform this repairs
            % eliminate trifurcations by adding short segments:
            tree = elimt_tree (tree);        
            % eliminate 0-length compartments:
            tree = elim0_tree (tree); 
            
            % perform reformatting
            [tree, err_status] = ...
                to_atlas2atlas_xform(tree, refi, refo, xDir, esuffix);

            if err_status
                % caught below, which reports the file name and goes on
                error('reformatting returned an error status')
            end

            if numel(tree.X) ~= 0

                tree.name = strrep(tree.name, ['_', refi], ['_', refo]);
                [filename, ~] = split_path(floatTrees{i});
                swc_tree(tree, [oDir, filesep, ...
                    strrep(filename, ['_', refi, '.swc'], ...
                    ['_', refo, esuffix, '.swc'])]);

            else

                fprintf('Empty tree - not saving\n')

            end

            clear tree floatIm ibpars
            
            
        catch ME

            fprintf(['tree: ', floatTrees{i}, ' failed: ', ME.message, '\n'])
            
        end
        
    else
        
        fprintf('Already corrected')
        
    end
    
    fprintf('\n')
    
end

end
