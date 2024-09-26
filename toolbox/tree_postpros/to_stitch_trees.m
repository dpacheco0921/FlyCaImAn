function tree = to_stitch_trees(tree)
% to_stitch_trees: stitch trees, it uses the shortest distance between
%   trees to stitch them toguether
%
% Usage:
%   tree = to_stitch_trees(tree)
%
% Args:
%   tree: tree
%
% Notes: see cat_tree

if isstruct(tree)
    tree_tp = [];
    for ii = 1:numel(tree)
        tree_tp{ii} = tree(ii);
    end
    tree = tree_tp;
    
end

tree_siz = cell2mat(cellfun(@(x) numel(x.X), ...
    tree, 'UniformOutput', false));
[tree_siz, idx_tree] = sort(tree_siz, 'descend');
tree = tree(idx_tree);

for ii = 1:numel(tree)
    idpar{ii} = idpar_tree(tree{ii});
end

% get distances of root to other tree
pairs2run = nchoosek(numel(tree):-1:1, 2);
min_d = [];
idx_d = [];

for ii = 1:size(pairs2run, 1)
    
    % run distances
    xyz_1 = getObj_xyz(tree{pairs2run(ii, 1)});
    xyz_2 = getObj_xyz(tree{pairs2run(ii, 2)});
    d_1to2 = pdist2(xyz_1, xyz_2(1, :), 'euclidean');
    d_2to1 = pdist2(xyz_2, xyz_1(1, :), 'euclidean');
    
    % check which direction is the shortest
    [min_d(1, ii), idx_d(1, ii)] = min(d_1to2);
    [min_d(2, ii), idx_d(2, ii)] = min(d_2to1);
    
end

tree_n = numel(tree);

figH = figure();

for ii = 1:tree_n - 1
    axH(ii) = subplot(tree_n - 1, 1, ii);
end

ii = 1;
while tree_n > 1
    
    % find the trees with the minimun distance
    %   between each other and concatenate them
    
    [direction_, pair_] = find(min_d == min(min_d(:)));
        
    if direction_ == 1
        
        fprintf(['concatenating subtree-', ...
            num2str(pairs2run(pair_, 2)), ...
            'to subtree-', num2str(pairs2run(pair_, 1)),'\n'])
        
        to_plottree_3d(tree{pairs2run(pair_, 1)}, ...
            [], axH(ii), 1, 'k', [], 1, []);
        
        tree{pairs2run(pair_, 1)} = ...
            cat_tree(tree{pairs2run(pair_, 1)}, ...
            tree{pairs2run(pair_, 2)}, ...
            idx_d(direction_, pair_), 1);
        
        % plot 
        to_plottree_3d(tree{pairs2run(pair_, 1)}, ...
            [], axH(ii), 1, 'b', [], 1, []);
        to_plottree_3d(tree{pairs2run(pair_, 2)}, ...
            [], axH(ii), 1, 'r', [], 1, []);
    
    else
        
        fprintf(['concatenating subtree-', ...
            num2str(pairs2run(pair_, 1)), ...
            'to subtree-', num2str(pairs2run(pair_, 2)),'\n'])
        
        to_plottree_3d(tree{pairs2run(pair_, 2)}, ...
            [], axH(ii), 1, 'k', [], 1, []);
        
        tree{pairs2run(pair_, 2)} = ...
            cat_tree(tree{pairs2run(pair_, 2)}, ...
            tree{pairs2run(pair_, 1)}, ...
            idx_d(direction_, pair_), 1);
        
        % only perform this repairs
        % eliminate trifurcations by adding short segments:
        tree{pairs2run(pair_, 2)} = ...
            elimt_tree (tree{pairs2run(pair_, 2)});        
        % eliminate 0-length compartments:
        tree{pairs2run(pair_, 2)} = ...
            elim0_tree (tree{pairs2run(pair_, 2)}); 
        
        % plot 
        to_plottree_3d(tree{pairs2run(pair_, 2)}, ...
            [], axH(ii), 1, 'b', [], 1, []);
        to_plottree_3d(tree{pairs2run(pair_, 1)}, ...
            [], axH(ii), 1, 'r', [], 1, []);
        
    end
    
    % overwrite min-dist
    min_d(direction_, pair_) = inf;
    tree_n = tree_n - 1;
    ii = ii + 1;
    
end

tree = tree{tree_siz == max(tree_siz)};

end
