function [clus] = hierarchicalClus(...
    Y, clus_n, pmethod, lmethod, ...
    clus_type, pdist_i, linkg_i)
% st_plot_rc_hist: Performs hierarchical clustering on Y
%
% Usage:
%   [clus] = hierarchicalClus(...
%       Y, clus_n, pmethod, lmethod, ...
%       clus_type, pdist_i, linkg_i)
%
% Args:
%   obj: SpaTemp object
%   Y: data [n, T], n = number of rois, T, time.
%   clus_n: number of clusters (when clus_type == 0, 1, 2) 
%       or distance (when clus_type == 3 or 4)
%       (default, 20)
%   pmethod: pairwise metric
%       ('euclidean', 'squaredeuclidean', 'seuclidean', ..
%       'cityblock', 'minkowski', 'chebychev', 'mahalanobis', ...
%       'cosine', 'correlation', 'spearman', 'hamming', 'jaccard')
%       (default, 'correlation')
%   lmethod: linkage method
%       ('average', 'centroid', 'complete', 
%       'median', 'single', 'ward', 'weighted'
%       'centroid' and 'median' methods are non-monotonic)
%       (default, 'ward')
%   clus_type: how to cluster depending on data size
%       for big matrices (matlab has issues loading so many leaves) use
%       clus_type = 2 or 4.
%       (0: find clus_n # clusters (display all leaves))
%       (1: find user define clus_n # clusters, manual input (display all leaves))
%       (2: find clus_n # clusters (for big matrices, dont show all leaves))
%       (3: find clusters using clus_n distance (show all leaves))
%       (4: find clusters using clus_n distance (for big matrices, dont show all leaves))
%       (default, 0)
%   pdist_i: predefined input distance metric
%       (provide only if previously calculated using pdist)
%       (default, [])
%   linkg_i: predefined input linkage values
%       (provide only if previously calculated using linkage)
%       (default, [])
%
% Returns:
%   clus: cluster structure
%       (clus.Y: pairwise distance)
%       (clus.Z: linkage)
%       (clus.lorder: rows sorted by cluster)
%       (clus.label: cluster labels)
%
% Notes
% see pdist linkage dendrogram
% to explore more: I = inconsistent(Z);

if ~exist('pmethod', 'var') || isempty(pmethod)
    pmethod = 'correlation';
end

if ~exist('lmethod', 'var') || isempty(lmethod)
    lmethod = 'complete';
end

if ~exist('clus_n', 'var') || isempty(clus_n)
    clus_n = 20;
end

if ~exist('clus_type', 'var') || isempty(clus_type)
    clus_type = 0;
end

if exist('pdist_i', 'var') && ~isempty(pdist_i)
    clus.Y = pdist_i;
else
    clus.Y = pdist(Y, pmethod);
end

if exist('linkg_i', 'var') && ~isempty(linkg_i)
    clus.Z = linkg_i;
else
    clus.Z = linkage(clus.Y, lmethod);
end

% Type of visualization
if clus_type == 0
    
    % 1) find clus_n clusters (display all leaves)
    
    % plot whole tree and decide on the number of clusters
    figure();
    [~, ~, clus.lorder] = dendrogram(clus.Z, 0);
    close(gcf);
    
    clus.label = cluster(clus.Z, 'MaxClust', clus_n);    
    
    % correct labels and order to reflect dendrogram plot
    clus = hc_correct_label_a(clus);
        
elseif  clus_type == 1
    
    % 2) find user define clusters, manual input (display all leaves)
    
    % plot whole tree and decide on the number of clusters
    figure();
    [~, ~, clus.lorder] = dendrogram(clus.Z, 0);
    clus_n = input('number of clusters to use : ');
    close(gcf)
    
    clus.label = cluster(clus.Z, 'MaxClust', clus_n);

    % correct labels and order to reflect dendrogram plot
    clus = hc_correct_label_a(clus);

elseif clus_type == 2
    
    % 3) find clus_n clusters (for big matrices, dont show all leaves)
  
    figure();
    [~, clus.label, clus.lorder] = dendrogram(clus.Z, clus_n);
    close(gcf);
    
    clus = hc_correct_label_b(clus);

    % order traces within cluster
    clus.lorder = order_within_cluster(Y, ...
        clus.label, pmethod, lmethod);
    
elseif clus_type == 3
    
    % 4) find clusters using distance (show all leaves)
    
    % plot whole tree and decide on the number of clusters
    figure();
    [~, ~, clus.lorder] = dendrogram(clus.Z, 0);
    close(gcf);
    
    clus.label = cluster(clus.Z, 'cutoff', clus_n, 'Criterion', 'distance');    
    
    % correct labels and order to reflect dendrogram plot
    clus = hc_correct_label_a(clus);
    
elseif clus_type == 4
    
    % 5) find clusters using distance (for big matrices, dont show all leaves)
    
    % get clusters
    clus.label = cluster(clus.Z, 'cutoff', clus_n, 'Criterion', 'distance');
    
    figure();
    [~, clus.label, clus.lorder] = ...
        dendrogram(clus.Z, numel(unique(clus.label)));
    close(gcf);
    
    clus = hc_correct_label_b(clus);

    % order traces within cluster
    clus.lorder = order_within_cluster(Y, ...
        clus.label, pmethod, lmethod);
    
end

clear Y

end

function lorder = order_within_cluster(...
    Y, label, pmethod, lmethod)
% order_within_cluster: cluster leaves within 
%   clustered data, providing an ordered data within clusters
%
% Usage:
%   lorder = order_within_cluster(...
%       Y, label, pmethod, lmethod)
%
% Args:
%   Y: data [n, T], n = number of rois, T, time.
%   label: cluster labels per Y row
%   pmethod: pairwise metric
%   lmethod: linkage method

% sort labels
[C, label_idx] = sort(label);

lorder = [];

for i = unique(C)'
    
    fprintf('*')

    idx2use = label_idx(C == i);

    Yi = pdist(Y(idx2use, :), pmethod);
    Zi = linkage(Yi, lmethod);
    [~, ~, lorder{i, 1}] = dendrogram(Zi, 0);
    lorder{i, 1} = idx2use(lorder{i});

    clear Yi Zi idx2use
    
end

lorder = cell2mat(lorder);

end
