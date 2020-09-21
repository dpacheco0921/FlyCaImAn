function [axH, subclus, Y_mean_trace, Y_mean_trace_sc, ...
    sub_clus_label, sub_clus_order] = ...
    hc_plot_mtrace_per_clus(Y, clus_label, ...
    link_ths, colorvect, y_gap, textgate, axH, ...
    textorder, idx2use4clus, colorvect_mean, ...
    idx2nan, indvidx, plottype, errortype, add_ref_trace)
% hc_plot_mtrace_per_clus: plot mean trace per cluster,
%   in addition it can get and plot the mean of subclusters 
%   (use link_ths), or plot error shades
%
% Usage:
%   [axH, subclus, Y_mean_trace, Y_mean_trace_sc, ...
%       sub_clus_label, sub_clus_order] = ...
%       hc_plot_mtrace_per_clus(Y, clus_label, ...
%       link_ths, colorvect, y_gap, textgate, axH, ...
%       textorder, idx2use4clus, colorvect_mean, ...
%       idx2nan, indvidx, plottype, errortype, add_ref_trace)
%
% Args:
%   Y: signal [n, T], where n == length(clus_label), T timepoints
%   clus_label: vector with cluster labels 
%   link_ths: linkage threshold for subcluster plotting
%       (default, [])
%   colorvect: matrix with colors per cluster
%       (default, jet(numel(unique(clus_label))))
%   y_gap: gap between traces (y axis)
%       (deafult, y_gap)
%   textgate: gate to show text
%       (default, 1)
%   axH: axes handle
%   textorder: ascend (1), descend (0)
%       (default, 1)
%   idx2use4clus: indeces of Y (in second dimension) to use for clustering
%       (default, [])
%   colorvect_mean: matrix with colors per cluster for mean value
%       (default, jet(numel(unique(clus_label))))
%   idx2nan: indeces of Y to make nan (for plotting)
%       (default, [])
%   indvidx: individual indeces, when data comes from 
%       many different individuals (size: size(Y, 1))
%       (default, [])
%   plottype: plot subclusters (using hierarchical clustering)
%       or just plot error bars (sem or sd)
%       (default, 0)
%   errortype: type of errorbar to plot (@steom, @std)
%       (default, @steom)
%   add_ref_trace: add a reference line (line of the mean across all timepoints)
%       to each plotted cluster
%       (default, 0)
%
% Notes

if ~exist('link_ths', 'var') || isempty(link_ths)
    link_ths = [];
end

if ~exist('axH', 'var') || isempty(axH)
    axH = subplot(1, 1, 1);
    axH.Position = [0 0 1 1];
end

if ~exist('colorvect', 'var') || isempty(colorvect)
    colorvect = jet(numel(unique(clus_label)));
end

if ~exist('y_gap', 'var') || isempty(y_gap)
    y_gap = 4.2;
end

if ~exist('textgate', 'var') || isempty(textgate)
    textgate = 1;
end

if ~exist('textorder', 'var') || isempty(textorder)
    textorder = 1;
end

if ~exist('idx2use4clus', 'var') || isempty(idx2use4clus)
    idx2use4clus = [];
end

if ~exist('colorvect_mean', 'var') || isempty(colorvect_mean)
    colorvect_mean = colorvect;
end

if ~exist('idx2nan', 'var') || isempty(idx2nan)
    idx2nan = [];
end

if ~exist('indvidx', 'var') || isempty(indvidx)
    indvidx = [];
end

if ~exist('plottype', 'var') || isempty(plottype)
    plottype = 0;
end

if ~exist('errortype', 'var') || isempty(errortype)
    errortype = @steom;
end
 
if ~exist('add_ref_trace', 'var') || isempty(add_ref_trace)
    add_ref_trace = 0;
end

% get clusters (remove 0 labels)
subclus = cell(numel(unique(clus_label)), 1);
clus2run = unique(clus_label)';
clus2run(isnan(clus2run)) = [];
clus2run = setdiff(clus2run, 0);

Y_mean_trace_sc = cell(max(clus2run), 1);
Y_mean_trace = nan([max(clus2run), size(Y, 2)]);

% zscore all traces
Y = zscorebigmem(Y);
sub_clus_order = [];

% collect new subclusters (same indexing as in Y)
sub_clus_label = zeros(size(Y, 1), 1);

k = 1;
clus_i = 1;
fileidx_n = nan(max(clus2run), 1);

for i = clus2run
    
    % generate mean Y trace for all Y indeces per cluster
    Y_idx = find(clus_label == i);
    temp_Y = Y(Y_idx, :);
    
    if ~isempty(indvidx)
        fileidx_n(i, 1) = numel(unique(indvidx(Y_idx)));
    end
    
    Y_mean_trace(i, :) = nanmean(temp_Y, 1);

    if ~isempty(link_ths) && ...
            size(temp_Y, 1) > 1 && plottype == 0
        
        if isempty(idx2use4clus)
            
            % get subclusters within cluster
            subclus{i} = hierarchicalClus(...
                zscorebigmem(temp_Y), ...
                link_ths, 'euclidean', 'ward', 3);
            
        else
            
            % use selected timestamps (Y indeces)
            subclus{i} = hierarchicalClus(...
                zscorebigmem(temp_Y(:, idx2use4clus)), ...
                link_ths, 'euclidean', 'ward', 3);
            
        end

        sub_clus_u = unique(subclus{i}.label);

        for j = sub_clus_u'

            % whiten extra traces:
            colorvect_ = colorGradient(colorvect(i, :), [1 1 1], 5);

            % generate mean Y trace for each subcluster
            Y_s_idx = find(subclus{i}.label == j);
            temp_Y_ = temp_Y(Y_s_idx, :);
            Y_mean_trace_sc{i}(j, :) = nanmean(temp_Y_, 1);
            
            % make nan gaps
            if ~isempty(idx2nan)
                Y_mean_trace_sc{i}(j, idx2nan) = nan;
            end
            
            plot(Y_mean_trace_sc{i}(j, :) + (i-1)*y_gap, ...
                'Color', colorvect_(3, :), 'Linewidth', 1)
            
            if i == 1
                hold(axH, 'on')
            end
            
            % Collect new clus label
            sub_clus_label(Y_idx(Y_s_idx), 1) = clus_i;
            sub_clus_order{clus_i, 1} = Y_idx(Y_s_idx);
            clus_i = clus_i + 1;
            
        end
        
    end
    
    % make nan gaps
    if ~isempty(idx2nan)
        Y_mean_trace(i, idx2nan) = nan;
    end
    
    if plottype == 0
        
        % plot mean trace
        plot(Y_mean_trace(i, :) + (i-1)*y_gap, ...
            'Color', colorvect_mean(i, :), 'Linewidth', 2) 
        
    else
        
        % plot mean trace and error bars
        if ~isempty(idx2nan)
            temp_Y(i, idx2nan) = nan;
        end
        
        if size(temp_Y, 1) > 1
            
            idx_nan = isnan(Y_mean_trace(i, :));
            
            if sum(idx_nan)
                
                time_st = 1:length(idx_nan);
                time_st(idx_nan) = nan;
                idx = all(isnan(time_st), 1);
                idy = 1 + cumsum(idx);
                idy(idx) = nan;
                segments_ = unique(idy(~idx));
                
                for j = 1:numel(segments_)
                   
                    shadedErrorBar(time_st(idy == segments_(j)), ...
                        temp_Y(:, idy == segments_(j)) + (i-1)*y_gap, ...
                        {@nanmean, errortype}, 'lineprops', ...
                        {'Color', colorvect_mean(i, :), 'Linewidth', 2}, ...
                        'patchSaturation', 0.4, 'Parent', axH)
                    
                end
                
            else
                
                shadedErrorBar([], temp_Y + (i-1)*y_gap, ...
                    {@nanmean, errortype}, 'lineprops', ...
                    {'Color', colorvect_mean(i, :), 'Linewidth', 2}, ...
                    'patchSaturation', 0.4, 'Parent', axH)
                
            end
            
        else
            
            plot(temp_Y + (i-1)*y_gap, ...
                'Color', colorvect_mean(i, :), ...
                'Linewidth', 2, 'Parent', axH)
            
        end
        
    end
    
    hold(axH, 'on')
    
    if add_ref_trace
        y_ = mean(Y_mean_trace(i, :) + (i-1)*y_gap);
        plot([0 size(Y_mean_trace, 2)], [y_ y_], ...
           '--k', 'Linewidth', 2) 
    end
    
    % add text with info on # of rows, and number of #files(individuals)
    if textgate
        if textorder
            text(axH, -80, (i-1)*y_gap, ...
                [num2str(i), ...
                '(',  num2str(size(temp_Y, 1)), ')', ...
                '(', num2str(fileidx_n(i)), ')']);
        else
            text(axH, -80, (i-1)*y_gap, ...
                [num2str(clus2run(end - k + 1)), ...
                '(', num2str(size(temp_Y, 1)), ')', ...
                '(', num2str(fileidx_n(i)), ')']);
        end
    end
    
    clear temp_Y_ temp_Y
    k = k + 1;
    
end

end
