function batch_plot_sMod_results(Filename, oDir, iparams)
% batch_plot_sMod_results: plots results from stimuli modulation test
%   (distribution of correlations, raw vs null, mod vs non-mod)
%
% Usage:
%   batch_plot_sMod_results(FileName, oDir, iparams)
%
% Args:
%   Filename: name pattern of files to use
%       (default, [])
%   oDir: output directory.
%   iparams: parameters to update
%       (fsuffix: suffix of raw data)
%           (default, '_rawdata.mat')
%       (fmetsuffix: suffix of metadata file)
%           (default, '_metadata.mat')
%       (hbins: bins for histogram)
%           (default, -1:0.01:1)
%       (hbinsp: bins for histogram pvals)
%           (default, 0:0.01:1)
%       (prct2use: bins for histogram)
%           (default, 30)
%       (fdr: bins for histogram)
%           (default, .01)
%       (mccor_method: multiple comparison 
%           correction to use: dep, pdep, bh)
%           see calculate_pval.m
%            (default, 'dep')
%       (dir_depth: depth of directory search)
%           (default, 0)
%
% Notes:
%
% ToDo:
%   add legends to pull plots and stitch ones

% default params
motpar = [];
motpar.fsuffix = '_prosroi.mat';
motpar.fmetsuffix = '_prosmetadata.mat';
motpar.hbins = -1:0.01:1;
motpar.hbinsp = 0:0.01:1.1;
motpar.prct2use = 30;
motpar.fdr = 0.01;
motpar.mccor_method = 'dep';
motpar.dir_depth = 0;

if ~exist('Filename', 'var')
    Filename = [];
end

if ~exist('oDir', 'var') || isempty(oDir)
    oDir = [pwd, filesep, 'smod'];
end

% update variables
if ~exist('iparams', 'var'); iparams = []; end
motpar = loparam_updater(motpar, iparams);

if ~isempty(oDir) && exist('oDir', 'var')
    mkdir(oDir)
end

% define files to use
if motpar.dir_depth == 0
    f2run = rdir(['.', filesep, '*', ...
        motpar.fsuffix]);
elseif motpar.dir_depth == 1
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', motpar.fsuffix]);
else
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', filesep, '*', ...
        motpar.fsuffix]);
end

f2run = {f2run.name}';
[filename, iDir] = split_path(f2run);

if ~isempty(Filename)
    f2run = find(contains(filename, Filename));
    filename = filename(f2run);
    iDir = iDir(f2run);
end

filename = strrep(filename, motpar.fsuffix, '');

fprintf(['Generating plots for ', ...
    num2str(numel(filename)), ' files\n'])

% plot sMod results
for i = 1:numel(filename)

    try
        plotperfly(filename{i}, iDir{i}, oDir, motpar);
    catch
        fprintf('/n ****************************************** /n')
        fprintf(['Failed file: ', filename{i}, '/n'])
        fprintf('****************************************** /n')
    end

end

fprintf('... Done\n')

end

function plotperfly(filename, iDir, oDir, motpar)
% plotperfly: plot results per file
%
% Usage:
%   plotperfly(filename, motpar)
%
% Args:
%   filename: file name to load
%   iDir: directory of filename
%   oDir: output directory
%   motpar: parameter variable

load([iDir, filesep, filename, motpar.fsuffix], ...
    'sMod', 'roi')
load([iDir, filesep, filename, motpar.fmetsuffix], ...
    'wDat')

% 1) ************ correlation related ************
% ************************************************

% correct correlation coefficients and generate pvalues

% load null
if isfield(sMod, 'CC_cs')
    corrcoef_shuffle = sMod.CC_cs;
end

if isfield(sMod, 'CC_as')
    corrcoef_shuffle = sMod.CC_as;
end

if isfield(sMod, 'CC_afs')
    corrcoef_shuffle = sMod.CC_afs;
end

% correct correlations
[corrcoef_stat, corrcoef_raw, corrcoef_shuffle] = ...
    correct_corrcoef(sMod.CC_raw, ...
    corrcoef_shuffle, motpar.prct2use);

% generate pvalues
[pval_raw, pvalc_dep, pvalc_pdep, ...
    pvalc_bh, ~] = calculate_pval(...
    corrcoef_stat, corrcoef_raw, corrcoef_shuffle, ...
    motpar.fdr, [1 1 1 0]);

% find stimuli modulated ROIs
pval_cor = [];
if strcmp(motpar.mccor_method, 'pdep')
    pval_cor = pvalc_pdep;
elseif strcmp(motpar.mccor_method, 'dep')
    pval_cor = pvalc_dep;
elseif strcmp(motpar.mccor_method, 'bh')
    pval_cor = pvalc_bh;                    
elseif strcmp(motpar.mccor_method, 'raw')
    pval_cor = pval_raw;
end
selIdx = pval_cor <= motpar.fdr;

% get histogram of stim mod vs non-stim mod
smod_hist = hist(flat_matrix(corrcoef_raw(selIdx, :)), ...
    motpar.hbins);
nsmod_hist = hist(flat_matrix(corrcoef_raw(~selIdx, :)), ...
    motpar.hbins);
smod_hist = bsxfun(@rdivide, ...
    smod_hist, sum(smod_hist, 2));
nsmod_hist = bsxfun(@rdivide, ...
    nsmod_hist, sum(nsmod_hist, 2));

% generate histograms per ROI
roi_n = size(corrcoef_raw, 1);
for roi_i = 1:roi_n
    
    corrcoef_hist(roi_i, :) = ...
        hist(corrcoef_raw(roi_i, :), motpar.hbins);
    corrcoef_hist_null(roi_i, :) = ...
        hist(corrcoef_shuffle(roi_i, :), motpar.hbins);
    
end
fprintf('\n')

% plot histogram of corrcoef from raw and null
fprintf('Plotting CC and pval\n')

plot_distribution_all(corrcoef_hist, ...
    corrcoef_hist_null, corrcoef_stat, ...
    smod_hist, nsmod_hist, pval_raw, ...
    pval_cor, motpar, filename, oDir)

% 2) ************ correlation related ************
% ************************************************

% get traces from roi
CaRaw = roi.filtered.dfof;

% load signal from ref channel
CaRef = roi.filtered_.df + roi.filtered_.F;

temp_stim = zeros(length(sMod.stim), ...
    size(sMod.stimM, 2));
temp_stim(sMod.tidx2use, :) = sMod.stimM;
CaPred = (zscorebigmem(temp_stim')'...
    *squeeze(mean(sMod.lFilter(:, :, :), 2)))';

tRmode = 0;
stim2use = [];
ITI = wDat.sTime(:, 1);
ITI = round(min(diff(ITI(:))));
dStim = max(wDat.sTime(:, 2) - wDat.sTime(:, 1));
time_range = [-ITI + dStim + 5, ITI - 5];

if contains(wDat.datatype, 'opto')
    tRmode = 1;
    stim2use = [1:3 5:7];
    time_range = [-8 10; -2 5; -2 5; ...
        nan nan; -2 5; -2 5; -2 5];
end

% generate traces per trial
if sum(selIdx)
    
    fprintf(['Processing ', num2str(sum(selIdx)), ' ROIs\n'])
    
    [wDat, CaRaw_zs_tl, CaPred_zs_tl, CaRef_zs_tl, ...
        stim_idx, stim_vect, time_tl, time_tl_rel] = ...
        trace2trial_avg(tRmode, stim2use, time_range, ...
        CaRaw(selIdx, :), CaPred(selIdx, :), ...
        CaRef(selIdx, :), wDat);
    
end

% plot smod traces
if sum(selIdx)
    
    fprintf('Plotting smod ROI responses\n')
    
    plot_med_sig_per_smodroi(CaRef_zs_tl, ...
        CaRaw_zs_tl, corrcoef_stat(selIdx), ...
        time_tl, stim_vect, stim_idx, wDat, ...
        motpar, filename, oDir)
end

% plot coverage videos
if sum(selIdx)
    
    fprintf('Plotting smod ROI coverage\n')
    
    plot_roi_coverage([filename, '_smod'], ...
        [1 1 0], wDat, roi, oDir, find(selIdx)')
    
end

end

function plot_distribution_all(corrcoef_hist, ...
    corrcoef_hist_null, corrcoef_stat, ...
    smod_hist, nsmod_hist, pval_raw, ...
    pval_cor, motpar, filename, oDir)
% plot_distribution_all: distribution of stimuli correlations
%
% Usage:
%   plot_distribution_all(corrcoef_hist, ...
%      corrcoef_hist_null, corrcoef_stat, ...
%      smod_hist, nsmod_hist, pval_raw, ...
%      pval_cor, motpar, filename)
%
% Args:
%   corrcoef_hist: distribution of correlation of raw data
%   corrcoef_hist_null: distribution of correlation of shuffle data
%   corrcoef_stat: selected percentile of corrcoef_raw 
%   smod_hist: correlation of raw data
%   nsmod_hist: correlation of shuffle data
%   pval_raw: pvalues
%   pval_cor: corrected pvalues
%   motpar: parameter variable
%   filename: file name
%   oDir: output directory

% generate histograms
roi_n = size(corrcoef_stat, 1);

% plot histogram
figH = figure('Position', ...
    genfigpos(1, 'center', [1600 700])); 
axH(1) = subplot(2, 4, 1);
axH(2) = subplot(2, 4, 5);

axH(3) = subplot(2, 4, 2);
axH(4) = subplot(2, 4, 6);

axH(5) = subplot(2, 4, 3);
axH(6) = subplot(2, 4, 7);

axH(7) = subplot(2, 4, 4);
axH(8) = subplot(2, 4, 8);

% collect all
corrcoef_hist_all = sum(corrcoef_hist, 1);
corrcoef_hist_null_all = sum(corrcoef_hist_null, 1);

% normalize per column
corrcoef_hist_all = bsxfun(@rdivide, ...
    corrcoef_hist_all, sum(corrcoef_hist_all, 2));
corrcoef_hist_null_all = bsxfun(@rdivide, ...
    corrcoef_hist_null_all, sum(corrcoef_hist_null_all, 2));

corrcoef_hist = bsxfun(@rdivide, ...
    corrcoef_hist, sum(corrcoef_hist, 2));
corrcoef_hist_null = bsxfun(@rdivide, ...
    corrcoef_hist_null, sum(corrcoef_hist_null, 2));

% plot histogram of corrcoef from raw and null
icolormap = [colorGradient([1 1 1], [1 0 1], 15); ...
    colorGradient([1 0 1], [0 1 1], 15)];

[~, idx_order] = sort(corrcoef_stat);

imagesc(motpar.hbins, 1:roi_n, ...
    corrcoef_hist_null(idx_order, :), ...
    'Parent', axH(1))
colormap(axH(1), icolormap)
caxis(axH(1), [0 .12])

imagesc(motpar.hbins, 1:roi_n, ...
    corrcoef_hist(idx_order, :), ...
    'Parent', axH(2))
colormap(axH(2), icolormap)
caxis(axH(2), [0 0.4])

axH(1).YDir = 'normal';
axH(2).YDir = 'normal';
axH(1).Title.String = 'CC Shuffle data';
axH(2).Title.String = 'CC Raw data';
axH(1).YLabel.String = 'ROI number';
axH(2).YLabel.String = 'ROI number'; 
axH(2).XLabel.String = 'pearson correlation';

cbH(1) = colorbar('peer', axH(1)); 
cbH(1).Label.String = 'Probability'; 
cbH(1).Label.FontSize = 10;

cbH(2) = colorbar('peer', axH(2)); 
cbH(2).Label.String = 'Probability'; 
cbH(2).Label.FontSize = 10;

% plot histogram smod vs ~smod
lineH(1) = plot(motpar.hbins, smod_hist, ...
    'k', 'Linewidth', 2, 'Parent', axH(3));
hold(axH(3), 'on')
lineH(2) = plot(motpar.hbins, nsmod_hist, ...
    'Color', [0.5 0.5 0.5], ...
    'Linewidth', 2, 'Parent', axH(3));

axH(3).Title.String = 'smod vs ~smod ROIs';
axH(3).YLabel.String = 'Probability'; 

legend(axH(3), lineH, {'smod', '~smod'}, ...
    'location', 'northwest')

% plot histogram shuffle vs raw
lineH(1) = plot(motpar.hbins, corrcoef_hist_null_all, ...
    'k', 'Linewidth', 2, 'Parent', axH(4));
hold(axH(4), 'on')
lineH(2) = plot(motpar.hbins, corrcoef_hist_all, ...
    'Color', [0.5 0.5 0.5], ...
    'Linewidth', 2, 'Parent', axH(4));

axH(4).Title.String = 'shuffle vs raw';
axH(4).XLabel.String = 'pearson correlation'; 
axH(4).YLabel.String = 'Probability'; 

legend(axH(4), lineH, {'shuffle', 'raw'}, ...
    'location', 'northwest')

% plot histogram shuffle vs raw
pval_raw_hist = hist(pval_raw, motpar.hbinsp);
pval_cor_hist = hist(pval_cor, motpar.hbinsp);
pval_raw_hist = bsxfun(@rdivide, ...
    pval_raw_hist, sum(pval_raw_hist, 2));
pval_cor_hist = bsxfun(@rdivide, ...
    pval_cor_hist, sum(pval_cor_hist, 2));

lineH(1) = plot(motpar.hbinsp, pval_raw_hist, ...
    'k', 'Linewidth', 2, 'Parent', axH(5));
hold(axH(5), 'on')
lineH(2) = plot(motpar.hbinsp, pval_cor_hist, ...
    'Color', [0.5 0.5 0.5], ...
    'Linewidth', 2, 'Parent', axH(5));

legend(axH(5), lineH, {'pval', 'corpval'}, ...
    'location', 'northeast')

plot(motpar.hbinsp, pval_raw_hist, ...
    'k', 'Linewidth', 2, 'Parent', axH(6));
hold(axH(6), 'on')
plot(motpar.hbinsp, pval_cor_hist, ...
    'Color', [0.5 0.5 0.5], ...
    'Linewidth', 2, 'Parent', axH(6));

axH(5).Title.String = 'pvalue';
axH(5).XLim = [0 0.9];
axH(5).XLabel.String = 'pvalue';
axH(5).YLabel.String = 'Probability';

axH(6).Title.String = 'pvalue zoom';
axH(6).XLim = [0 0.1];
axH(6).XLabel.String = 'pvalue';
axH(6).YLabel.String = 'Probability';

legend(axH(6), lineH, {'pval', 'corpval'}, ...
    'location', 'northeast')

plot(pval_raw, corrcoef_stat, ...
    '.k', 'MarkerSize', 20, 'Parent', axH(7));
hold(axH(7), 'on')
plot(pval_raw(pval_cor <= motpar.fdr), ...
    corrcoef_stat(pval_cor <= motpar.fdr), ...
    '.c', 'MarkerSize', 20, 'Parent', axH(7));

plot(pval_cor, corrcoef_stat, ...
    '.k', 'MarkerSize', 20, 'Parent', axH(8));

axH(7).Title.String = ['pvalue vs cc (', ...
    num2str(sum(pval_raw <= motpar.fdr)), ', ', ...
    num2str(sum(pval_raw <= motpar.fdr)*100/roi_n), ')'];
axH(7).XLim = [0 0.05];
axH(7).XLabel.String = 'pvalue';
axH(7).YLabel.String = 'pearson correlation';
line(.01*ones(1, 2), ylim(axH(7)), 'color', 'r', 'Parent', axH(7))

axH(8).Title.String = ['pvalue-cor vs cc (', ...
    num2str(sum(pval_cor <= motpar.fdr)), ', ', ...
    num2str(sum(pval_cor <= motpar.fdr)*100/roi_n), ')'];
axH(8).XLim = [0 0.05];
axH(8).XLabel.String = 'pvalue';
axH(8).YLabel.String = 'pearson correlation';
line(.01*ones(1, 2), ylim(axH(8)), 'color', 'r', 'Parent', axH(8))

fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;

figformat = [1 0 0 0 0 0 0 0 1 0 1];
save_edit_fig_int(axH, figH, oDir, ...
    [filename, '_smod_results'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH)

end

function plot_med_sig_per_smodroi(CaRef_zs_tl, ...
    CaRaw_zs_tl, corrcoef_stat, time_tl, stim_vect, ...
    stim_idx, wDat, motpar, filename, oDir)
% plot_med_sig_per_smodroi: plot traces of smod rois
%
% Usage:
%   plot_med_sig_per_smodroi(CaRef_zs_tl, ...
%       CaRaw_zs_tl, corrcoef_stat, time_tl, stim_vect, ...
%       wDat, motpar, filename)
%
% Args:
%   CaRef_zs_tl: z-scored median traces of reference channel
%   CaRaw_zs_tl: z-scored median traces of raw channel
%   corrcoef_stat: selected percentile of corrcoef_raw 
%   time_tl: timepoints
%   stim_vect: stimulus vector
%   stim_idx: stimulus vector
%   wDat: metadata parameter
%   motpar: parameter variable
%   filename: file name
%   oDir: output directory

[~, idx_order] = sort(corrcoef_stat);
color_vect = cool(numel(unique(stim_idx)));

figH = figure('Position', ...
    genfigpos(1, 'center', [1600 700])); 
axH(1) = subplot(1, 3, 1);
axH(2) = subplot(1, 3, 2);

if numel(idx_order) > 1
    
    imagesc(time_tl, 1:numel(idx_order), ...
        CaRaw_zs_tl(idx_order, :), 'Parent', axH(1))
    colormap(axH(1), 'parula')
    hold(axH(1), 'on')
    caxis(axH(1), [-1 1])

    imagesc(time_tl, 1:numel(idx_order), ...
        CaRef_zs_tl(idx_order, :), 'Parent', axH(2))
    colormap(axH(2), 'parula')
    hold(axH(2), 'on')
    caxis(axH(2), [-1 1])
    
    axH(1).YLabel.String = 'ROI number';
    axH(1).YDir = 'normal';
    axH(2).YDir = 'normal';
    
else
    
    plot(time_tl, CaRaw_zs_tl(idx_order, :), ...
        'Parent', axH(1))
    hold(axH(1), 'on')

    plot(time_tl, CaRef_zs_tl(idx_order, :), ...
        'Parent', axH(2))
    hold(axH(2), 'on')
        
end

axH(1).XLabel.String = 'Time (s)'; 
axH(2).XLabel.String = 'Time (s)';

k = 1;
for i = unique(stim_idx)

    try
        sEtn = [find(stim_vect == i, 1, 'first') ...
            find(stim_vect == i, 1, 'last')];
        linH(k) = line(time_tl(sEtn([1 1])), ylim(axH(1)), ...
           'color', color_vect(i, :), 'linewidth', 2, ...
           'Parent', axH(1));
        line(time_tl(sEtn([2 2])), ylim(axH(1)), ...
           'color', color_vect(i, :), 'linewidth', 2, ...
           'Parent', axH(1))
        legend_str{k} = strrep(wDat.sName{i}, '_', '-');
        k = k + 1;
    end
    
end

lgH = legend(axH(2), linH, legend_str);
lgH.Position = [0.67 0.45 0.27 0.15];

fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;

figformat = [1 0 0 0 0 0 0 0 1 0 1];
save_edit_fig_int(axH, figH, oDir, ...
    [filename, '_smod_roi_results'], figformat, ...
    fitsize, axcolor, figcolor, xyzcolor, ...
    tickgate, [], fontsiz)

close(figH)

end
