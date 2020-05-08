function batch_plot_stim_ln_of_motor_results(...
    Foldername, Filename, oDir, iparams)
% batch_plot_stim_ln_of_motor_results: plots results of stimuli modulation of
%   motor variables
%
% Usage:
%   batch_plot_stim_ln_of_motor_results(...
%       Foldername, FileName, oDir, iparams)
%
% Args:
%   Foldername: name pattern of folders to use
%       (default, [])
%   Filename: name pattern of files to use
%       (default, [])
%   oDir: output directory.
%   iparams: parameters to update
%       (fmetsuffix: suffix of metadata file)
%           (default, '_metadata.mat')
%       (hbins: bins for histogram)
%           (default, -1:0.01:1)
%       (hbinsp: bins for histogram pvals)
%           (default, 0:0.01:1)
%       (fdr: false discovery rate)
%            (default, 0.01)
%       (mccor_method: multiple comparison correction to use: dep, pdep, bh)
%            (default, 'dep')
%       (dir_depth: depth of directory search)
%           (default, 0)
%
% Notes:

% default params
motpar = [];
motpar.fmetsuffix = '_prosmetadata.mat';
motpar.hbins = -1:0.01:1;
motpar.hbinsp = 0:0.1:1.1;
motpar.prct2use = 100;
motpar.fdr = 0.01;
motpar.mccor_method = 'dep';
motpar.dir_depth = 0;
motpar.coeffname = 'pearson correlation';
motpar.coeffrange_im = {[0 .12], [0 .4]};

if ~exist('Filename', 'var')
    Filename = [];
end

if ~exist('oDir', 'var') || isempty(oDir)
    oDir = [pwd, filesep, 'sti_ln_mot'];
end

% update variables
if ~exist('iparams', 'var'); iparams = []; end
motpar = loparam_updater(motpar, iparams);

if ~exist(oDir, 'dir')
    mkdir(oDir)
end

% define files to use
if motpar.dir_depth == 0
    f2run = rdir(['.', filesep, '*', ...
        motpar.fmetsuffix]);
elseif motpar.dir_depth == 1
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', motpar.fmetsuffix]);
else
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', filesep, '*', ...
        motpar.fmetsuffix]);
end

f2run = {f2run.name}';
[filename, iDir] = split_path(f2run);

if ~isempty(Foldername)
    f2run = find(contains(iDir, Foldername));
    filename = filename(f2run);
    iDir = iDir(f2run);
end

if ~isempty(Filename)
    f2run = find(contains(filename, Filename));
    filename = filename(f2run);
    iDir = iDir(f2run);
end

filename = strrep(filename, motpar.fmetsuffix, '');

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

load([iDir, filesep, filename, motpar.fmetsuffix], ...
    'wDat', 'sti_ln_mot')

% 1) ************ correlation related ************
% ************************************************

% correct correlation coefficients and generate pvalues

% load null correlation distribution eVar
if isfield(sti_ln_mot, 'eVar_cs')
    eVar_shuffle = sti_ln_mot.eVar_cs;
end

if isfield(sti_ln_mot, 'eVar_as')
    eVar_shuffle = sti_ln_mot.eVar_as;
end

if isfield(sti_ln_mot, 'eVar_afs')
    eVar_shuffle = sti_ln_mot.eVar_afs;
end

% correct correlations
[eVar_stat, eVar_raw, eVar_shuffle] = ...
    correct_corrcoef(sti_ln_mot.eVar, ...
    eVar_shuffle, 1);

% generate pvalues
[pval_raw, pvalc_dep, pvalc_pdep, ...
    pvalc_bh, ~] = calculate_pval(...
    eVar_stat, eVar_raw, eVar_shuffle, ...
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
smod_hist = hist(flat_matrix(eVar_raw(selIdx, :)), ...
    motpar.hbins);
nsmod_hist = hist(flat_matrix(eVar_raw(~selIdx, :)), ...
    motpar.hbins);
smod_hist = bsxfun(@rdivide, ...
    smod_hist, sum(smod_hist, 2));
nsmod_hist = bsxfun(@rdivide, ...
    nsmod_hist, sum(nsmod_hist, 2));

% generate histograms per ROI
roi_n = size(eVar_raw, 1);
for roi_i = 1:roi_n
    
    eVar_hist(roi_i, :) = ...
        hist(eVar_raw(roi_i, :), motpar.hbins);
    eVar_hist_null(roi_i, :) = ...
        hist(eVar_shuffle(roi_i, :), motpar.hbins);
    
end
fprintf('\n')

% plot histogram of corrcoef from raw and null
fprintf('Plotting CC and pval\n')

plot_distribution_all(eVar_hist, ...
    eVar_hist_null, eVar_stat, ...
    smod_hist, nsmod_hist, pval_raw, ...
    pval_cor, motpar, filename, oDir)

% plot traces of predicted and raw motvar
var2plot = find(pval_cor <= motpar.fdr);

if ~isempty(var2plot)
    
    plot_traces_of_modulated_vars(var2plot, ...
        wDat, sti_ln_mot, filename, oDir)
    
end

end

function plot_traces_of_modulated_vars(var2plot, ...
    wDat, sti_ln_mot, filename, oDir)
% plot_traces_of_modulated_vars: plot traces of raw and predicted motor
%   variables
%
% Usage:
%   plot_traces_of_modulated_vars(var2plot, ...
%      wDat, sti_ln_mot, oDir)
%
% Args:
%   var2plot: indeces of motor variables that are significantly modulated
%   wDat: metadata structure
%   sti_ln_mot: motor modulation structure
%   oDir: output directory

oDir_ = [oDir, filesep, filename,  'stim_mod'];
mkdir(oDir_)

fitsize = 0;
axcolor = 'none';
figcolor = 'none';
xyzcolor = 'k';
tickgate = 'on';
fontsiz = 10;
figformat = [1 0 0 0 0 0 0 0 1 0 1];

Y_pred = get_predicted_signal(sti_ln_mot.train_idx, ...
    sti_ln_mot.lFilter, sti_ln_mot.stimM);

for i = 1:numel(var2plot)
    figH = figure('Position', genfigpos(1, 'center', [1000 300])); 
    axH(1) = subplot(1, 1, 1);

    plot(wDat.fTime, zscorebigmem(sti_ln_mot.Y(var2plot(i), :)), ...
        'k', 'Linewidth', 2, 'Parent', axH(1));
    hold(axH(1), 'on')
    plot(wDat.fTime, Y_pred(var2plot(i), :), ...
        'r', 'Linewidth', 2, 'Parent', axH(1));

    axH(1).XLim = wDat.fTime([1 end]);
    axH(1).XLabel.String = 'Time (s)';
    axH(1).YLabel.String = 'Signal (SD)';

    save_edit_fig_int(axH, figH, oDir_, ...
        ['motvar_', num2str(var2plot(i))], figformat, ...
        fitsize, axcolor, figcolor, xyzcolor, ...
        tickgate, [], fontsiz)

    close(figH)
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
caxis(axH(1), motpar.coeffrange_im{1})

imagesc(motpar.hbins, 1:roi_n, ...
    corrcoef_hist(idx_order, :), ...
    'Parent', axH(2))
colormap(axH(2), icolormap)
caxis(axH(2), motpar.coeffrange_im{2})

axH(1).YDir = 'normal';
axH(2).YDir = 'normal';
axH(1).Title.String = 'CC Shuffle data';
axH(2).Title.String = 'CC Raw data';
axH(1).YLabel.String = 'ROI number';
axH(2).YLabel.String = 'ROI number'; 
axH(2).XLabel.String = motpar.coeffname;
axH(1).XLim = [-0.2 0.2];
axH(2).XLim = [-0.2 0.2];

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
axH(3).XLim = [-1 1];

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
axH(4).XLabel.String = motpar.coeffname; 
axH(4).YLabel.String = 'Probability'; 
axH(4).XLim = [-1 1];

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
    num2str(sum(pval_raw <= motpar.fdr)*100/roi_n), ') (#, %)'];
axH(7).XLim = [0 0.05];
axH(7).XLabel.String = 'pvalue';
axH(7).YLabel.String = motpar.coeffname;
line(.01*ones(1, 2), ylim(axH(7)), 'color', 'r', 'Parent', axH(7))

axH(8).Title.String = ['pvalue-cor vs cc (', ...
    num2str(sum(pval_cor <= motpar.fdr)), ', ', ...
    num2str(sum(pval_cor <= motpar.fdr)*100/roi_n), ') (#, %)'];
axH(8).XLim = [0 0.05];
axH(8).XLabel.String = 'pvalue';
axH(8).YLabel.String = motpar.coeffname;
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
