function batch_plotShiftPerStack(FolderName, FileName, iparams)
% batch_plotShiftPerStack: plots shifts for timeseries and stitching
%
% Usage:
%   batch_plotShiftPerStack(iparams)
%
% Args:
% Args:
%   FolderName: name of folders to load
%   FileName: name of files to load
%   iparams: parameters to update
%       (fsuffix: suffix of raw data)
%           (default, '_rawdata.mat')
%       (metsuffix: suffix of metadata)
%           (default, '_metadata.mat')
%       (fo2reject: folders to reject)
%           (default, {'.', '..', 'preprocessed', 'BData'})
%       (fi2reject: files to reject)
%           (default, [])
%       (colorm: color for shift trace)
%           (default,[1 0 0; 0 0 0; 0 0 1; 0 0.5 0])
%       (cgate: gate to plot timeseries motion correction)
%           (default, 1)
%       (sgate: gate to plot stitching correction)
%           (default, 0)
%       (ctype: plot each stack: 0 == separately, 1 == together)
%           (default, 0)
%       (stype: plot each stack: 0 == separately, 1 == together)
%           (default, 1)
%       (screen: screen to use for plotting)
%           (default, 1)
%       (oDir: output folder to save figures)
%           (default, [])
%
% Notes:
%
% ToDo:
% add legends

motpar = [];
motpar.fsuffix = '_rawdata.mat';
motpar.metsuffix = '_metadata.mat';
motpar.fo2reject = {'.', '..', 'preprocessed', 'BData'};
motpar.fi2reject = []; 
motpar.colorm = [1 0 0; 0 0 0; 0 0 1; 0 0.5 0];
motpar.cgate = 1;
motpar.sgate = 0;
motpar.ctype = 0;
motpar.stype = 1;
motpar.screen = 1;
motpar.oDir = [];

if ~exist('FileName', 'var'); FileName = []; end
if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('iparams', 'var'); iparams = []; end
motpar = loparam_updater(motpar, iparams);

% Selecting folders
cDir = pwd;
f2run = dir;
f2run = str2match(FolderName, f2run);
f2run = str2rm(motpar.fo2reject, f2run);
f2run = {f2run.name};

fprintf('Reading motion from metadata files\n')
fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

% plot
if motpar.cgate
    for i = 1:numel(f2run)
        
        fprintf(['Running folder : ', f2run{i}, '\n']); 
        cd(f2run{i});
        runperfolder(FileName, motpar);
        cd(cDir)
        
    end
end

if motpar.sgate
    runperfly_st(FolderName, FileName, motpar);
end

fprintf('... Done\n')

end

function runperfolder(FileName, motpar)

% Plot all flies per folder
[f2plot, ~] = suffixgen(FileName, motpar.fsuffix, motpar.fi2reject);
f2plot = unique(f2plot);
fprintf(['Plotting n-flies : ', num2str(numel(f2plot)), '\n'])

for fly_i = 1:numel(f2plot)
	runperfly_mc(f2plot{fly_i}, motpar)
end

end

function [filename, repnum] = suffixgen(basename, fsuffix, fi2reject)
% suffixgen: select input files
% 
% Usage:
%   [filename, repnum] = suffixgen(basename, fsuffix, fi2reject)
%
% Args:
%   basename: pattern to select input files
%   fsuffix: suffix to find matching files
%   fi2reject: file names to reject

f2run = rdir(['.', filesep, '*', fsuffix]);
f2run = str2rm(fi2reject, f2run);
f2run = str2match(basename, f2run);
f2run = {f2run.name};
f2run = strrep(f2run, ['.', filesep], '');
filename = cell(1, numel(f2run));  
repnum = zeros(1, numel(f2run));

for FNum = 1:numel(f2run)
    
    TempS = strsplit2(f2run{FNum}, '_');
    filename{1, FNum} = [TempS{1}, '_', TempS{2}];
    repnum(1, FNum) = str2double(TempS{3});
    
end

end

function runperfly_st(FolderName, FileName, motpar)
% runperfly_mc: plot stitching motion per inputfile
% 
% Usage:
%   runperfly_st(motpar)
%
% Args:
%   FolderName: pattern to select input files
%   FileName: pattern to select input files
%   motpar: parameters

f2run = rdir(['.', filesep, '*', filesep, '*', motpar.fsuffix]);
f2run = str2rm(motpar.fi2reject, f2run);
f2run = str2match(FolderName, f2run);
f2run = str2match(FileName, f2run);
f2run = {f2run.name};

if motpar.stype
    
    motpar.figH = figure('position', genfigpos(1, 'center', [578 700]));
    motpar.AxH(1) = subplot(311);
    motpar.AxH(2) = subplot(312);
    motpar.AxH(3) = subplot(313);
    motpar.colorm = jet(numel(f2run));
    
end

for i = 1:numel(f2run)
    
    fprintf('*')
    
    if ~isempty(strfind(motpar.fsuffix, 'prosdata'))
        
        load(strrep(f2run{i}, motpar.fsuffix, motpar.metsuffix), 'wDat')
        
        if ~motpar.stype
            
            motpar.figH = figure('position', genfigpos(1, 'center', [578 700]), ...
                'Name', strrep(strrep(f2run{f2r}, '_rawdata', ''), ['.', filesep], ''));
            motpar.AxH(1) = subplot(311);
            motpar.AxH(2) = subplot(312);
            motpar.AxH(3) = subplot(313);
            motpar.colorm(i, :) = [0 0 0];
            
        end
        
        ploteachtrace_st(wDat.Zstitch.Xshift*wDat.XYZres{2}(1), ...
            wDat.Zstitch.Yshift*wDat.XYZres{2}(2), ...
            wDat.Zstitch.Zshift(:, 1)*wDat.XYZres{2}(3), ...
            motpar.colorm(i, :), motpar.AxH)
    
    end
    
end

fprintf('Done\n')

end

function runperfly_mc(inputfile_str, motpar)
% runperfly_mc: plot timeseries motion per inputfile
% 
% Usage:
%   runperfly_mc(inputfile_str, motpar)
%
% Args:
%   inputfile_str: pattern to select input files
%   motpar: parameters

f2run = rdir(['.', filesep, '*', motpar.fsuffix]);
f2run = str2match(inputfile_str, f2run);
f2run = {f2run.name};

for f2r = 1:numel(f2run)
    
    fprintf('*')
    
    if ~isempty(strfind(motpar.fsuffix, 'rawdata'))
        
        load(strrep(f2run{f2r}, motpar.fsuffix, motpar.metsuffix), 'mcDat', 'iDat')
        
        if ~motpar.ctype || f2r == 1
            
            motpar.figH = figure('position', genfigpos(1, 'center', [578 700]), ...
                'Name', strrep(strrep(strrep(f2run{f2r}, '_rawdata', ''), ...
                    ['.', filesep], ''), '.mat', ''));
            
            if size(mcDat.rigid, 1) == 3
                motpar.AxH(1) = subplot(311);
                motpar.AxH(2) = subplot(312);
                motpar.AxH(3) = subplot(313);
            else
                motpar.AxH(1) = subplot(211);
                motpar.AxH(2) = subplot(212);  
            end
            
        end
        
        if ~isempty(mcDat.rigids)
            
            ploteachtrace_mc(mcDat.rigids(1, :)*iDat.MetaData{3}(1), ...
                mcDat.rigids(2, :)*iDat.MetaData{3}(2), ...
                mcDat.rigids(3, :)*iDat.MetaData{3}(3), ...
                motpar.colorm(2, :), motpar.AxH)
        end
        
        ploteachtrace_mc(mcDat.rigid(1, :)*iDat.MetaData{3}(1), ...
            mcDat.rigid(2, :)*iDat.MetaData{3}(2), ...
            mcDat.rigid(3, :)*iDat.MetaData{3}(3), ...
            motpar.colorm(1, :), motpar.AxH)
    
    elseif ~isempty(strfind(motpar.fsuffix, 'prosdata'))
        
        load(strrep(f2run{f2r}, motpar.fsuffix, motpar.metsuffix), 'wDat')
        
        if motpar.ctype
            
            motpar.figH = figure(...
                'position', genfigpos(1, 'center', [578 700]), ...
                'Name', strrep(strrep(strrep(f2run{f2r}, '_rawdata', ''), ...
                    ['.', filesep], ''), '.mat', ''));
            
            if isfield(wDat.MotCor, 'Zshift')
                motpar.AxH(1) = subplot(311);
                motpar.AxH(2) = subplot(312); 
                motpar.AxH(3) = subplot(313);
            else
                motpar.AxH(1) = subplot(211);
                motpar.AxH(2) = subplot(212);
            end
            
            motpar.colorm = jet(size(wDat.MotCor.sXshift, 2));
            
        end
        
        for tIdx = 1:size(wDat.MotCor.sXshift, 2)
            
            if ~motpar.ctype
                
                motpar.figH = figure(...
                    'position', genfigpos(1, 'center', [578 700]), ...
                    'Name', strrep(strrep(strrep(f2run{f2r}, '_rawdata', ''), ...
                    ['.', filesep], ''), '.mat', ''));
                
                if isfield(wDat.MotCor, 'Zshift')
                    motpar.AxH(1) = subplot(311);
                    motpar.AxH(2) = subplot(312);
                    motpar.AxH(3) = subplot(313);
                else
                    motpar.AxH(1) = subplot(211);
                    motpar.AxH(2) = subplot(212);
                end
                
                motpar.colorm(tIdx, :) = [0 0 0];
                
            end
            
            ploteachtrace_mc(wDat.MotCor.sXshift(:, tIdx)*wDat.XYZres{2}(1), ...
                wDat.MotCor.sYshift(:, tIdx)*wDat.XYZres{2}(2), ...
                wDat.MotCor.Zshift(:, tIdx)*wDat.XYZres{2}(3), ...
                motpar.colorm(tIdx, :), motpar.AxH)
            
            if ~motpar.ctype
                figEdit(motpar.AxH, motpar.figH);
            end
            
        end
        
    end
    
    figEdit(motpar.AxH, motpar.figH);
    
    if ~isempty(motpar.oDir)
        
        if ~exist(motpar.oDir, 'dir')
            mkdir(motpar.oDir)
        end
        
        figname = strrep(strrep(f2run{f2r}, motpar.fsuffix, ''), ['.', filesep], '');
        figformat = [0 0 0 0 0 0 0 0 1];
        resolution_ = '-r300';
        
        savefig_int(motpar.figH, motpar.oDir, figname, figformat, resolution_);
        close(motpar.figH)
        
    end
    
end

fprintf('\n')

end

function ploteachtrace_mc(shifts_x, shifts_y, shifts_z, intCorvect, axH)
% ploteachtrace_st: plots timeseries shifts_x, shifts_y, shifts_z on axes axH
% 
% Usage:
%   ploteachtrace_st(shifts_x, shifts_y, shifts_z, intCorvect, ax)
%
% Args:
%   shifts_*: shifts on X, Y or Z axis
%   intCorvect: color vector
%   axH: axes handle

ivect = {'shifts_x', 'shifts_y', 'shifts_z'};
axS = {'X', 'Y', 'Z'};

for i = 1:numel(axH)
    
    eval(['plot(', ivect{i},', ''color'', intCorvect, ''linewidth'', 2, ''Parent'', axH(i));']);
    hold(axH(i), 'on');
    title(axH(i), ['displacements along ', axS{i}]);
    xlabel(axH(i), 'time (stacks)');
    ylabel(axH(i), 'motion um')
    
end

end

function ploteachtrace_st(shifts_x, shifts_y, shifts_z, intCorvect, axH)
% ploteachtrace_st: plots stitching shifts_x, shifts_y, shifts_z on axes axH
% 
% Usage:
%   ploteachtrace_st(shifts_x, shifts_y, shifts_z, intCorvect, ax)
%
% Args:
%   shifts_*: shifts on X, Y or Z axis
%   intCorvect: color vector
%   axH: axes handle

ivect = {'shifts_x', 'shifts_y', 'shifts_z'};
axS = {'X', 'Y', 'Z'};

for i = 1:numel(axH)
    
    eval(['plot(', ivect{i},', ''color'', intCorvect, ''linewidth'', 2, ''Parent'', axH(i));']);
    hold(axH(i), 'on');
    title(axH(i), ['displacements along ', axS{i}]);
    xlabel(axH(i), 'Z (stacks)');
    ylabel(axH(i), 'motion um')
    
end

end