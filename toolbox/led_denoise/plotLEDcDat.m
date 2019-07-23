function plotLEDcDat(FolderName, FileName)
% plotLEDcDat: plot results of LED correction for selected folders and
% files
%
% Usage:
%   plotLEDcDat(FolderName, FileName)
%
% Args:
%   FolderName: Folder name to load
%   FileName: File name to load

% Default params
if ~exist('FolderName','var'); FolderName = []; end
if ~exist('FileName','var'); FileName = []; end

% Selecting folders
ipars.cDir = pwd;
ipars.file2reject = 'Zstack';
ipars.filesuffix = '_metadata';
ipars.folder2reject = {'.', '..', 'preprocessed', 'BData'};

fol2run = dir;
fol2run = str2match(FolderName, fol2run);
fol2run = str2rm(ipars.folder2reject, fol2run);
fol2run = {fol2run.name};
fprintf(['Running n-folders : ', num2str(numel(fol2run)), '\n'])

for i = 1:numel(fol2run)
    
    fprintf(['Running folder : ', fol2run{i}]); 
    cd(fol2run{i});
    plotperfolder(FileName, ipars);
    cd(ipars.cDir)
    
end

fprintf('... Done\n')

end

function plotperfolder(FileName, ipars)
% plotperfolder: find files to plot
%
% Usage:
%   runperfolder(FileName, params)
%
% Args:
%   FileName: image related field
%   ipars: internal parameters

sep = filesep;
f2run = rdir(['.', sep, '*', ipars.filesuffix, '.mat']);
f2run = str2match(FileName, f2run);
f2run = str2rm(ipars.file2reject, f2run);
f2run = {f2run.name};
fprintf([' , files (', num2str(numel(f2run)), ') '])

for file2run = 1:numel(f2run)
    
    fprintf('*')
    
    try
        
        load(f2run{file2run}, 'cDat')
        axH = figuregenerator(cDat.Iter, ...
            strrep(f2run{file2run}, ['.', sep], ''));
        
        for Iter = 1:cDat.Iter
            resultplotter(cDat, axH([Iter:cDat.Iter:cDat.Iter*3]), Iter)
        end
        
    catch
        display(f2run{file2run})
    end
    
    clear cDat Axeshandle
    
end

fprintf('Done\n')

end

function axH = figuregenerator(iternum, fname)

if iternum == 1
    figure('position', [220 633 1440 331], 'Name', fname)
elseif iternum == 2
    figure('position', [248 460 1440 623], 'Name', fname)
elseif iternum == 3
    figure('position', [248 138 1440 945], 'Name', fname)
end

k = 1;

for i = 1:iternum
    for ii = 1:3
    	eval(['axH(', num2str(i),', ', num2str(ii), ...
            ') = subplot(', num2str(iternum), ', 3, ', num2str(k), ');'])
        k = k + 1;
    end
end

end

function resultplotter(cDat, axH, Iter)

% LTA (light triggered avegare)
cDat.Iter = Iter;
plot(cDat.LTM(cDat.Iter, :), 'b', 'Parent', axH(1)); hold(axH(1), 'on'); 
plot(cDat.LTA(cDat.Iter, :), 'r', 'Parent', axH(1))
plot(cDat.LTMCor(cDat.Iter, :), 'c', 'Parent', axH(1)); 
plot(cDat.LTACor(cDat.Iter, :), 'm', 'Parent', axH(1))
title(axH(1), 'Light triggered avegare')
xlabel(axH(1), 'Time (pixels units)')
ylabel(axH(1), 'Fluorescence (a.u)')
xlim(axH(1), [1 size(cDat.LTACor, 2)])
ylim(axH(1), [min([min(cDat.LTA)*1.1, -2]) max([max(cDat.LTA)*1.1, 2])])

% 2D Pre
imagesc(cDat.XYprojRaw{cDat.Iter}, 'Parent', axH(2))
title(axH(2), 'Raw movie')
xlabel(axH(2), 'Time (frames)')
ylabel(axH(2), 'YX projection (pixels)')
caxis(axH(2), [0 max(cDat.XYprojRaw{cDat.Iter}(:))])
colorbar(axH(2))

% 2D Post
imagesc(cDat.XYprojRawCor{cDat.Iter}, 'Parent', axH(3))
title(axH(3), 'Corrected movie')
xlabel(axH(3), 'Time (frames)')
ylabel(axH(3), 'YX projection (pixels)')
colorbar(axH(3))
caxis(axH(3), [0 max(cDat.XYprojRawCor{cDat.Iter}(:))])

end
