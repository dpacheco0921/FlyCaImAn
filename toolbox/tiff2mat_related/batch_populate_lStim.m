function batch_populate_lStim(FolderName, FileName, iparams)
% batch_populate_lStim: populates lStim with missing fields
%
% Usage:
%   batch_populate_lStim(FolderName, FileName, iparams)
%
% Args:
%   FolderName: name of folders to load
%   FileName: name of files to load
%   iparams: parameters to update
%       (cDir: current directory)
%       (fo2reject: folders to reject)
%       (fi2reject: files to reject)
%       (fsuffix: suffix of files with Data variable)
%           (default, '_rawdata.mat')
%       (fmetsuffix: suffix of metadata file)
%           (default, '_metadata.mat') 
% 
% Notes

% default params
metpars.cDir = pwd;
metpars.fo2reject = {'.', '..', 'preprocessed', 'BData'};
metpars.fi2reject = [];
metpars.fsuffix = '_rawdata.mat';
metpars.fmetsuffix = '_metadata.mat';

% update variables
if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('FileName', 'var'); FileName = []; end
if ~exist('iparams', 'var'); iparams = []; end
metpars = loparam_updater(metpars, iparams);

% find folders
fo2run = dir;
fo2run = str2match(FolderName, fo2run);
fo2run = str2rm(metpars.fo2reject, fo2run);
fo2run = {fo2run.name};

fprintf(['Running n-folders : ', num2str(numel(fo2run)), '\n'])

for i = 1:numel(fo2run)
    
    fprintf(['Running folder : ', fo2run{i}, '\n']);
    cd(fo2run{i}); 
    runperfolder(FileName, metpars);
    cd(metpars.cDir)
    
end

fprintf('... Done\n')

end

function runperfolder(fname, metpars)
% runperfolder: function collects timestamps and stimulus from bin file 
%   with "metpars.fsuffix" suffix and fDat info
%
% Usage:
%   runperfolder(fname, foname, metpars)
%
% Args:
%   fname: file name template string
%   metpars: parameters

% directory params
fi2run = rdir(['.', filesep, '*', metpars.fsuffix]);
fname = addsuffix(fname, metpars.fsuffix);
fi2run = str2match(fname, fi2run);
fi2run = str2rm(metpars.fi2reject, fi2run);
fi2run = {fi2run.name};
fi2run = strrep(fi2run, ['.', filesep], '');
fi2run = strrep(fi2run, '_rawdata.mat', '');

for F2R_idx = 1:numel(fi2run)
    
    display(['Running file : ', fi2run{F2R_idx}])   
    runperfile(fi2run{F2R_idx}, metpars)
    
end

fprintf('****** Done ******\n')

end

function runperfile(filename, metpars)
% runperfolder: function that compiles 
%   missing fields of lStim
%
% Usage:
%   runperfolder(filename, ip)
%
% Args:
%   filename: name of file to load
%   metpars: parameters

% Load main variables 'lStim', 'sDat'
load(['.', filesep, filename, '_metadata.mat'], ...
    'lStim')
load(['.', filesep, filename, '.mat'], 'sDat');

lStim.fs = sDat.fs;
lStim.trialn = sDat.trial;
lStim.fStrain = sDat.fStrain;
lStim.channels = sDat.channels;

lStim.lstEn = optostim_init_end(lStim.trace, ...
    sDat);

% collect extra metadata
lStim.sPars.freq = ...
    repmat(sDat.freq, [size(lStim.lstEn, 1), 1]); 
lStim.sPars.width = ...
    repmat(sDat.width, [size(lStim.lstEn, 1), 1]); 
lStim.sPars.int = ...
    repmat(sDat.intensity, [size(lStim.lstEn, 1), 1]); 
lStim.sPars.sr = ...
    repmat(sDat.fs, [size(lStim.lstEn, 1), 1]); 
lStim.sPars.basPre = ...
    repmat(sDat.silencePre, [size(lStim.lstEn, 1), 1]); 
lStim.sPars.basPost = ...
    repmat(sDat.silencePost, [size(lStim.lstEn, 1), 1]);

lStim.sPars.order = 1:size(lStim.lstEn, 1);

if isfield(sDat, 'stimFileName')
    lStim.sPars.name = ...
        repmat({sDat.stimFileName}, [size(lStim.lstEn, 1), 1]);
else
    lStim.sPars.name = ...
        repmat({'OPTO'}, [size(lStim.lstEn, 1), 1]);                
end

save(['.', filesep, filename, '_metadata.mat'], ...
    'lStim', '-append')

end
