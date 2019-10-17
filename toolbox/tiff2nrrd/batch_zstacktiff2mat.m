function batch_zstacktiff2mat(FolderName, FileName, iparams)
% batch_zstacktiff2mat: parse and converts tiff files (in a zstack
% structure) to mat files (files with *Zstack pattern)
%
% Usage:
%   batch_zstacktiff2mat(FolderName, FileName, iparams)
%
% Args:
%   FolderName: folders to use (name)
%   FileName: files to use (name)
%   iparams: parameters
%       (FieldOfView: 768 um, hard coded constant (microscope dependent variable)
%       (FileName: files to use)
%       (suffix: pattern of files to use '_Zstack_')
%       (ch2save: channels to save ([1 2], red and green))
%       (SpMode: type of data (3DxT))
%       (ths: intensity threshold, to detect frames where PMTs went off)
%       (unitres: number of decimals to keep (10^4))
%       (format: nrrd format raw vs gzip ('gzip'))
%       (redo: gate to overwrite old files)
%       (zres: resolution in z axis (1))
%       (pixelsym: pixel size is symmetric or not, 0 = asymmetric, 1 = symmetric)
%       (integrate_flag: number of timepoints to integrate as single timepoint)
%
% Notes:
% FieldOfView needs to be define for each setup (see batch_tiff2mat.m)

zt2m = []; 
zt2m.cDir = pwd;
zt2m.FieldOfView = 768;
zt2m.FileName = []; 
zt2m.suffix = 'Zstack';
zt2m.ch2save = [1 2]; 
zt2m.SpMode = '3DxT'; 
zt2m.ths = 55;
zt2m.fo2reject = {'.', '..', 'colldata', 'BData', 'preprocessed'};
zt2m.unitres = 10^4;
zt2m.format = 'gzip';
zt2m.redo = 0;
zt2m.zres = 1;
zt2m.pixelsym = 0;
zt2m.integrate_flag = 1;

if ~exist('iparams', 'var'); iparams = []; end
zt2m = loparam_updater(zt2m, iparams);

if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('FileName', 'var'); FileName = []; end

% Finding folders and filtering out data that is not selected

fo2run = dir;
fo2run = str2match(FolderName, fo2run);
fo2run = str2rm(zt2m.fo2reject, fo2run);
fo2run = {fo2run.name};

fprintf('Running Tiff2Mat\n');
fprintf(['Running n-folders : ', num2str(numel(fo2run)), '\n'])

% Checking files inside folder

for folder_i = 1:numel(fo2run)
    
    zt2m.Folder2Run = fo2run{folder_i}; 
    cd(fo2run{folder_i}); 
    
    % NameSplitter: select all 'year_filenum_ 'zt2m.suffix' _*_*.tiff' names
    [filename_all, ~, ~] = NameSplitter([], FileName, zt2m.suffix);
    
    filename_u = unique(filename_all); 
    clear filename_all
    
    fprintf(['Running Folder : ', num2str(fo2run{folder_i}), ...
        ' (n-exp-types, ', num2str(numel(filename_u)),'),'])
    
    % running each exptype (u_fname)
    % get unique filenum and file-trials per basename
    for filename_i = 1:numel(filename_u)
        
        [~, flyi_all, trial_all] = NameSplitter(filename_u(filename_i), FileName, zt2m.suffix);
        
        for file_i = unique(flyi_all) 
            
            % get the number of trials per unique filenum and basename            
            zt2m.figH = figure();
            zt2m.AxH(1) = subplot(2, 1, 1);
            zt2m.AxH(2) = subplot(2, 1, 2);
            zt2m.figH.Name = [filename_u{filename_i}, '-', num2str(file_i)];
            triali_all = trial_all(flyi_all == file_i);
            
            if isempty(triali_all)
                
                fprintf(['Does not have any file with flynumber :', num2str(file_i), ' \n'])
                
            else
                
                fprintf([' (n-trials, ', num2str(unique(triali_all)), ')\n']); 
                
                for trial_i = unique(triali_all)
                    
                    % Collapsing all timepoints and z-slices to 1 mat file
                    rootname = [filename_u{filename_i}, '_', ...
                        num2str(file_i), '_Zstack_', num2str(trial_i), '_'];
                    
                    if ~exist([rootname(1:end-1), '.nrrd'], 'file') || zt2m.redo
                        
                        tic;
                        collalltrials(rootname, zt2m);
                        clear rootname;
                        toc
                        
                    else
                        
                        fprintf(['File ', rootname(1:end-1), '.nrrd already exist - skipping \n'])
                        
                    end
                    
                end
                
            end
            
            clear triali_all
            
        end
        
        clear flyi_all trial_i
        
    end
    
    cd(zt2m.cDir)
    fprintf('\n')
    
end

fprintf('\n ********** Done **********\n')

end

function [basename, filenum, repnum] = ...
    NameSplitter(foldername, filename, isuffix)
% NameSplitter: splits name into: basename, file-number and rep/trial-number
% 
% Usage:
%   [basename, flynum, repnum] = ...
%       NameSplitter(foldername, filename, isuffix)
%
% Args:
%   foldername: folders to use (name)
%   filename: filename to use (name)
%   isuffix: suffix to select files
%
% Notes:

if ~exist('foldername', 'var'); foldername = []; end

% hardcoded for current format naming
alltiff = rdir(['*', isuffix, '_*_*.tif']);

alltiff = str2match(filename, alltiff);
alltiff = str2match(foldername, alltiff);
alltiff = {alltiff.name};

basename = cell(1, numel(alltiff)); 
filenum = zeros(1, numel(alltiff)); 
repnum = zeros(1, numel(alltiff));

for FNum = 1:numel(alltiff)
    
    TempS = strsplit2(alltiff{FNum}, '_');
    basename{1, FNum} = TempS{1};
    filenum(1, FNum) = str2double(TempS{2});
    repnum(1, FNum) = str2double(TempS{4});
    
end

end

function collalltrials(rootname, zt2m)
% collalltrials: parse all tiffs with 'rootname*.tif', 
%   average trials, and generate a 3D matrix
%
% Usage:
%   collalltrials(rootname, zt2m)
%
% Args:
%   rootname: suffix of tif files to use
%   zt2m: input params
%       (FieldOfView: 768 um, hard coded constant (microscope dependent variable)
%       (FileName: files to use)
%       (Folder2Run: current folder)
%       (ch2save: channels to save ([1 2], red and green))
%       (SpMode: type of data (3DxT))
%       (ths: intensity threshold)
%       (unitres: number of decimals to keep (10^4))
%       (format: nrrd format raw vs gzip ('gzip'))
%       (zres: resolution in z axis (1))
%
% Notes:

% collapsing files with the same fly and trial number to one mat files
t_fname = rdir([rootname, '0.tif']);
t_fname = {t_fname.name};
t_fname = t_fname{1};
t_fname = strprune(t_fname, '_', 4);
t_num = numel(rdir([rootname, '*.tif']));
clear rootname

% importing tiff files
Data = []; 
fprintf(['n-repts = ', num2str(t_num), '\n']) 

for rep_i = 1:t_num
    
    % load and concatenate all files
    try
        
        [tempdata, ImMeta] = ...
            tiff2mat_scanimage([t_fname, '_', num2str(rep_i - 1)], zt2m.SpMode, 1);
        
        % tiff2mat_scanimage output is always 4D (X, Y, frame, pmt)
        Data = cat(3, Data, tempdata);
        clear tempdata;
        
    end
    
    fprintf('*');
    
    if mod(rep_i, 60) == 0 || rep_i == t_num
        fprintf('\n');
    end
    
end

clear rep_i t_num
fprintf(['Channels imported: ',num2str(ImMeta.ChNum), '\n'])

% update ImMeta, prune Data (ch)
% Selecting channel to save

ch_num = size(Data, 4);

if ch_num > 1
    
    if length(zt2m.ch2save) < 2
        Data = squeeze(Data(:, :, :, zt2m.ch2save)); 
        fprintf([' collecting just channel ', num2str(zt2m.ch2save)])
    end
    
end

clear ch_num
fprintf(['Data original size: ', num2str(size(Data)), '\n'])

% shift distribution by 5 values to the right and make negative values 0
Data = Data + 5; Data(Data < 0) = 0;
ImMeta.RepeatNum = floor(size(Data, 3)/ImMeta.Z); 
ImMeta.FrameNum = ImMeta.Z;

% resize Data
siz = size(Data);

try
    siz(5) = siz(4); 
catch
    siz(5) = 1; 
end

siz(3) = ImMeta.FrameNum; 
siz(4) = ImMeta.RepeatNum;
Data = Data(:, :, 1:siz(4)*siz(3), :);

% integrate signal over time
if zt2m.integrate_flag > 1
    fprintf('integrate signal over time\n')

    new_vol_num = floor(siz(4)/zt2m.integrate_flag);

    preData = Data(:, :, ...
        1:(siz(3)*zt2m.integrate_flag*new_vol_num), :);
    preData = reshape(preData, [siz(1:2), siz(3), ...
        zt2m.integrate_flag, new_vol_num, siz(5)]);
    preData = squeeze(sum(preData, 4));

    Data = preData;
    siz(4) = new_vol_num;
    Data = reshape(Data, [siz(1:2), siz(4)*siz(3), siz(5)]);
    
end

% get average signal
fprintf('Averaging volumes\n')
Data = avbrain(Data, zt2m.ths, siz, zt2m.AxH);

Data = permute(Data, [1 2 4 3]);
ImMeta.FrameNum = size(Data, 4);

fprintf(['Data final size: ',num2str(size(Data)), '\n'])
eval(['Data = ', ImMeta.Imclass, '(Data);'])

premaxRed = max(Data(:, :, :, 1), [], 3); 
try 
    preMaxGreen = max(Data(:, :, :, 2), [], 3);
catch
    preMaxGreen = nan;
end

display(['Max val per channel ', ...
    num2str([max(premaxRed(:)), max(preMaxGreen(:))])])

% generate metadata

[fDat, iDat] = generatemetadata(t_fname, ImMeta, zt2m);

siz = size(Data); 
Data = reshape(Data, [siz(1:2), prod(siz(3:4))]);

% save data as nrrd image
nrrdWriter([t_fname, '.nrrd'], mat2uint16(Data, 0), ...
    iDat.MetaData{3}, [0 0 0], zt2m.format);

% save metadata
save([zt2m.cDir, filesep, zt2m.Folder2Run, filesep, ...
    fDat.FileName, '_metadata.mat'], 'fDat', 'iDat')

clear Data iDat fDat ImMeta

end

function [fDat, iDat] = generatemetadata(t_fname, ImMeta, zt2m)
% generatemetadata: Loading existing metadata, Load matfile asociated with this tiff file
% Usage:
%   [fDat, iDat] = generatemetadata(t_fname, ImMeta, zt2m)
%
% Args:
%   t_fname: folders to use (name)
%   ImMeta: image intensity threshold
%   zt2m: input params
%       (FieldOfView: 768 um, hard coded constant (microscope dependent variable)
%       (Folder2Run: current folder)
%       (SpMode: type of data (3DxT))
%       (ths: intensity threshold)
%       (unitres: number of decimals to keep (10^4))
%       (zres: resolution in z axis (1))
%
% Notes:

try
    load([zt2m.cDir, filesep, zt2m.Folder2Run, ...
        filesep, t_fname, '.mat'], 'fDat', 'iDat');
end

% fDat

fDat.FileName = strprune(t_fname, '_', 3); 
fDat.FolderOrigin = zt2m.cDir;
fDat.FolderTrace = zt2m.cDir((max(strfind(zt2m.cDir, filesep)) + 1):end);
fDat.DataType = zt2m.SpMode; 
fDat.fName = []; 

iDat = struct('MetaData', [], 'ZoomFactor', ImMeta.Zoom);

% pixel size (symmetric or not)
if zt2m.pixelsym == 0
    
    iDat.MetaData = {'voxelsize', 'x y z', ...
        round([zt2m.FieldOfView/(ImMeta.Zoom*ImMeta.X), ...
        zt2m.FieldOfView/(ImMeta.Zoom*ImMeta.Y), ...
        zt2m.zres]*zt2m.unitres)/zt2m.unitres};
    
else
    
    iDat.MetaData = {'voxelsize', 'x y z', ...
        round([zt2m.FieldOfView/(ImMeta.Zoom*ImMeta.X), ...
        zt2m.FieldOfView/(ImMeta.Zoom*ImMeta.X), ...
        zt2m.zres]*zt2m.unitres)/zt2m.unitres};
    
end

iDat.FrameN = ImMeta.FrameNum; 
iDat.StackN = ImMeta.RepeatNum;
iDat.FrameSize = [ImMeta.Y, ImMeta.X]; %[lines, pixels per line, slices] [y, x, z]
iDat.Power = ImMeta.Power; 
iDat.f_ths = zt2m.ths;

end

function avgim = avbrain(im, im_ths, isiz, axH)
% avbrain: get the average image of green and red channel
%
% Usage:
%   avgIm = avbrain(flyname, segment_n, zt2m)
%
% Args:
%   im: folders to use (name)
%   im_ths: image intensity threshold
%   isiz: image size
%   axH: axis handle
%
% Notes:

% get max per frame (both channels)
maxpertime_r = squeeze(max(max(im(:, :, :, 1), [], 1), [], 2));
if isiz(end) > 1
    maxpertime_g = squeeze(max(max(im(:, :, :, 2), [], 1), [], 2)); 
end

im = reshape(im, isiz);

maxpertime_r = reshape(maxpertime_r, isiz([3 4]));
if isiz(end) > 1
    maxpertime_g = reshape(maxpertime_g, isiz([3 4]));
end

% remove flyback planes
im(:, :, 1, :, :) = [];
maxpertime_r(1, :) = [];
try maxpertime_g(1, :) = []; end

% plot max's
plot(maxpertime_r(:), 'r', 'Parent', axH(2));
try plot(maxpertime_g(:), 'g', 'Parent', axH(1)); end
hold(axH(1), 'on'); hold(axH(2), 'on')

% generate average brain
avgim = inf(size(im, 1), size(im, 2), size(im, 3), size(im, 5));

for z = 1:size(im, 3)
    
    if sum((maxpertime_r(z, :)) > im_ths) >= 2
        avgim(:, :, z, 1) = mean(squeeze(im(:, :, z, ...
            (maxpertime_r(z, :)) > im_ths, 1)), 3);
    end
    
    if isiz(end) > 1
        if sum((maxpertime_g(z, :)) > im_ths) >= 2
            avgim(:, :, z, 2) = mean(squeeze(im(:, :, z, ...
                (maxpertime_g(z, :)) > im_ths, 2)), 3);
        end
    end
    
end

% find Inf planes and debug them

M = squeeze(max(max(avgim, [], 1), [], 2));

if sum(isinf(M(:))) ~= 0 || (sum(M(:) > 2060) ~= 0)
    
    % for dealing with this exceptions use:
    
    edit zstacktiff2mat_maninput.m
    keyboard
    
end

end