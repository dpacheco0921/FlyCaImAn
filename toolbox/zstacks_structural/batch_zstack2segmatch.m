function batch_zstack2segmatch(FolderName, FileName, iparams)
% batch_zstack2segmatch: generates zstacks used for segment to fly and 
%   fly to atlas registration. Whole fly image has the pattern (_w) and
%   fly segment image has the pattern (_s). It also does some preprocessing
%   of images (see notes).
%
% Usage:
%   batch_zstack2segmatch(FolderName, FileName, iparams)
%
%   FolderName: Folder name to load
%   FileName: File name to load
%   ip: parameters to update
%       (cDir: current directory)
%       (fisuffix: file suffix)
%       (fo2reject: folder to reject)
%       (redo: redo)
%           (0, default)
%       (refcha: channel to save as nrrd, for wholebrain)
%           (1 = red channel, default)
%       (refcha_seg: channel to save as nrrd, for seg)
%           (1 = red channel, default)
%           (2 = green channel)
%       (nchannels: number of channels)
%           (2, default)
%       %%% smoothing settings %%%
%       (sig: sigma)
%           (2, default)
%       (size: size)
%           (3, default)
%       %%% resampling settings %%%
%       (oXYZres: output spatial resolution (width, height, depth))
%           ([1.2 1.2 1] um, default)
%       (fres: number to use for rounding spatial resolution)
%           (10^4, default)
%       %%% padding settings %%%
%       (padgate: gate to pad)
%           (0, default)
%       (padnum: number of planes, lines to use for padding)
%           (10, default)
%       (fshift: shift distribution of F)
%           (0, default)
%       (zflipgate: gate to flip orientation in Z)
%           (1, default)
%       (segnplanes: original number of planes per segment)
%           (9, default)
%       (iDir: files origin)
%       (zcorrect: correct slicing in Z)
%           (1, default)
%       (flyname: '0' each file has its own brain image, '1', many files have the same fly brain)
%           (0, default)
%
% Notes:
% 1) For segments it reads the red channel from wDat and saves it as nrrd image
%   (uses pixel resolution provided in wDat)
% 2) For whole brains it gets the resolution from the nrrd image
% 3) Performs some processing to whole brains:
%   *reshape
%   *choose reference channel
%   *flip Z axes
%   *blur
%   *resample
%   *mirror
%   *output image is saved as nrrd
%
% To do:
% fix bug related to imprepros_correctpath(wDat.cDir)

% Default iparams
z2spars.cDir = pwd; 
z2spars.fisuffix = 'prosmetadata.mat';
z2spars.fo2reject = {'.', '..', 'local-T'};
z2spars.redo = [0 0];
z2spars.refcha = 1;
z2spars.refcha_seg = 1;
z2spars.refcha_suffix = '_01';
z2spars.nchannels = 2;
z2spars.sig = 2;
z2spars.size = 3;
z2spars.oXYZres = [1.2 1.2 1];
z2spars.fres = 10^4;
z2spars.padgate = 0; 
z2spars.padnum = 10;
z2spars.fshift = 0;
z2spars.zflipgate = 1;
z2spars.segnplanes = 9;
z2spars.iDir = [];
z2spars.zcorrect = 1;
z2spars.flyname = 0;

% update variables
if ~exist('FolderName', 'var') || isempty(FolderName); FolderName = []; end
if ~exist('FileName', 'var') || isempty(FileName); FileName = []; end
if ~exist('iparams', 'var'); iparams = []; end
z2spars = loparam_updater(z2spars, iparams);

% Selecting folders
fo2run = dir; 
fo2run = str2match(FolderName, fo2run);
fo2run = str2rm(z2spars.fo2reject, fo2run); 
fo2run = {fo2run.name};

fprintf(['Running n-folders : ', num2str(numel(fo2run)), '\n'])

for i = 1:numel(fo2run)
    
    fprintf(['Running folder : ', fo2run{i}, '\n']); 
    cd(fo2run{i});
    gennrrd(FileName, z2spars);
    cd(z2spars.cDir)
    
end

fprintf('... Done\n')

end

function gennrrd(FileName, z2spars)
% gennrrd: run per folder
%
% Usage:
%   gennrrd(FileName, z2spars)
%
% Args:
%   FileName: files to run
%   z2spars: input parameters

fly2run = rdir(['.', filesep, '*', z2spars.fisuffix]); 
fly2run = str2match(FileName, fly2run);
fly2run = {fly2run.name}; 
fly2run = strrep(fly2run, ['.', filesep], '');
fly2run = strrep(fly2run, z2spars.fisuffix, '');

for fi = 1:numel(fly2run)
    
    fprintf(['Running fly :', strrep(fly2run{fi}, '_', ' '), '\n'])
    genrrdIm(fly2run{fi}, z2spars)
    
end

end

function genrrdIm(flyname, z2spars)
% genrrdIm: per file name, generate _s and _w nrrd images
%
% Usage:
%   genrrdIm(flyname, z2spars)
%
% Args:
%   flyname: file to run
%   z2spars: input parameters

% 1) Load segment and generate the nrrd image
% 1.1) load segment res and side of segment associated to this fly
flyname = flyname(1:end-1);
load([flyname, '_', strrep(z2spars.fisuffix, '.mat', ''), '.mat'], 'wDat')
    
% 1.2) Make nrrd of the segment which is assumed to be already flipped if necessary
if ~exist([flyname, '_s', z2spars.refcha_suffix, '.nrrd'], 'file') || z2spars.redo(1) == 1
    
    if z2spars.zcorrect
        % In this case it passes the resolution described in "wDat.XYZres"
        if z2spars.refcha_seg == 1
            [~, cIm, oRes] = wDat_Z_correctslicing(...
                wDat, z2spars.segnplanes, wDat.RedChaMean);
        elseif z2spars.refcha_seg == 2
            [~, cIm, oRes] = wDat_Z_correctslicing(...
                wDat, z2spars.segnplanes, wDat.GreenChaMean);
        end
    else
        oRes = wDat.XYZres{2}(1:3);
        if z2spars.refcha_seg == 1
            cIm  = wDat.RedChaMean;
        elseif z2spars.refcha_seg == 2
            cIm  = wDat.GreenChaMean;
        end
    end
    
    cIm  = cIm + z2spars.fshift;
    fprintf(['Image resolution ', ...
            num2str(oRes(1)), ' ', num2str(oRes(2)), ' ', num2str(oRes(3)), '\n']);
    nrrdWriter([flyname, '_s', z2spars.refcha_suffix, '.nrrd'], ...
        mat2uint16(cIm, 0), oRes, [0 0 0], 'raw');
    clear iXYZres    
    
else
    
    fprintf('File _s_ already generated\n')
    
end

% 2) Load whole brain and generate the nrrd image
% 2.1) load whole brain res of this fly
% provide correct path to data directory

if ~isempty(z2spars.iDir)
    origCdir = z2spars.iDir;
else
    origCdir = imprepros_correctpath(wDat.cDir);
end

if z2spars.flyname
    
    % remove last underline
    TempS = strsplit2(flyname, '_');
    flyname = [TempS{1}, '_', TempS{2}];
    
end

if ~exist([flyname, ['_w', z2spars.refcha_suffix], '.nrrd'], 'file') || z2spars.redo(2) == 1
    
    % 2.1) reshape Data (necessary when it has 2 channels, all the cases)
    fprintf('reshaping data, ')
    [Data, meta] = nrrdread(fullfile([origCdir, filesep, flyname, '_Zstack','.nrrd']));
    iXYZres = round(nrrdread_res(meta)*z2spars.fres)/z2spars.fres;
    siz = size(Data);
    Data = double(reshape(Data, ...
        [siz(1:2), z2spars.nchannels, prod(siz(3:end))/(z2spars.nchannels)]));
    Data = permute(Data, [1 2 4 3]);
    
    % reduce number of channels to ref channel
    Data = Data(:, :, :, z2spars.refcha);

    % 2.2) Correct for z order
    if z2spars.zflipgate
        Data = flip(Data, 3);
    end
    
    % 2.3) smooth just up to the 2nd dimension (in this case this cleans up the resampling)
    if ~isempty(z2spars.sig)
        % hardcoded-pre smoothing:
        fprintf('smoothing data, ');
        Data = imblur(Data, z2spars.sig, 3, 2);
    end

    % 2.4) mirror if side is on the right (information comes from segment)
    fprintf('mirroring image, ');
    if ~isempty(strfind(wDat.bSide, 'R'))
        Data = flip(Data, 2);
    end

    % 2.5) resample (z2spars.oXYZres: (width, height, depth))
    if sum(z2spars.oXYZres == iXYZres) ~= 3
        
        fprintf(['resample image, old(', ...
            num2str(iXYZres(1)), ' ', num2str(iXYZres(2)), ...
            ' ', num2str(iXYZres(3)), ...
            ') new(', num2str(z2spars.oXYZres(1)), ...
            ' ', num2str(z2spars.oXYZres(2)), ' ', ...
            num2str(z2spars.oXYZres(3)), ') ']);
        Data = interp3DxT(Data, iXYZres, z2spars.oXYZres, 3);
        %interp3DxT(Data3DxT, XYZinit, XYZfinal, LastDim)
        
    end
    
    % 2.6) smooth just up to the 2nd dimension (in this case this blurs the image)
    if ~isempty(z2spars.sig) && ~isempty(z2spars.size)
        
        fprintf(', smoothing data, ')
        Data = imblur(Data, z2spars.sig, z2spars.size, ...
            max([length(z2spars.size), length(z2spars.sig)]));
        
    end
    
    % 2.7) pad image in Z (adds zero frames above and below)
    if z2spars.padgate
        
        fprintf('padding image, ')
        if length(size(Data)) == 3
            Data = padarray(Data, [0, 0, z2spars.padnum]); 
        else
            Data = padarray(Data, [0, 0, z2spars.padnum, 0]); 
        end
        
    end

    % 2.8) make and save nrrd of the segment
    fprintf('saving ... ')
    nrrdWriter([flyname, '_w', z2spars.refcha_suffix, '.nrrd'], ...
        mat2uint16(Data, 0), z2spars.oXYZres, ...
        [0 0 0], 'raw');
    fprintf('Done\n\n')
    
else
    
    fprintf('File _w_ already generated\n')
    
end

end
