function batch_zstack2segmatch_t(FolderName, FileName, iparams)
% batch_zstack2segmatch_t: generates zstacks used for segment to fly and 
%   fly to atlas registration. Whole fly image has the pattern (_w) and
%   fly segment image has the pattern (_s). It also does some preprocessing
%   of images (see notes).
%
% Usage:
%   batch_zstack2segmatch_t(FolderName, FileName, iparams)
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
%       (zflipgate: gate to flip orientation in Z)
%           (1, default)
%       (mirror_flag: flag to also save mirror image (_mw))
%           (1, default)
%       (oDir: output directory)
%           (pwd, default)
%
% Notes:
% 1) For segments it reads the red channel from wDat and saves it as nrrd image
%   (uses pixel resolution provided in wDat)
% 2) For whole brains it gets the resolution from the metadata.mat 
%   file and performs some corrections:
%   *reshape, *choose ref channel, *flip Z axes, *blur, *resample, *mirror and *save as nrrd

% Default iparams
z2spars = [];
z2spars.cDir = pwd;
z2spars.fisuffix = '_Zstack.nrrd';
z2spars.fo2reject = {'.', '..', 'BData', 'preprocessed'};
z2spars.redo = 0;
z2spars.refcha = [1, 2];
z2spars.nchannels = 2;
z2spars.sig = 2; 
z2spars.size = 3;
z2spars.oXYZres = [1.2 1.2 1];
z2spars.fres = 10^4;
z2spars.padgate = 0;
z2spars.padnum = 10;
z2spars.zflipgate = 1; 
z2spars.mirror_flag = 1;
z2spars.subfolder = [];
z2spars.oDir = pwd;
%[~, ~, z2spars.oDir] = ...
%    imprepros_pathdiradd(z2spars.cDir);

% update variables
if ~exist('FolderName', 'var') || isempty(FolderName)
    FolderName = [];
end
if ~exist('FileName', 'var') || isempty(FileName)
    FileName = [];
end
if ~exist('iparams', 'var'); iparams = []; end
z2spars = loparam_updater(z2spars, iparams);

% Selecting folders
fo2run = dir;
fo2run = str2match(FolderName, fo2run);
fo2run = str2rm(z2spars.fo2reject, fo2run);
fo2run = {fo2run.name};
fprintf(['Running n-folders : ', ...
    num2str(numel(fo2run)), '\n'])

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
    
    fprintf(['Running fly :', ...
        strrep(fly2run{fi}, '_', ' '), '\n']);
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

% generate target folder
if ~isempty(z2spars.oDir)
    targetDir = [z2spars.oDir, filesep, ...
        z2spars.subfolder, filesep];
else
    targetDir = ['.', filesep];
end

if ~exist(targetDir, 'dir')
    mkdir(targetDir)
end

if ~exist([targetDir, flyname, '_w_01.nrrd'], 'file') ...
        || ~exist([targetDir, flyname, '_w_02.nrrd'], 'file') ...
        || ~exist([targetDir, flyname, '_wm_01.nrrd'], 'file') ...
        || ~exist([targetDir, flyname, '_wm_02.nrrd'], 'file') ...
        || z2spars.redo

    % 1.1) reshape Data (necessary when it has 2 channels, all the cases)
    fprintf('reshaping data, ')
    [Data, meta] = nrrdread(fullfile([flyname, z2spars.fisuffix]));
    iXYZres = nrrdread_res(meta);
    siz = size(Data);
    Data = double(reshape(Data, ...
        [siz(1:2), z2spars.nchannels, ...
        prod(siz(3:end))/(z2spars.nchannels)]));

    Data = permute(Data, [1 2 4 3]);
    
    % reduce number of channels to ref channel
    Data = Data(:, :, :, z2spars.refcha);
   
    % 1.2) Correct for z order (flip in z axis)
    if z2spars.zflipgate
        Data = flip(Data, 3);
    end
    
    % 1.3) smooth just up to the 2nd dimension (in this case this cleans up the resampling)
    if ~isempty(z2spars.sig) % hardcoded-pre smoothing
        fprintf('smoothing data, ');
        Data = imblur(Data, z2spars.sig(1), 3, 2);
    end
    
    if isempty(z2spars.oXYZres)
        z2spars.oXYZres = iXYZres;
    end
    
    % 1.4) resample (z2spars.oXYZres: (width, height, depth))
    if sum(z2spars.oXYZres == iXYZres) ~= 3
        
        fprintf(['resample image, old(', ...
            num2str(iXYZres(1)), ' ', ...
            num2str(iXYZres(2)), ' ', num2str(iXYZres(3)), ...
            ') new(', num2str(z2spars.oXYZres(1)), ...
            ' ', num2str(z2spars.oXYZres(2)), ' ', ...
            num2str(z2spars.oXYZres(3)), ')\n']);
        Data = interp3DxT(Data, iXYZres, z2spars.oXYZres, 3);
        %interp3DxT(Data3DxT, XYZinit, XYZfinal, LastDim)

    end

    % 1.5) smooth just up to the 2nd dimension (in this case this blurs the image)
    if ~isempty(z2spars.sig) && ~isempty(z2spars.size)
        
        fprintf(', smoothing data, ')
        Data = imblur(Data, z2spars.sig, z2spars.size, ...
            max([length(z2spars.size), length(z2spars.sig)]));
        
    end
   
    % 1.6) pad image in Z (adds zero frames above and below)
    if z2spars.padgate
        
        fprintf('padding image, ')
        if length(size(Data)) == 3
            Data = padarray(Data, [0, 0, z2spars.padnum]); 
        else
            Data = padarray(Data, [0, 0, z2spars.padnum, 0]); 
        end
        
    end

    % 1.7) make and save nrrd of the segment
    % assumes 1st channel is template and saves it accordingly
    for t = 1:size(Data, 4)
        
        tData = Data(:, :, :, t);
        
        % original
        lname = [flyname, '_w_0', num2str(t), '.nrrd'];
        nrrdWriter(lname, mat2uint16(tData, 0), ...
            z2spars.oXYZres, [0 0 0], 'gzip');
        
        if ~strcmp(targetDir, ['.', filesep])
            movefile(['.', filesep, lname], ...
                [targetDir, lname])
        end
        
        % mirror
        if z2spars.mirror_flag
            tData = flip(tData, 2);
            lname = [flyname, '_wm_0', num2str(t), '.nrrd'];
            nrrdWriter(lname, mat2uint16(tData, 0), ...
                z2spars.oXYZres, [0 0 0], 'gzip');
            clear tData
            
            if ~strcmp(targetDir, ['.', filesep])
                movefile(['.', filesep, lname], ...
                    [targetDir, lname])
            end
            
        end
        
    end
    
else
    
    fprintf('File already generated\n\n')
    
end

end
