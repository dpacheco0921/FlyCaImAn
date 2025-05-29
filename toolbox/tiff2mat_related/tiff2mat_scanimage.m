function [Im, ImMeta] = tiff2mat_scanimage(tifname, datatype, verbose, timeavg_flag)
% tiff2mat_scanimage: reads metadata and data from a tiff file
%
% Usage:
%   [Im, ImMeta] = tiff2mat_scanimage(tifname, datatype, verbose, timeavg_flag)
%
% Args:
%   tifname: name of tiff file to load
%   datatype: type of data
%   verbose: verbose
%   timeavg_flag: flag to load and average signal per plane directly from tiff instead of loading full matrix
%           (default: 0)
%
% Returns:
%   Im: image from tiff files
%   ImMeta: metadata
%       (Y: Lines per frame or height)
%       (X: Pixels per line or Width)
%       (Z: planes per volume)
%       (ChNum: number of channels collected)
%       (Zoom: zoom factor)
%       (Power: power)
%       (framerate: frame rate)
%       (volumerate: volume rate)
%       (sympixels: flag for symetric pixels)
%       (Imclass: class of encoding)
%       (Z_stepsiz: step size in Z)

if ~exist('verbose', 'var') || isempty(verbose)
    verbose = 0;
end

if ~exist('verbose', 'var') || isempty(verbose)
    verbose = 0;
end

if ~exist('timeavg_flag', 'var') || isempty(timeavg_flag)
    timeavg_flag = 0;
end

try
    
    info = imfinfo(fullfile([tifname, '.tif']));
    Imclass = getClass(info);
    
    if contains(datatype, {'2DxT_single', 'old'})
        [Y, X, Channels, Zoom, Power, Z, ...
            framerate, volumerate, sympixels] = ...
            TiffMetadataOld(info);
    else
        [Y, X, Channels, Zoom, Power, Z, ...
            framerate, volumerate, sympixels, Z_stepsiz] = ...
            TiffMetadata(info);
    end

    % getting frame offset
    StartOffset = cell2mat({info.StripOffsets}');
    StartOffset = StartOffset(:, 1);
    Frames = numel(StartOffset)/Channels;
    RepeatNum = Frames/Z;

    if Channels > 1
        StartOffset_per_cha = {StartOffset(1:2:end), StartOffset(2:2:end)};
    else
        StartOffset_per_cha = {StartOffset};
    end

    if timeavg_flag

        % load data per slice and average instead of loading full
        %   matrix
        
        StartOffset_per_cha{1} = reshape(StartOffset_per_cha{1}, [Z, RepeatNum]);
        if Channels > 1
            StartOffset_per_cha{2} = reshape(StartOffset_per_cha{2}, [Z, RepeatNum]);
        end
    
        % prealocating data
        Im = single(zeros(Y, X, Z, Channels));
        
        % for each plane get the average image
        for Z_idx = 1:Z
    
            fprintf(['Running plane: ', num2str(Z_idx), '\n'])

            for r_i = 1:RepeatNum
    
                % odd, Green
                temp_im_ch1(:, :, r_i, 1) = ...
                    frame2mat(tifname, Y, X, ...
                    StartOffset_per_cha{1}(Z_idx, r_i), Imclass);
    
                if Channels > 1
                    
                    % even, Red        
                    temp_im_ch2(:, :, r_i, 1) = ...
                        frame2mat(tifname, Y, X, ...
                        StartOffset_per_cha{2}(Z_idx, r_i), Imclass);        
                            
                end            
    
            end
    
            % generate the average
            Im(:, :, Z_idx, 1) = mean(temp_im_ch1, 3);
            if Channels > 1
                Im(:, :, Z_idx, 2) = mean(temp_im_ch2, 3);            
            end
    
        end

    else
             
        % load full matrix and do processing after
        
        % prealocating data
        Im = single(zeros(Y, X, Frames, Channels));
    
        for FrameIdx = 1:Frames
            
            % odd, Green
            Im(:, :, FrameIdx, 1) = ...
                frame2mat(tifname, Y, X, ...
            StartOffset_per_cha{1}(FrameIdx), Imclass);             
            
            if Channels > 1

                % even, Red
                Im(:, :, FrameIdx, 1) = ...
                    frame2mat(tifname, Y, X, ...
                    StartOffset_per_cha{2}(FrameIdx), Imclass);
                
            end

        end
    end

    ImMeta.X = X;
    ImMeta.Y = Y;
    ImMeta.Z = Z;
    ImMeta.ChNum = Channels;
    ImMeta.Zoom = Zoom;
    ImMeta.Power = Power;
    ImMeta.framerate = framerate;
    ImMeta.volumerate = volumerate;
    ImMeta.sympixels = sympixels;    
    ImMeta.Imclass = Imclass;
    
    if exist('Z_stepsiz', 'var') && ...
            ~isempty(Z_stepsiz)
        ImMeta.Z_stepsiz = Z_stepsiz;
    end
    
catch error
    
    if verbose == 1
        fprintf(['could not run file ', strrep(tifname, '\', ''), '\n']);
        keyboard
        display(error);
    end
    
end

end

function Im = frame2mat(tifname, Y, X, StartOffset, Imclass)
% frame2mat: reads frames from tiff file
%
% Usage:
%   Im = frame2mat(BaseName, Y, X, StartOffset, Imclass)
%
% Args:
%   tifname: name of tiff file to load
%   Y: dim
%   X: dim
%   StartOffset: scanimage variable
%   Imclass: type of data encoding
%
% Returns:
%   Im: image from tiff files
%
% Notes:
% this might change when scanimage is fixed
% "predata = double(reshape(predata, X, Y)); "
% however it does not matter if X and Y are the same value

fid = fopen(fullfile([tifname, '.tif']));
fseek(fid, StartOffset, 'bof');

% We knew this from BitsPerSample
Im = fread(fid, Y*X, ['*', Imclass]);

fclose(fid);
Im = single(reshape(Im, X, Y)');

end

function Imclass = getClass(info)
% getClass: get class of encoding
%
% Usage:
%   Imclass = getClass(info)
%
% Args:
%   info: tiff metadata
% 
% Returns:
%   Imclass: tiff data class

microscopetxt = info(1).ImageDescription;

if ~contains(microscopetxt, 'DataType = ')
    microscopetxt = info(1).Software;
end

if contains(microscopetxt, 'int16')
    Imclass = 'int16';
elseif contains(microscopetxt, 'uint16')
    Imclass = 'uint16';
elseif contains(microscopetxt, 'int8')
    Imclass = 'int8';
elseif contains(microscopetxt, 'uint8')
    Imclass = 'uint8';
end

end

function [Y, X, Channels, Zoom, Power, Z, ...
    framerate, volumerate, sympixels] = TiffMetadataOld(info)
% TiffMetadataOld: get metadata from tiff files
%
% Usage:
%   [Y, X, Channels, Zoom, Power, Z] = TiffMetadataOld(info)
%
% Args:
%   info: tiff metadata (using imfinfo)
% 
% Returns:
%   Y: Lines per frame or height
%   X: Pixels per line or Width
%   Channels: number of channels collected
%   Zoom: zoom factor
%   Power: power
%   Z: planes per volume
%   framerate: frame rate
%   volumerate: volume rate
%   sympixels: flag for symetric pixels

params = strfind(info(1,1).ImageDescription, 'state.');

for pIdx = 2:numel(params)
    if pIdx == numel(params)
        eval([info(1,1).ImageDescription(params(pIdx):end-1),';'])
    else
        eval([info(1,1).ImageDescription(params(pIdx):params(pIdx+1)-2),';'])
    end
end

Y = state.acq.linesPerFrame;
X = state.acq.pixelsPerLine;
Channels = state.acq.numberOfChannelsAcquire;
Zoom = state.acq.zoomFactor;
Power = [];
Z = [];
framerate = [];
volumerate = [];
sympixels = 0;

end

function [Y, X, Channels, Zoom, Power, Z, ...
    framerate, volumerate, sympixels, Z_stepsiz] = TiffMetadata(info)
% TiffMetadata: get metadata from tiff files
%
% Usage:
%   [Y, X, Channels, Zoom, Power, Z, ...
%       framerate, volumerate, sympixels, Z_stepsiz] = TiffMetadata(info)
%
% Args:
%   info: tiff metadata (using imfinfo)
% 
% Returns:
%   Y: Lines per frame or height
%   X: Pixels per line or Width
%   Channels: number of channels collected
%   Zoom: zoom factor
%   Power: power
%   Z: planes per volume
%   framerate: frame rate
%   volumerate: volume rate
%   sympixels: flag for symetric pixels
%   Z_stepsiz: step size in Z
%
% Notes:
% So far compatible with old and new scanimage

params = strfind(info(1,1).ImageDescription, 'scanimage.');

if ~isempty(params)
    % works for scanimage 2014
    
    for pIdx = 2:numel(params)
        try
            if pIdx == numel(params) || pIdx == 103 ...
                    || pIdx == 183 || pIdx == 184 || pIdx == 185
                %eval([info(1,1).ImageDescription(params(pIdx):end-1),';'])
            else
                eval([info(1,1).ImageDescription(params(pIdx):params(pIdx+1)-2), ';'])
            end
        catch
            fprintf([num2str(pIdx), ' '])
        end
    end
    
else
    
    % works for scanimage 2016 and later versions
    params = strfind(info(1,1).Software, 'SI.');

    for pIdx = 1:numel(params)
        try
            if pIdx == 223
            else
                eval([info(1,1).Software(params(pIdx):params(pIdx+1)-2), ';'])
            end
        catch
            fprintf([num2str(pIdx), ' '])
        end
    end
    
    scanimage.SI = SI;
    
end

% parameters to collect
Y = scanimage.SI.hRoiManager.linesPerFrame;
X = scanimage.SI.hRoiManager.pixelsPerLine;

% determine which mode of zstack was used
stackmode = [];
Z = 1;
try 
    stackmode = scanimage.SI.hStackManager.stackMode;
end

if ~isempty(stackmode)
    if contains(stackmode, 'slow')
        Z = scanimage.SI.hStackManager.actualNumSlices;
    elseif contains(stackmode, 'fast')
        Z = scanimage.SI.hStackManager.numFramesPerVolumeWithFlyback;
    end
else
    try
        Z = scanimage.SI.hFastZ.numFramesPerVolume;
    catch
        if isfield(scanimage.SI, 'hStackManager') && ...
            scanimage.SI.hStackManager.enable
            Z = scanimage.SI.hStackManager.numSlices;
        end
    end
end

Channels = sum(scanimage.SI.hChannels.channelSave > 0);
Zoom = scanimage.SI.hRoiManager.scanZoomFactor;
Power = scanimage.SI.hBeams.powers;
framerate = scanimage.SI.hRoiManager.scanFrameRate;
volumerate = scanimage.SI.hRoiManager.scanVolumeRate;
sympixels = scanimage.SI.hRoiManager.forceSquarePixels;

Z_stepsiz = [];
if isfield(scanimage.SI, 'hStackManager')
    Z_stepsiz = scanimage.SI.hStackManager.stackZStepSize;
end

end
