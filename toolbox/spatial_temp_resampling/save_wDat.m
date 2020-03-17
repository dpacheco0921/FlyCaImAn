function save_wDat(filename, datatype, ...
    imdirection, bkgate, fshift, blowcap)
% save_wDat: load all metadata variables and
%   pass relevant fields to wDat and save 
%   wDat (for compatibility with postprocessing steps like ROIseg)
%
% Usage:
%   save_wDat(ffilename, datatype, ...
%       imdirection, bkgate, fshift, blowcap)
%
% Args:
%   filename: file name
%   datatype: data type (2DxT, or 3DxT)
%   imdirection: Z axis direction
%   bkgate: gate to do custom background substraction
%   fshift: value to shift distribution
%   blowcap: minimun pixel value (it replaces values less than this)
%
% Notes:
% It assumes that either Greenchannel or Redchannel exist and is not empty

mcDat = [];
cDat = [];
wDat = [];

% load all metadata variables
load(filename, 'iDat', 'fDat', 'mcDat', 'lStim');
load(filename, 'vidDat');

if ~exist('vidDat', 'var')
    vidDat = [];
end

if contains(fDat.DataType, 'opto') && ...
        ~contains(fDat.DataType, 'prv')
    load(filename, 'cDat');
end

% stitch related
wDat.Zstitch.Zshift = [];
wDat.Zstitch.Zend = [];
wDat.Zstitch.Zidx = []; 
wDat.Zstitch.Xshift = [];
wDat.Zstitch.Yshift = [];

% motion correction
wDat.MotCor.Zshift = [];
wDat.MotCor.sXshift = [];
wDat.MotCor.sYshift = [];

% video related
wDat.vid.varNames = [];
wDat.vid.var = [];
wDat.vid.fstEn = [];

% flip z-axis (3D data only)
if contains(datatype, '3DxT')
    
    if contains(imdirection, 'invert')
        
        if ~isfield(iDat, 'RedChaMean') || ...
                ~isfield(iDat, 'RedChaMean')
            iDat.RedChaMean = [];
            iDat.GreenChaMean = [];
        else
            iDat.RedChaMean = flip(iDat.RedChaMean, 3);
            iDat.GreenChaMean = flip(iDat.GreenChaMean, 3);
        end
        
    end
    
end

% generate mean-channels
try
    wDat.RedChaMean = iDat.RedChaMean;
catch
    wDat.RedChaMean = [];
end

try
    wDat.GreenChaMean = iDat.GreenChaMean;
catch
    wDat.GreenChaMean = [];
end

if isempty(wDat.GreenChaMean) ...
        && isempty(wDat.RedChaMean)
    
    fprintf('***********************************************************************\n')
    fprintf(['File probably failed at motion correction ', ...
        '(did not generate mean fields (RedChaMean GreenChaMean) in iDat)\n'])
    fprintf('***********************************************************************\n')
    
else
    
    % background correction
    if bkgate
        wDat.RedChaMean = wDat.RedChaMean - iDat.bs(1); 
        wDat.GreenChaMean = wDat.GreenChaMean - iDat.bs(2); 
    end

    if isfield(iDat, 'bs')
        wDat.bs = iDat.bs;
    else
        wDat.bs = [];
    end

    wDat.RedChaMean = wDat.RedChaMean + fshift(1); 
    wDat.RedChaMean(wDat.RedChaMean < blowcap) = blowcap;
    wDat.GreenChaMean = wDat.GreenChaMean + fshift(2);
    wDat.GreenChaMean(wDat.GreenChaMean < blowcap) = blowcap;

    if bkgate
        wDat.prepros.bsSubs = 1;
    end
    
end

% get stimuli info
wDat = getStimInfo(wDat, iDat, fDat, lStim, cDat, ...
    mcDat, vidDat, strrep(filename, '_metadata', ''), ...
    1, 1, [], iDat.FrameN);

if isempty(wDat.GreenChaMean) ...
        && isempty(wDat.RedChaMean)
    
    fprintf('***********************************************************************\n')
    fprintf(['File probably failed at motion correction ', ...
        '(did not generate mean fields (RedChaMean GreenChaMean) in iDat)\n'])
    fprintf('***********************************************************************\n')
    
else
    
    % find nan-slices
    if iDat.FrameN > 1

        try
            siz = size(wDat.GreenChaMean);
            nan_pix = isnan(wDat.GreenChaMean);
        catch
            siz = size(wDat.RedChaMean);
            nan_pix = isnan(wDat.RedChaMean);        
        end

        plane2keep = sum(reshape(nan_pix, ...
            [prod(siz([1 2])) siz(3)])) ~= prod(siz([1 2]));

    end

    % tag and remove nan-xy pixels
    if iDat.FrameN > 1
        if ~isempty(wDat.RedChaMean)
            wDat.mask = max(isnan(...
                wDat.RedChaMean(:, :, plane2keep)), [], 3);
        else
            wDat.mask = max(isnan(...
                wDat.GreenChaMean(:, :, plane2keep)), [], 3);
        end
    else
        if ~isempty(wDat.RedChaMean)
            wDat.mask = max(isnan(wDat.RedChaMean), [], 3);
        else
            wDat.mask = max(isnan(wDat.GreenChaMean), [], 3);
        end    
    end

    % bmask, use all pixels (2D data only)
    if iDat.FrameN == 1
        wDat.bMask = ones(size(wDat.GreenChaMean(:, :, 1)));
    end

    % update size
    if ~isempty(wDat.RedChaMean)
        wDat.fSize = [size(wDat.RedChaMean, 1), ...
            size(wDat.RedChaMean, 2)];
    else
        wDat.fSize = [size(wDat.GreenChaMean, 1), ...
            size(wDat.GreenChaMean, 2)];
    end

    if iDat.FrameN > 1

        if ~isempty(wDat.RedChaMean)
            wDat.vSize = [wDat.fSize, ...
                size(wDat.RedChaMean, 3)];
        else
            wDat.vSize = [wDat.fSize, ...
                size(wDat.GreenChaMean, 3)];
        end    
        wDat.vOrient = imdirection;

    else

        wDat.vSize = [wDat.fSize, 1];
        wDat.vOrient = [];

    end
end

save(filename, 'wDat', '-append');

end
