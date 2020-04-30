function [Y, imtime, vidtime, SVD_Dat] = ...
    load_edit_video_SVD(wDat, filename, sig, siz, tres_flag)
% load_edit_video_SVD: load SVD decomposition of video motion energy
%
% Usage:
%   [Y, imtime, vidtime, SVD_Dat] = ...
%       load_edit_video_SVD(wDat, filename, sig, siz, tres_flag)
%
% Args:
%   wDat: metadata variable
%   filename: filename
%   sig: std of gaussian kernel
%   siz: size of kernel
%   tres_flag: flag to resample to imaging time resolution

if ~exist('sig', 'var')
    sig = [];
end

if ~exist('siz', 'var')
    siz = [];
end

if ~exist('tres_flag', 'var') || isempty(tres_flag)
    tres_flag = 0;
end

if exist(filename, 'file')
    load(filename, 'proc')
else
    fprintf('*_proc.mat fiel does not exist - do facemap\n')
    return
end

% get imaging timepoints
imtime = wDat.fTime;   
imtime = round(imtime*100)/100;
vidtime = wDat.vid.fstEn{1}(:, 2);

% get temporal components
Y = proc.motSVD{1};

% reduce timepoints
Y = Y(1:numel(vidtime), :);

% smooth traces
if ~isempty(sig) && ~isempty(siz)
   Y = imblur(Y, sig, siz, 1); 
end

Y = Y';

% pass whole variable
SVD_Dat = proc;

if tres_flag
    
    % resample to imaging resolution defined in wDat.fTime
    Y = interp1(vidtime, Y', imtime, 'linear');
    Y = framegapfill(find(isnan(Y(:, 1))), Y');
    
end

end
