function brainmask = wDat_generatemask(wDat, minsize, ...
    image2replace, rmsmall_flag)
% wDat_generatemask: generates a logical brain mask based 
%   on wDat fields manually provided
%
% Usage:
%   brainmask = wDat_generatemask(wDat, minsize, ...
%       image2replace, rmsmall_flag)
%
% Args:
%   wDat: main metadata structure
%   minsize: min size of connected components to keep
%   image2replace: input image to use instead of wDat.GreenChaMean
%   rmsmall_flag: remove small connected components 
%
% Notes:
% wDat.bF_zbin defines the intervals of planes to be used for each
%   intensity thresholding

if ~exist('minsize', 'var') || isempty(minsize)
    minsize = 0;
end

if ~exist('image2replace', 'var') || isempty(image2replace)
    image2replace = [];
end

if ~exist('rmsmall_flag', 'var') || isempty(rmsmall_flag)
    rmsmall_flag = 1;
end

% replace GreenChannel
if ~isempty(image2replace)
    wDat.GreenChaMean = image2replace;
end

% generate mask (simple intensity thresholding)
for i = 1:numel(wDat.bF_zbin)-1
    fprintf('*')
    subv = wDat.GreenChaMean(:, :, ...
        wDat.bF_zbin(i)+1: wDat.bF_zbin(i+1));
    brainmask(:, :, ...
        wDat.bF_zbin(i)+1: wDat.bF_zbin(i+1)) = ...
        subv > wDat.bF_ths(1, i);
    clear subv hDat
end

% remove connected voxels that are small
if rmsmall_flag
    
    % Prune small components
    l_cs = bwconncomp(brainmask);
    pix2del = find(cellfun(@numel, l_cs.PixelIdxList) >= minsize);
   
    % reconstruct the mask
    brainmask = false(wDat.vSize);
    for i =  pix2del
        fprintf('.')
        brainmask(l_cs.PixelIdxList{i}) = 1;
    end
    
end

% replace edges by preceading planes
if size(brainmask, 3) > 1
    brainmask(:, :, 1) = brainmask(:, :, 2);
    if size(brainmask, 3) > 2
        brainmask(:, :, end) = brainmask(:, :, end-1);
    end
end

end
