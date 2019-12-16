function [Im] = split_im_into_cc(Im, siz_ths, bin_flag)
% split_im_into_cc: splits image into connected components
%   and prunes a binary image based on size of connected
%   components
%
% Usage:
%   [Im] = split_im_into_cc(Im, siz_ths, bin_flag)
%
% Args:
%   Im: input image
%   siz_ths: size threshold to prune connected components
%   bin_flag: homogenize connected component index (to a single label)

if ~exist('siz_ths', 'var') || ...
        isempty(siz_ths)
    siz_ths = 0;
end

if ~exist('bin_flag', 'var') || ...
        isempty(bin_flag)
    bin_flag = 1;
end

l_cs = bwconncomp(Im);
pix2del = find(cellfun(@numel, l_cs.PixelIdxList) >= siz_ths);
Im = zeros(size(Im));

for i = 1:numel(pix2del)
    Im(l_cs.PixelIdxList{pix2del(i)}) = i;
end

% convert to a single component
if bin_flag
    Im = logical(Im);
end

fprintf('Done\n')

end
