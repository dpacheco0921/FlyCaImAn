function plot_roi_coverage(filename, ...
    vid2plot, wDat, roi, oDir, roi2sel, ...
    cha2use, com_flag, indroi_color_flag)
% plot_roi_coverage: plot roi coverage
%   it overlays both sum of max norm weights 
%   and binarized voxels to defined channel wDat.GreenChaMean/wDat.RedChaMean
%
% Usage:
%   plot_roi_coverage(filename, ...
%       vid2plot, wDat, roi, oDir, roi2sel, ...
%       cha2use, com_flag, indroi_color_flag)
%
% Args:
%   filename: file name
%   vid2plot: type of video to plot
%       [1 0 0] energy normalized, [0 1 0] max normalized, [0 0 1] binarized
%       (default, [1 1 1])
%   wDat: metadata parameter
%   roi: roi parameter
%   oDir: output directory
%   roi2sel: indeces of ROIs to use
%   channel2use: define the channel to use
%   com_flag: flag to overlay ROI number using its center of mass
%   indroi_color_flag: flag to use different color for each ROI

if ~exist('vid2plot', 'var') || isempty(vid2plot)
   vid2plot = [1 1 1];
end

if ~exist('roi2sel', 'var') || isempty(roi2sel)
   roi2sel = [];
end

if ~exist('oDir', 'var') || isempty(oDir)
   oDir = pwd;
end

if ~exist('cha2use', 'var') || isempty(cha2use)
   cha2use = 'wDat.GreenChaMean';
end

if ~exist('com_flag', 'var') || isempty(com_flag)
   com_flag = 0;
end

if ~exist('indroi_color_flag', 'var') || isempty(indroi_color_flag)
   indroi_color_flag = 0;
end

if ~exist('com_text_cmap', 'var')
    com_text_cmap = 'cyan';
end

% deafult video settings
vidpar.xyzres = wDat.XYZres{2};
vidpar.overcor = [1 0 0];
vidpar.lag = 0.01;
vidpar.vgate = 1;
        
% set brackground image and intensity range
eval(['vidpar.Y2 = double(', cha2use,' )']); 
vidpar.Y2range = [0 prctile(vidpar.Y2(:), 95)];
vidpar.sizY = size(vidpar.Y2);

if com_flag
    
    % add ROI center of mass info and text
    vidpar.Y3text = cf(@(x) num2str(x), ...
        chunk2cell(1:size(roi.A, 2), 1));
    vidpar.Y3textcmap = com_text_cmap;
    vidpar.Y3 = round(com(roi.A, wDat.vSize(1), ...
        wDat.vSize(2), wDat.vSize(3)));
    
end

% select ROIs
if ~isempty(roi2sel)
   roi.A = roi.A(:, roi2sel); 
end

% 1) energy normalized
A2use = [roi.A, roi.b];
K = size(A2use, 2);
nA = sqrt(sum(A2use.^2))';
A2use = A2use/spdiags(nA, 0, K, K);
A2use = A2use(:, 1:end-1);

% collapse all ROIs into a single one
if indroi_color_flag
    
    % plot each ROI independently (this is very memory demanding only do for tens of ROIs
    %   or cluster them prior)
    A2use = double(full(A2use));
    vidpar.overcor = jet(size(A2use, 2));
    
else
    
    A2use = double(full(sum(A2use, 2)));
    
end

if vid2plot(1)
    
    % plot video
    temp_A = A2use(:);
    vidpar.range = [min(temp_A(:))...
        prctile(temp_A(temp_A ~= 0), 95)];
    clear temp_A
    
    vidpar.vname = [oDir, filesep, filename, '_roicov_energyn'];
    slice3Dmatrix(A2use, vidpar);
    
end

% 2) max normalized
A2use = roi.A;
A2use = bsxfun(@times, A2use, 1./max(A2use, [], 1));

% collapse all ROIs into a single one
if indroi_color_flag
    
    % plot each ROI independently (this is very memory demanding only do for tens of ROIs
    %   or cluster them prior)
    A2use = double(full(A2use));
    vidpar.overcor = jet(size(A2use, 2));
    
else
    
    A2use = double(full(max(A2use, [], 2))); 
    
end

if vid2plot(2)
    
    % plot video
    vidpar.range = [min(A2use(:))...
        prctile(A2use(A2use ~= 0), 95)];
    vidpar.vname = [oDir, filesep, filename, '_roicov_maxn'];
    slice3Dmatrix(A2use, vidpar);
    
end

% 3) binarized
A2use = double(A2use > 0); 

if vid2plot(3)
    
    % plot video
    vidpar.range = [0 2];
    vidpar.vname = [oDir, filesep, filename, '_roicov_bin'];
    slice3Dmatrix(A2use, vidpar);
    
end

end
