function plot_fused_image_vid(...
    float_im_1, float_im_2, filename, ...
    oDir, iparams)
% plot_fused_image_vid: plot a video of fused input images
%
% Usage:
%   plot_fused_image_vid(...
%       float_im_1, float_im_2, filename, oDir)
%
% Args:
%   float_im_1: floating image 1
%   float_im_2: floating image 2
%   filename: name of output video
%   oDir: output folder to save figures
%   iparams: parameters to update
%       (vquality: video quality)
%           (default, 100)
%       (frate: video rate, frames per second)
%           (default, 10)
%
% Notes:
% Writes a video file named oDir/filename (no return value). Each video
%   frame is one plane (3rd dimension) with image 1 in red and image 2 in green.

% default video settings (user values in iparams take precedence)
ip.vquality = 100;
ip.frate = 10;
if ~exist('iparams', 'var'); iparams = []; end
ip = loparam_updater(ip, iparams);

% Normalize images, scale to [0 1]
min_ = prctile([float_im_1(:); float_im_2(:)], 1);
max_ = prctile([float_im_1(:); float_im_2(:)], 99);

float_im_1 = float_im_1 - min_;
float_im_1 = float_im_1/max_;

float_im_2 = float_im_2 - min_;
float_im_2 = float_im_2/max_;

figH = figure('Position', ...
    genfigpos(1, 'center', [500 500]));
axH = subplot(1, 1, 1);

% video related
vHandle = VideoWriter([oDir, filesep, filename]);
vHandle.Quality = ip.vquality;
vHandle.FrameRate = ip.frate;
open(vHandle);

for i = 1:size(float_im_1, 3)
   
        im_1 = mat2gray(float_im_1(:, :, i), [0 1]);
        im_2 = mat2gray(float_im_2(:, :, i), [0 1]);

        C = imfuse(im_1, im_2, ...
           'falsecolor', 'Scaling', ...
           'joint', 'ColorChannels', ...
           [1 2 0]);
        imshow(C, 'Parent', axH)
        
        % append frames to video
        writeVideo(vHandle, getframe(axH));
        
end

% close video object
close(vHandle); 
close(gcf);

end
