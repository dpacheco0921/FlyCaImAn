function gen_fictrac_calibration_transform(...
    sphere_cx_px, sphere_cy_px, ...
    sphere_radius_px, width_pix, ...
    height_pix, vfov, oDir)
% gen_fictrac_calibration_transform: This function creates a new calibration_transform file
%   based on the position of the ball in the image (center, radius).
%
% Usage:
%	gen_fictrac_calibration_transform(...
%       sphere_cx_px, sphere_cy_px, ...
%       sphere_radius_px, width_pix, height_pix, vfov)
%
% Args:
%   sphere_cx_px: sphere points x
%   sphere_cy_px: sphere points y
%   sphere_radius_px: spehere radious
%   width_pix: width of mask in pixels
%   height_pix: height of mask in pixels
%   vfov: vertical field of view
%   oDir: output directory
% 
% Example
%   sphere_cx_px = 254.72;
%   sphere_cy_px = 347.14;
%   sphere_radius_px = 175.52;
%
% Notes
%   vfov, to calculate this see:
%       1) https://www.reddit.com/r/fictrac/comments/e71ida/how_to_get_the_right_vfov/
%       2) Use a ruler and:
%           2.1) calculate distance from camera to ruler
%               distance_ = 18; % cm;
%           2.2) measure the height of the FOV (using the ruler, see camera_fov_callibration.jpg)
%               height_ = 0.82; % cm;
%       3) get the angle
%           vfov = radtodeg(2*atan(height_/(2*distance_))) = 2.6097;
%       This example comes from a ROI of 882 x 992 [height width]
%       if using full FOV 1024 x 1280 [height width], then height_ = 0.9520 % cm
%           and vfov = 3.0296
%       if cropping FOV update calculation
%           (put ruler in focal point and measure distance between pixels)
%       if using the same camera the height_ might change depending on final zoom used with lenses.   

% Default Parameters
if ~exist('width_pix', 'var') || isempty(width_pix)
    width_pix = 480;
end

if ~exist('height_pix', 'var') || isempty(height_pix)
    height_pix = 480;
end

if ~exist('vfov', 'var') || isempty(vfov)
    vfov = deg2rad(2.15);
end

if ~exist('oDir', 'var') || isempty(oDir)
    oDir = [pwd, filesep, 'calibration-transform.dat'];
end

% delete calibration-transform
if exist(oDir, 'file')
    delete(oDir)
end

% example: lines 1-3 with z down (during an experiment). Need to round.
line1 = '-1 0 0';
line2 = '0 0 1';
line3 = '0 1 0';

%line 4 is unused, so leave as
line4 = '0 0 0';

%line 5 (after normalization). It's a vector (3 coordinates) -
PixToVec = @(x, y) [(x+.5) - width_pix*.5, ...
    (y+.5) - height_pix*.5 , ...
    height_pix * 0.5 / tan(vfov * 0.5)];
sphere_centre = PixToVec(sphere_cx_px, sphere_cy_px);
sphere_centre = sphere_centre / norm(sphere_centre);
line5 = num2str(sphere_centre);

%line 6 -
sphere_circum = PixToVec(sphere_cx_px + sphere_radius_px, sphere_cy_px);
sphere_circum = sphere_circum / norm(sphere_circum);
sphere_fov = acos(dot(sphere_centre, sphere_circum))*2;
line6 = num2str(sphere_fov);

%lines 7-8 (for black and white ball) -
line7 = '0 0 0';
line8 = '0 0 0';

%line 9 - vfov
line9 = num2str(vfov, 6);

%line 10
line10 = '0';

%write to file
fid = fopen(oDir, 'wt');

if fid == -1
    disp('wrong FicTracConfig file name')
    return
end

fprintf(fid, '%s\n', line1);
fprintf(fid, '%s\n', line2);
fprintf(fid, '%s\n', line3);
fprintf(fid, '%s\n', line4);
fprintf(fid, '%s\n', line5);
fprintf(fid, '%s\n', line6);
fprintf(fid, '%s\n', line7);
fprintf(fid, '%s\n', line8);
fprintf(fid, '%s\n', line9);
fprintf(fid, '%s\n', line10);

fclose(fid);
 
end
