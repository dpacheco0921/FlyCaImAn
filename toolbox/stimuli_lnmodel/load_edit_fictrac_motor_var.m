function [Y, imtime, vidtime] = ...
    load_edit_fictrac_motor_var(wDat, ball_radious, ...
    var2gen, sig, siz, tres_flag)
% load_edit_fictrac_motor_var: gather motor variables to use 
%   (resample to imaging sampling rate), smooth, and resample
%
% Usage:
%   [Y, imtime, vidtime] = load_edit_fictrac_motor_var(wDat, ...
%       ball_radious, var2gen, sig, siz, tres_flag)
%
% Args:
%   wDat: metadata variable
%   ball_radious: radious of the ball (mm)
%   var2gen: cell array of motor variables to use, any of 'speed', 'fV',
%       'lV', 'yaw', 'pitch', 'roll' (an unknown name gives an error)
%   sig: std of gaussian kernel (default, [] = no smoothing)
%   siz: size of kernel (default, [] = no smoothing)
%   tres_flag: flag to resample to imaging time resolution (default, 0)
%
% Returns:
%   Y: one row per entry of var2gen [numel(var2gen), T]. Ground velocities
%       (speed, fV, lV) are in mm/s. Rotational rows are the fictrac
%       per-frame rotation times the frame rate (called *_deg_s in the
%       code, plotted as rad/s elsewhere; unit not verified).
%       T = video frames, or imaging frames if tres_flag.
%   imtime: imaging frame times in seconds, rounded to 10 ms
%   vidtime: video frame times in seconds
%
% Notes:
%   this setup x = yaw, y = pitch, z = roll (this needs to be tested for every setup)
%   this setup y = forward and x = lateral (this needs to be tested for every setup)

if ~exist('sig', 'var')
    sig = [];
end

if ~exist('siz', 'var')
    siz = [];
end

if ~exist('tres_flag', 'var') || isempty(tres_flag)
    tres_flag = 0;
end

framerate = round(1/mean(diff(wDat.vid.fstEn{1}(:, 1))));
imtime = wDat.fTime;   
imtime = round(imtime*100)/100;
vidtime = wDat.vid.fstEn{1}(:, 2);

% calculate velocity (forward/lateral)
% in this setup y = forward
fv_mm_s = diff(wDat.vid.var{1}(:, 16)*ball_radious*framerate);
fv_mm_s = [fv_mm_s(1); fv_mm_s];

% in this setup x = lateral
lv_mm_s = diff(wDat.vid.var{1}(:, 15)*ball_radious*framerate);
lv_mm_s = [lv_mm_s(1); lv_mm_s];

% calculate rotational velocity (pitch, yaw, roll)
% in this setup x = yaw, y = pitch, z = roll
yaw_deg_s = (wDat.vid.var{1}(:, 6))*framerate;
pitch_deg_s = (wDat.vid.var{1}(:, 7))*framerate;
roll_deg_s = (wDat.vid.var{1}(:, 8))*framerate;

% smooth traces
if ~isempty(sig) && ~isempty(siz)
   fv_mm_s = imblur(fv_mm_s, sig, siz, 1); 
   lv_mm_s = imblur(lv_mm_s, sig, siz, 1); 
   yaw_deg_s = imblur(yaw_deg_s, sig, siz, 1); 
   pitch_deg_s = imblur(pitch_deg_s, sig, siz, 1); 
   roll_deg_s = imblur(roll_deg_s, sig, siz, 1); 
end

% compute extra: speed and aceleration
speed_mm_s = sqrt(fv_mm_s.^2 + lv_mm_s.^2);

% generate locomotor variables
for i = 1:numel(var2gen)

    % track whether this name matched, otherwise the row would silently
    %   stay zero (or be missing if it is the last one)
    matched = false;

    if contains(var2gen{i}, 'speed')
        Y(i, :) = speed_mm_s;
        matched = true;
    end

    if contains(var2gen{i}, 'fV')
        Y(i, :) = fv_mm_s;
        matched = true;
    end

    if contains(var2gen{i}, 'lV')
        Y(i, :) = lv_mm_s;
        matched = true;
    end

    if contains(var2gen{i}, 'yaw')
        Y(i, :) = yaw_deg_s;
        matched = true;
    end


    if contains(var2gen{i}, 'pitch')
        Y(i, :) = pitch_deg_s;
        matched = true;
    end


    if contains(var2gen{i}, 'roll')
        Y(i, :) = roll_deg_s;
        matched = true;
    end

    if ~matched
        error('load_edit_fictrac_motor_var:unknownVar', ...
            'Unknown motor variable: %s', var2gen{i});
    end

end

if tres_flag
    
    % resample to imaging resolution defined in wDat.fTime
    Y = interp1(vidtime, Y', imtime, 'linear');
    Y = framegapfill(find(isnan(Y(:, 1))), Y');
    
end

end
