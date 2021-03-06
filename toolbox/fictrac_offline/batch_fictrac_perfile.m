function batch_fictrac_perfile(FolderName, FileName, iparams)
% batch_fictrac_perfile: runs fictrac on each video 
%
% Usage:
%	batch_fictrac_perfile(FolderName, FileName, iparams)
%
% Args:
%   FolderName: name of folder to run
%   FileName: name of file to run
%   iparams: parameters to update
%       (redo: flag to redo tracking)
%           (default, 0)
%       (cDir: current directory)
%           (default, pwd)
%       (cDir: current directory)
%           (default, pwd)
%       (do_search_flag: flag to do search)
%           (default, pwd)
%       (load_template_flag: flag to load external template)
%           (default, pwd)
%       (template_fn: template image to use)
%           (default, [])
%       (oDir: directory to save figures)
%           (default, [pwd, filesep, 'fictrac'])
%       (ball_radious: radious of ball)
%           (default, 4.5 mm)
%       (hbins_speed: bins to use for histogram)
%           (default, [-20:0.1:60])
%
% Notes:
%   requires fictrac to be installed.
%       for windows use https://github.com/murthylab/fic-trac-win
%       If you are using a Point Grey/FLIR camera, make sure the FlyCapture SDK is installed. 
%           Copy FlyCapture2_C.dll from the Point Grey directory
%           (it is in the bin folder - for instance, 
%           C:\Program Files\Point Grey Research\FlyCapture2\bin64) 
%           and place it in your FicTrac directory. 
%           If it is named FlyCapture2_C_v100.dll rename it.
%           some times it requires both: FlyCapture2_v100 and FlyCapture2_C_v100
%       see https://github.com/murthylab/fly-vr for some additional
%           instalation info
%
%   requires a config_file template in the current directory
%       (generated by batch_gen_fictrac_input_files.m)
%
%   From command line fictrac can be run using the following lines
%       FicTrac-PGR FicTracPGR_config.txt
%       FicTrac FicTrac_config.txt
%
%   this setup x = yaw, y = pitch, z = roll (this needs to be tested for every setup)
%   this setup y = forward and x = lateral (this needs to be tested for every setup)
%   if this order changes then edit fictrac_txt2mat.m accordingly to match
%       the above.

% default params
fictracpars.redo = 0;
fictracpars.cDir = pwd;
fictracpars.do_search_flag = 0;
fictracpars.load_template_flag = 0;
fictracpars.template_fn = 'template_im.jpg';
fictracpars.oDir = [pwd, filesep, 'fictrac'];
fictracpars.ball_radious = 4.5;
fictracpars.hbins_speed = -20:0.1:60;

if ~exist('FolderName', 'var')
   FolderName = [];
end

if ~exist('FileName', 'var')
   FileName = [];
end

if ~exist('iparams', 'var'); iparams = []; end
fictracpars = loparam_updater(fictracpars, iparams);

fo2reject = {'.', '..', 'preprocessed', ...
    'BData', 'motcor', 'rawtiff', 'stitch'};

if ~isempty(fictracpars.oDir) && ...
        ~exist(fictracpars.oDir, 'dir')
    mkdir(fictracpars.oDir)
end

% finding folders and filtering out data that is not selected
fo2run = dir;
fo2run = str2match(FolderName, fo2run);
fo2run = str2rm(fo2reject, fo2run);
fo2run = {fo2run.name};

fprintf(['Running n-folders : ', num2str(numel(fo2run)), '\n'])

for folder_i = 1:numel(fo2run)
    
    fprintf(['Running folder : ', fo2run{folder_i}, '\n']);
    cd(fo2run{folder_i}); 
    runperfolder(FileName, fictracpars);    
    cd(fictracpars.cDir)
    
end

fprintf('... Done\n')

end

function runperfolder(fname, fictracpars)
% runperfolder: function that runs all files per directory
%
% Usage:
%   runperfolder(fname, fictracpars)
%
% Args:
%   fname: file name
%   fictracpars: input parameters

% find all videos to process
prefixstr = '.mp4';
vid2run = rdir(['*', prefixstr]);
vid2run = {vid2run.name}';

if isempty(vid2run)
    prefixstr = '.avi';
    vid2run = rdir(['*', prefixstr]);
    vid2run = {vid2run.name}';
    vid2run = str2rm({'-debug.avi'}, vid2run);
end

vid2run = str2match(fname, vid2run);
vid2run = strrep(vid2run, ['_vid', prefixstr], '');

fprintf(['Running n-videos : ', num2str(numel(vid2run)), '\n'])

% default config file name
config_file = 'FicTrac_config.txt';
temp_config_file = 'temp_config.txt';

for i = 1:numel(vid2run)
    
    % run fictrac
    %   this generates the following outout files:
    %       *_fictrac.txt'
    %       *_maskim.tiff'
    %       *_calibration-transform.dat'
    %       *_template.jpg'
    
    if ~exist([vid2run{i}, '_fictrac.txt'], 'file') || ...
            fictracpars.redo
        % edit config file to include file name
        edit_FicTrac_config(config_file, ...
            vid2run{i}, fictracpars.do_search_flag, ...
            fictracpars.load_template_flag, ...
            fictracpars.template_fn, prefixstr)

        command2run = ['FicTrac ', temp_config_file];

        % execute fictrac
        if ispc
            dos(command2run);
        else
            unix(command2run);
        end

        % delete temp config file
        delete(temp_config_file)
    end
    
    % add fictrac variable to '_vDat.mat' file
    fictrac_txt2mat([vid2run{i}, '_fictrac.txt'])
    
    % plot basic fictrac results
    fictrac_plot_results([vid2run{i}, '_vDat.mat'], fictracpars)
    
end

end

function edit_FicTrac_config(config_file, ...
    vid_name, do_search_flag, load_template_flag, ...
    template_name, prefixstr)
% edit_FicTrac_config: edit config_file to run each video file
%
% Usage:
%	edit_FicTrac_config(config_file, ...
%       vid_name, do_search_flag, load_template_flag, ...
%       template_name)
%
% Args:
%   config_file: config file to use as template
%   vid_name: video name
%   do_search_flag: flag to do search
%   	(default, pwd)
%   load_template_flag: flag to load external template
%       (default, pwd)
%   template_fn: template image to use
%       (default, [])
%   prefixstr: video prefix

% load text file
fid = fopen(config_file, 'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;

while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

% replace field
A{2} = ['input_vid_fn        ', vid_name, '_vid' prefixstr];
A{3} = ['output_fn           ', vid_name, '_fictrac.txt'];
A{4} = ['mask_fn             ', vid_name, '_maskim.tiff'];
A{5} = ['transform_fn        ', vid_name, '_calibration-transform.dat'];

if load_template_flag
    A{6} = ['template_fn         ', template_name];
    A{7} = 'load_template        1';
    A{37} = 'do_update            0';
else
    A{6} = ['template_fn         ', vid_name, '_template.jpg'];
end

if do_search_flag
    A{26} = 'do_search           1';
end

% save txt file
fid = fopen('temp_config.txt', 'w');

for i = 1:numel(A)
    
    if A{i+1} == -1
        fprintf(fid, '%s', A{i});
        break
    else
        fprintf(fid, '%s\n', A{i});
    end
    
end

fclose(fid);

end

function fictrac_plot_results(fname, fictracpars)
% fictrac_plot_results: function to plot results
%
% Usage:
%	fictrac_plot_results(fname, fictracpars)
%
% Args:
%   fname: name if file containing 'fictDat'
%   fictracpars: video name
%   fictracpars: input parameters
%
% Notes
%   this setup x = yaw, y = pitch, z = roll (this needs to be tested for every setup)
%   this setup y = forward and x = lateral (this needs to be tested for every setup)

% load fictDat
load(fname, 'fictDat')

% plot regular motion stats
figH = figure('Position', genfigpos(1, 'center', [900 600]));
axH(1) = subplot(2, 2, 1);
axH(2) = subplot(2, 2, 2);
axH(3) = subplot(2, 2, 3:4);

x_mm = fictDat.var(:, 15)*fictracpars.ball_radious;
y_mm = fictDat.var(:, 16)*fictracpars.ball_radious;
speed_mm_s = fictDat.var(:, 19)*fictracpars.ball_radious;
lv_mm_s = diff(x_mm);
fv_mm_s = diff(y_mm);
time_stamps = 1:numel(fictDat.var(:, 1));
    
% plot XY
plot(x_mm, y_mm, 'Color', 'k', 'Parent', axH(1))

[y_, ~] = hist(speed_mm_s, fictracpars.hbins_speed);

lineH(1) = plot(fictracpars.hbins_speed, y_/sum(y_), ...
    'Color', 'b', 'Parent', axH(2));

hold(axH(2), 'on')

[y_, ~] = hist(fv_mm_s, fictracpars.hbins_speed);

lineH(2) = plot(fictracpars.hbins_speed, y_/sum(y_), ...
    'Color', 'r', ...
    'Parent', axH(2));

[y_, ~] = hist(lv_mm_s, fictracpars.hbins_speed);

lineH(3) = plot(fictracpars.hbins_speed, y_/sum(y_), ...
    'Color', 'm', ...
    'Parent', axH(2));

% plot instantaneous speed
plot(time_stamps, speed_mm_s, 'Color', 'b', ...
    'Parent', axH(3))
hold(axH(3), 'on')

% plot instantaneous lateral speed
plot(time_stamps(2:end), lv_mm_s, 'Color', 'm', ...
    'Parent', axH(3))

% plot instantaneous forward speed
plot(time_stamps(2:end), fv_mm_s, 'Color', 'r', ...
    'Parent', axH(3))

axH(1).Title.String = 'Distance traveled';
axH(1).XLabel.String = 'X distance (mm)';
axH(1).YLabel.String = 'Y distance (mm)';

axH(2).Title.String = 'Distribution of Speed/velocity';
axH(2).XLabel.String = 'speed (mm/s)';
axH(2).YLabel.String = 'Probability';

axH(3).Title.String = 'Speed/velocity over time';
axH(3).XLabel.String = 'Time (frames)';
axH(3).YLabel.String = 'Instantaneous speed (mm/frame)';
axH(3).XLim = [1 numel(time_stamps)];

legend(axH(3), lineH, {'speed', 'f-velocity', 'l-velocity'}, 'Location', 'northeast')

savefig_int(figH, fictracpars.oDir, [strrep(fname, '_vDat.mat', ''), '_fictrac_stats'], ...
    [1 0 0 0 0 0 0 0 1])
close(figH)

end
