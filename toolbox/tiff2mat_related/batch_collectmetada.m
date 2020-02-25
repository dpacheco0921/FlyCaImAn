function batch_collectmetada(FolderName, FileName, iparams)
% batch_collectmetada: collect metadata from *.mat files (from Ephys/Stimuli
%   delivery PC) generated by prv or LEDcontroler
%
% Usage:
%   batch_collectmetada(FolderName, FileName, iparams)
%
% Args:
%   FolderName: name of folders to load
%   FileName: name of files to load
%   iparams: parameters to update
%       (Ygalvo_Ch: Y-axis of galvo or resonant channel)
%           (default, 2)
%       (stimcopy_Ch: channel receiving a copy of stimuli delivered in bin file)
%           (default, 1)
%       (stim_oCh: analog ouput channel use to deliver stimuli)
%           (default, 1)
%       (buffer: time to chop from recording end 0 (ms))
%           (default, 0)
%       (pgate: plot frame width time for each file per folder)
%           (default, 0)
%       (pgates: plot raw Y trace, start and end of each frame)
%           (default, 0)
%       (cDir: current directory)
%       (fo2reject: folders to reject)
%       (fi2reject: files to reject)
%       (fsuffix: suffix of files with Data variable)
%           (default, '_rawdata.mat')
%       (fmetsuffix: suffix of metadata file)
%           (default, '_metadata.mat')
%       (mode: two modes:
%           (0, readout resuls, default)
%           (1, generate and save)
%       (minframet: minimun frame duration in tens of miliseconds)
%           (default, 900 (90 ms))
%       (findstim: gate to use stimuli trace to find stimuli start and end)
%           (default, 0)
%       (minwidth: minimun width of stimuli (ms))
%           (default, [])
%       (stimths: voltage threshold to find stimuli)
%           (default, 0)      
%       (vid_time_init: initial time to use for mean substration)
%           (default, 100 timepoints)      
%       (vid_volt_ths: voltage threshold to find pulses)
%           (default, [2.5 2.5])      
%       (vid_ttl_width: expected width of pulses)
%           (default, 20 timepoints)      
% 
% Notes
% this function updates: lStim and iDat
% lStim: stimuli related (auditory/opto) metadata structure
%   adds the following fields:
%       lStim.trace: stimulus raw trace
%       lStim.lstEn: stimulus onset and offset  (sound and opto)
%       lStim.sPars.order: stimulus order (sound and opto)
%       lStim.sPars.name: name of each stimulus (sound and opto)
%       lStim.sPars.int: stimulus intensity (sound and opto)
%       lStim.sPars.sr: sampling rate (sound and opto)
%       lStim.sPars.basPre: silence pre stimulus (sound and opto)
%       lStim.sPars.basPost: silence post stimulus (sound and opto)
%       lStim.sPars.freq: frequency of pulse stimulus (opto)
%       lStim.sPars.width: width of pulse stimulus (opto)
%
% iDat: image metadata structure
%   adds the following fields:
%       iDat.StackN: updates depending on the length of time stamps
%           collected
%       iDat.fstEn: frame timestamps [onset offset]
%
% 2019-02-22: bleedthrough is adding a stronger DC change than expected and
%   affecting Frame time estimation, so now it runs calculate the mode frame
%   width and then uses it to eliminate peaks in between. need to work on
%   denoising this opto-related signal
% 2019-07-26:
%   1) now it is compatible with cases where no Y-galvo info is
%   provided or extra stimuli info
%   2) compatible with new setup (reads Y-galvo from data_*.mat files)

% default params
metpars.Ygalvo_Ch = 2;
metpars.stimcopy_Ch = 1;
metpars.stim_oCh = 1;
metpars.buffer = 0;
metpars.pgate = 0;
metpars.pgates = 0;
metpars.cDir = pwd;
metpars.fo2reject = {'.', '..', 'preprocessed', 'BData'};
metpars.fi2reject = [];
metpars.fsuffix = '_rawdata.mat';
metpars.fmetsuffix = '_metadata.mat';
metpars.mode = 1;
metpars.minframet = 900;
metpars.findstim = 0;
metpars.minwidth = [];
metpars.stimths = 0;
metpars.vid_ttl_width = 20;
metpars.vid_time_init = 200;
metpars.vid_volt_ths = [2 2];

% update variables
if ~exist('FolderName', 'var'); FolderName = []; end
if ~exist('FileName', 'var'); FileName = []; end
if ~exist('iparams', 'var'); iparams = []; end
metpars = loparam_updater(metpars, iparams);

% find folders
fo2run = dir;
fo2run = str2match(FolderName, fo2run);
fo2run = str2rm(metpars.fo2reject, fo2run);
fo2run = {fo2run.name};

fprintf(['Running n-folders : ', num2str(numel(fo2run)), '\n'])

for i = 1:numel(fo2run)
    
    fprintf(['Running folder : ', fo2run{i}, '\n']);
    cd(fo2run{i}); 
    runperfolder(FileName, fo2run{i}, metpars);
    cd(metpars.cDir)
    
end

fprintf('... Done\n')

end

function runperfolder(fname, foname, metpars)
% runperfolder: function collects timestamps and stimulus from bin file 
%   with "metpars.fsuffix" suffix and fDat info
%
% Usage:
%   runperfolder(fname, foname, metpars)
%
% Args:
%   fname: file name template string
%   foname: figure name
%   metpars: parameters

% directory params
fi2run = rdir(['.', filesep, '*', metpars.fsuffix]);
fname = addsuffix(fname, metpars.fsuffix);
fi2run = str2match(fname, fi2run);
fi2run = str2rm(metpars.fi2reject, fi2run);
fi2run = {fi2run.name};
fi2run = strrep(fi2run, ['.', filesep], '');
fi2run = strrep(fi2run, '_rawdata.mat', '');

% plot times
if metpars.pgate == 1
    metpars.figH = figure('Name', foname);
    metpars.AxH = subplot(1, 1, 1);
end

for F2R_idx = 1:numel(fi2run)
    
    display(['Running file : ', fi2run{F2R_idx}])   
    runperfile(fi2run{F2R_idx}, metpars)
    
end

fprintf('****** Done ******\n')

end

function runperfile(filename, metpars)
% runperfolder: function collects timestamps and stimulus from bin file 
% with "metpars.fsuffix" suffix and fDat info
%
% Usage:
%   runperfolder(filename, ip)
%
% Args:
%   filename: name of file to load
%   metpars: parameters

% Load main variables 'lStim', 'iDat', 'fDat'
display(['Running file : ', filename])
load(['.', filesep, filename, '_metadata.mat'], ...
    'lStim', 'iDat', 'fDat')

% copy datatype (2D or 4D)
datatype = fDat.DataType;

% input data to load
stim_file2load = [];
if contains(datatype, 'song') || contains(datatype, 'prv') ...
        && ~contains(datatype, 'fict')
    
    stim_file2load = 'prv';
    metpars.Ygalvo_Ch = 2;
    metpars.stimcopy_Ch = 1;
    
    if contains(datatype, 'snd_odor')
        metpars.Ygalvo_Ch = 3;
    elseif contains(datatype, 'diego_noball')
        metpars.stimcopy_Ch = [];
    elseif contains(datatype, 'diego_FlyOnTheBall')
        metpars.stimcopy_Ch = [3 4];
    end    
    
    if contains(datatype, 'opto')
        metpars.stim_oCh = 2;
    end
    
elseif contains(datatype, 'opto') && ~contains(datatype, 'prv') ...
        && ~contains(datatype, 'fict')
    
    stim_file2load = 'LEDcontroler';
    
elseif contains(datatype, 'fict')
    
    stim_file2load = 'fict';
    metpars.Ygalvo_Ch = 1;
    metpars.stimcopy_Ch = 2;
    
else
    
   stim_file2load = 'nostim';
    
end

% Notes:
% cases for 3DxT new data, sometimes the Y movement stops, giving you a
%   wrong frame end for the last frame, which is then deleted

% collect galvo and stimuli copy channel(s)
Ch = load_galvo_and_stim_channels(...
    stim_file2load, filename);

% if metadata iDat variable already has framerate information, then use it
if isfield(iDat, 'framerate') && ...
        ~isempty(iDat.framerate)
    % iDat.framerate is in Hzs
    metpars.minframet = floor(10000/iDat.framerate);
end

% collect timestamps using galvo channels
frame_num = iDat.FrameN*iDat.StackN;
iDat.fstEn = collect_frame_timestamps(Ch, ...
    metpars.Ygalvo_Ch, metpars.minframet, ...
    frame_num, stim_file2load);

% plot stim trace
if metpars.pgates == 1 && ...
        ~isempty(Ch) && ~isempty(iDat.fstEn)

    figure('Name', filename);
    AxHs = subplot(1, 1, 1);
    plot(Ch(metpars.Ygalvo_Ch, :), 'Parent', AxHs);
    hold(AxHs, 'on');
    plot(iDat.fstEn(:, 1), ...
        Ch(metpars.Ygalvo_Ch, iDat.fstEn(:, 1)), ...
        'o', 'Parent', AxHs)
    plot(iDat.fstEn(:, 2), ...
        Ch(metpars.Ygalvo_Ch, iDat.fstEn(:, 2)), ...
        'go', 'Parent', AxHs)
    xlabel(AxHs, 'Time 0.1 ms');
    ylabel(AxHs, 'V')

end

% update lStim field 'fName'
if contains(datatype, 'song') || contains(datatype, 'prv')
    lStim.fName = [lStim.fName, '_bin.mat'];
elseif contains(datatype, 'opto') && ~contains(datatype, 'prv')
    lStim.fName = [lStim.fName, '.bin'];
else
    lStim.fName = 'no stimulus file';
end

if ~isempty(iDat.fstEn)
    fprintf(['first frame ', num2str(iDat.fstEn(1, :)), ...
        ' second ', num2str(iDat.fstEn(2, :)), '\n'])
end

% homogenize frames and frametimes
%   (chops either image stacks or stimuli vectors)
%   edits iDat fields iDat.fstEn, iDat.FrameN & iDat.StackN
if ~isempty(iDat.fstEn)
    iDat = homogenize_frame_frametime(...
        filename, frame_num, iDat, ...
        datatype, metpars.mode);
end

% plot frame-diff
if metpars.pgate == 1 && ...
        ~isempty(iDat.fstEn)

    plot(iDat.fstEn(:, 2) - iDat.fstEn(:, 1), ...
        'parent', metpars.AxH);
    hold(metpars.AxH,'on')
    xlabel(metpars.AxH, 'Frame');
    ylabel(metpars.AxH, 'Frame width 0.1 ms')

end

% define time before 1st frame and last frame to include
%   (default is just from 1st to last frame, but it could include some buffer
%   see metpars.buffer)
if ~isempty(iDat.fstEn)
    
    start_end = [iDat.fstEn(1, 1) - metpars.buffer*lStim.fs, ...
        (iDat.fstEn(end, 2) + metpars.buffer*lStim.fs)];
    % set first frame start to index 1
    iDat.fstEn = iDat.fstEn - start_end(1) + 1;
    
else
    
    start_end = [1 numel(Ch(1, :))];
    iDat.fstEn = [1 numel(Ch(1, :))];
    
end

% load video related variables (frames and fictrac vars)
if exist([filename, '_vid.mp4'], 'file')
    
    vidDat = extract_video_timestamps(filename, Ch, ...
        metpars.vid_time_init, metpars.vid_volt_ths, ...
        metpars.vid_ttl_width, metpars.pgate, start_end);
        
end

% chop and pass Ch (stimuli trace) to lStim.trace
if ~isempty(Ch)
    Ch = Ch(metpars.stimcopy_Ch, :);

    if start_end(2) > numel(Ch)
        lStim.trace = Ch(start_end(1):end);
        lStim.trace = [lStim.trace, zeros(1, start_end(2)-numel(Ch))];
    else
        lStim.trace = Ch(start_end(1):start_end(2));
    end

else
    lStim.trace = zeros(1, iDat.fstEn(end, 2));          
end

% get stim onset and offset and extra metadata
if contains(stim_file2load, 'prv')

    % load rDat
    load(['.', filesep, filename, '_vDat.mat'], 'rDat')

    [lStim.lstEn, lStim.trace] = ...
        songstim_init_end_rDat(rDat, start_end, ...
        metpars.stim_oCh, metpars.findstim, ...
        metpars.minwidth, metpars.stimths);

    % collect extra metadata
    lStim.sPars.order = rDat.stimOrder;
    lStim.sPars.name = rDat.ctrl.stimFileName;

    all_int = cell2mat(rDat.ctrl.intensity(:, 1));

    if contains(datatype, 'song')
        lStim.sPars.int = all_int(:, 1);
    elseif contains(datatype, 'opto')
        lStim.sPars.int = all_int(:, 2);
    end

    lStim.sPars.sr = rDat.ctrl.rate;
    lStim.sPars.basPre = rDat.ctrl.silencePre;
    lStim.sPars.basPost = rDat.ctrl.silencePost;
    clear rDat

elseif contains(stim_file2load, 'LEDcontroler')
    
    load(['.', filesep, filename, '.mat'], 'sDat');
    lStim.lstEn = optostim_init_end(lStim.trace, ...
        sDat);

    % collect extra metadata
    lStim.sPars.freq = ...
        repmat(sDat.freq, [size(lStim.lstEn, 1), 1]); 
    lStim.sPars.width = ...
        repmat(sDat.width, [size(lStim.lstEn, 1), 1]); 
    lStim.sPars.int = ...
        repmat(sDat.intensity, [size(lStim.lstEn, 1), 1]); 
    lStim.sPars.sr = ...
        repmat(sDat.fs, [size(lStim.lstEn, 1), 1]); 
    lStim.sPars.basPre = ...
        repmat(sDat.silencePre, [size(lStim.lstEn, 1), 1]); 
    lStim.sPars.basPost = ...
        repmat(sDat.silencePost, [size(lStim.lstEn, 1), 1]);

    lStim.sPars.order = 1:size(lStim.lstEn, 1);

    if isfield(sDat, 'stimFileName')
        lStim.sPars.name = ...
            repmat({sDat.stimFileName}, [size(lStim.lstEn, 1), 1]);
    else
        lStim.sPars.name = ...
            repmat({'OPTO'}, [size(lStim.lstEn, 1), 1]);                
    end

    clear sDat
    
elseif contains(stim_file2load, 'fict')
    
    % read text file with stimulus info
    sDat = parse_opto_fictrac_txtfile(...
        ['.', filesep, filename, '.txt']);
    
    lStim.lstEn = find_stim_int(lStim.trace, ...
        metpars.minwidth, metpars.stimths);
        
    % collect extra metadata
    lStim.sPars.freq = []; 
    lStim.sPars.width = []; 
    lStim.sPars.int = sDat.intensity(1:size(lStim.lstEn, 1)); 
    lStim.sPars.sr = sDat.rate(1:size(lStim.lstEn, 1)); 
    lStim.sPars.basPre = sDat.silencePre(1:size(lStim.lstEn, 1)); 
    lStim.sPars.basPost = sDat.silencePost(1:size(lStim.lstEn, 1));
    lStim.sPars.order = 1:size(lStim.lstEn, 1);

    if isfield(sDat, 'stimFileName')
        lStim.sPars.name = ...
            sDat.stimFileName(1:size(lStim.lstEn, 1));
    else
        lStim.sPars.name = ...
            repmat({'OPTO'}, [size(lStim.lstEn, 1), 1]);                
    end

    clear sDat
    
else

    % generate arbitrary stimulus information
    lStim.sPars.freq = 1;
    lStim.sPars.width = 1; 
    lStim.sPars.int = 0; 
    lStim.sPars.sr = 10^4; 
    lStim.sPars.basPre = 0; 
    lStim.sPars.basPost = 0;
    lStim.sPars.order = 1;
    lStim.sPars.name = {'nostim'};
    lStim.lstEn = [1 2];

end

if metpars.mode
    save(['.', filesep, filename, '_metadata.mat'], ...
        'iDat', 'lStim', '-append')
end

if exist('vidDat', 'var')
    save(['.', filesep, filename, '_metadata.mat'], ...
        'vidDat', '-append')
end

clear FrameInit FrameEnd order2chop Data NewStart lStim

clear iDat fDat Ch lStim ROI lStim

end

function vidDat = extract_video_timestamps(filename, Ch, ...
    vid_time_init, vid_volt_ths, vid_ttl_width, plot_flag, start_end)
% extract_video_timestamps: get frame timestamps from trigger and shutter
%
% Usage:
%   vidDat =  extract_video_timestamps(filename, Ch, ...
%       vid_time_init, vid_volt_ths, vid_ttl_width, plot_flag)
% 
% Args:
%   filename: file name
%   Ch: loaded channels
%   metpars: metadata parameters
%   vid_time_init: initial time to use for mean substration
%   vid_volt_ths: voltage threshold to find pulses
%   vid_ttl_width: expected width of pulses
%   plot_flag: flag to plot results
%   start_end: start and end timepoint to use

% load prv-related variable rDat and fictrac-related variable "fictDat"
vidDat = [];
try
    load(['.', filesep, filename, '_vDat.mat'], ...
        'rDat', 'fictDat')
    vidDat = fictDat;
    clear fictDat
end

% readout video frames info
vid_obj = VideoReader([filename, '_vid.mp4']);
read(vid_obj, inf);
vid_frame_num = vid_obj.NumberOfFrames;

% readout trigger and shutter channels
if sum(rDat.vidframeCh) || ...
        sum(rDat.vidtriggerCh)
    
    fprintf('video related channels found\n')

    % extract triggered frames
    trigger_onset_offset = ...
        get_ttl_start_end(Ch(rDat.vidtriggerCh, :), ...
        vid_volt_ths(2), vid_time_init, ...
        vid_ttl_width, rDat.Fs);

    % extract recorded frames        
    shutter_onset_offset = ...
        get_ttl_start_end(Ch(rDat.vidframeCh, :), ...
        vid_volt_ths(1), vid_time_init, ...
        vid_ttl_width, rDat.Fs);

end

% print comparison
fprintf(['frames in mp4: ', num2str(vid_frame_num), '\n'])
fprintf(['frames triggered: ', ...
    num2str(size(trigger_onset_offset, 1)), '\n'])
fprintf(['frames readout (from shutter opening): ', ...
    num2str(size(shutter_onset_offset, 1)), '\n'])

% print video frame rate
fprintf(['frame rate from trigger: ', ...
    num2str(rDat.Fs/mean(diff(trigger_onset_offset(:, 1)))), '\n'])
fprintf(['frame rate from readout (from shutter opening): ', ...
    num2str(rDat.Fs/mean(diff(shutter_onset_offset(:, 1)))), '\n'])

% chop times outside frame timpoints
l2chop = find(shutter_onset_offset(:, 1) > start_end(2));

if ~isempty(l2chop)
    shutter_onset_offset = ...
        shutter_onset_offset(1:(l2chop - 1), :);
    trigger_onset_offset = ...
        trigger_onset_offset(1:(l2chop - 1), :);
end

% make timpoints units relative to start point
shutter_onset_offset = ...
    shutter_onset_offset - start_end(1) + 1;
trigger_onset_offset = ...
    trigger_onset_offset - start_end(1) + 1;

if isfield(vidDat, 'var') && ~isempty(vidDat.var)
    fprintf('fictrac file found\n')

    % reducing timestamps to actual number of frames processed
    if size(trigger_onset_offset, 1) < vid_frame_num || ...
            size(shutter_onset_offset, 1) < vid_frame_num
        
        % chop vidDat.var
        idx2use = find(vidDat.var(:, 1) <= min([...
            size(trigger_onset_offset, 1), ...
            size(trigger_onset_offset, 1)]));
        vidDat.var = vidDat.var(idx2use, :);
        
    end
    
    shutter_onset_offset = ...
        shutter_onset_offset(vidDat.var(:, 1), :);
    trigger_onset_offset = ...
        trigger_onset_offset(vidDat.var(:, 1), :);

end

if plot_flag && ...
        ~isempty(shutter_onset_offset) && ...
        ~isempty(trigger_onset_offset)
    
    % display camera
    figure()
    axH = subplot(1, 1, 1);
    plot(Ch(rDat.vidtriggerCh, 1:2*10^4), 'k', 'Parent', axH)
    hold(axH, 'on')
    plot(Ch(rDat.vidframeCh, 1:2*10^4), 'b', 'Parent', axH)

    plot(shutter_onset_offset(:, 1), ...
            Ch(rDat.vidframeCh, shutter_onset_offset(:, 1)), ...
            'o', 'Parent', axH)
    plot(shutter_onset_offset(:, 2), ...
        Ch(rDat.vidframeCh, shutter_onset_offset(:, 2)), ...
        'go', 'Parent', axH)
    axH.XLim = [0 2*10^4];
    
end

% append shutter times
vidDat.fstEn = shutter_onset_offset;

end

function fstEn = collect_frame_timestamps(Ch, ...
    Ygalvo_Ch, minframet, frame_num, stim_file2load)
% collect_frame_timestamps: function that estimates timestamps from galvo
%   channel readout.
%
% Usage:
%   fstEn = collect_frame_timestamps(Ch, ...
%     Ygalvo_Ch, minframet, frame_num, stim_file2load)
% 
% Args:
%   Ch: channels extracted from *.bin / *.mat files
%   Ygalvo_Ch: galvo channel
%   minframet: minimun frame time
%   frame_num: frames imaged (obtained from tif)
%   stim_file2load: string with setup/data acquisition function info

% collect timestamps
if ~isempty(Ch) && ~isempty(frame_num)

    [FrameInit, FrameEnd] = ...
        colecttimestamp(Ch(Ygalvo_Ch, :), ...
        frame_num, minframet, stim_file2load);
    
else

    if isempty(frame_num)
        % no images were collected
        fprintf('No stacks imaged generating empty fields\n')
        FrameInit = [];
        FrameEnd = [];
        
    else
        % images but not galvo channel was collected
        fprintf('No galvo info provided generating arbitrary timestamps\n')

        % generate arbitrary timestamps
        %   (in case you dont have a readout of the Y galvo)
        FrameInit = (0:(frame_num-1))*10^3;
        FrameInit(1) = 1;
        FrameEnd = FrameInit + 900;
        
    end
    
end

fstEn = [FrameInit' FrameEnd'];

end

function Ch = load_galvo_and_stim_channels(...
    stim_file2load, filename)
% load_galvo_and_stim_channels: function that loads galvo and stimulus copy
%   channel(s)
%
% Usage:
%   Ch = load_galvo_and_stim_channels(...
%       stim_file2load, filename)
%
% Args:
%   stim_file2load: string with setup/data acquisition function info
%   filename: name of file to load

if contains(stim_file2load, 'prv')
    
    % stimuli delivered using prv code
    bin_file_name = ['.', filesep, filename, '_bin.mat'];

    if exist(bin_file_name, 'file')

        load(bin_file_name, 'data', 'dataScalingFactor')
        Ch = data';
        Ch = double(Ch)/dataScalingFactor;
        clear data

    else

        fprintf('No stimulus & timestamp-related bin file found\n')
        Ch = [];

    end
    
    clear bin_file_name
    
elseif contains(stim_file2load, 'LEDcontroler')

    % stimuli delivered using LEDcontroler (old setup)
    Ch = local_binread(lStim);

elseif contains(stim_file2load, 'fict')

    % stimuli delivered using **
    data = h5load(['.', filesep, filename, '.h5']);
    eval(['Ch = ', 'data.input.samples(:,3);']);
    Ch = double(Ch)';
    clear data
    
    if contains(datatype, 'opto')

        % read stimulus from *.h5 file
        data = h5load(['.', filesep, filename, '.h5']);
        eval(['Ch(2, :) = ', 'data.input.samples(:, 2);']);
        clear data
    
    else
        Ch(2, :) = zeros(size(Ch));
    end

else
    fprintf('No stimulus & timestamp-related bin file found\n')
    Ch = [];
end

end

function iDat = homogenize_frame_frametime(...
    filename, frame_n, iDat, datatype, imode)
% homogenize_frame_frametime: chop frametimes or 
%   Data depending on length, correct time stamp 
%   for last frame.
% it edits iDat.fstEn, iDat.FrameN, iDat.StackN
%
% Usage:
%   homogenize_frame_frametime(frame_n, iDat, datatype, imode, Data)
%
% Args:
%   filename: name of file to load
%   frame_n: number of frames
%   iDat: image metadata variable
%   datatype: type of data
%   imode: two modes
%       (0, readout resuls)
%       (1, generate and save)

if imode
    % load data
    load(['.', filesep, filename, '_rawdata.mat'], 'Data');
end

if frame_n < size(iDat.fstEn, 1)

    % chop stamps
    fprintf('Chopping time stamps\n')
    iDat.fstEn = iDat.fstEn(1:frame_n, :);

elseif frame_n > size(iDat.fstEn, 1)

    % chop data
    if imode

        if contains(datatype, '2DxT')

            fprintf('Chopping Data matrix\n')
            Data = Data(:, :, 1:size(iDat.fstEn, 1), :);
            iDat.FrameN = 1; 
            iDat.StackN = size(Data, 3);

        elseif contains(datatype, '3DxT')

            if frame_n - size(iDat.fstEn, 1) < iDat.FrameN*0.5

                fprintf('increase fstEn\n')
                
                % extrapolate the frame time
                for fi = (size(iDat.fstEn, 1) + 1):frame_n
                    iDat.fstEn(fi, :) = iDat.fstEn(end, :) + ...
                        iDat.fstEn(end, :) - iDat.fstEn(end-1, :);
                end

            else
                
                fprintf('Chopping Data matrix\n')
                % prune n-volume
                z_end = floor((size(iDat.fstEn, 1))/iDat.FrameN);
                Data = Data(:, :, :, 1:z_end, :);
                iDat.FrameN = size(Data, 3);
                iDat.StackN = size(Data, 4);
                frame_n = iDat.FrameN*iDat.StackN;
                
                % prune fstEn
                iDat.fstEn = iDat.fstEn(1:frame_n, :);
                
            end

        end

        save(['.', filesep, filename, '_rawdata.mat'], 'Data', '-v7.3')

    end

elseif frame_n == size(iDat.fstEn, 1)

    % sometimes last frame ending is not counted right
    iDat.fstEn(end, 2) = iDat.fstEn(end, 1) + ...
        iDat.fstEn(end-1, 2) - iDat.fstEn(end-1, 1);

end

end

function [Frame_Init, Frame_End] = ...
    colecttimestamp(mirror_y_trace, ...
    frame_num, min_frame_time, stim_file2load)
% colecttimestamp: calculate frame onset anf offset based on the movement
% of the Y galvo.
%
% Usage:
%   [Frame_Init, Frame_End] = ...
%       colecttimestampNew(mirror_y_trace, frame_num, min_frame_time)
%
% Args:
%   mirror_y_trace: mirror y trace from 2PM
%   frame_num: number of frames (from Data)
%   min_frame_time: minimun time between frames
%   stim_file2load: stimuli delivery setup
%
% Notes:
% regular vs resonant galvo

if size(mirror_y_trace, 1) > 1; mirror_y_trace = mirror_y_trace'; end

% edit trace for resonant galvo
if contains(stim_file2load, 'fict')
    idx_start = find(mirror_y_trace < -1, 1, 'last');
    mirror_y_trace(1, 1:idx_start) = mirror_y_trace(idx_start + 1);
end

% calculate the diff of the smooth vector
tracet = -diff(diff(smooth(mirror_y_trace(1, :), 10)));

% amplitude threshold
max_pred = prctile(tracet, 99.8)*0.5;

% find frame ends & delete peaks that are shorter than expected frame interval
[~, Frame_End] = findpeaks(tracet, 'MinPeakHeight', max_pred);
Frame_End(diff(Frame_End) < min_frame_time*0.8) = [];
Frame_End = Frame_End';

% correct for artifact peak at the beginning
if Frame_End(1) > min_frame_time*0.8
    % frame end is as expected
    Frame_End = Frame_End(1:end);
else
    % frame end is too short, so delete this frametime
    % this could be due to a later start (happens for different settings).
    Frame_End = Frame_End(2:end);
end

% amplitude threshold
max_pred = prctile(-tracet, 99.8)*0.5;

% find frame inits delete peaks that are shorter than expected frame interval
[~, Frame_Init] = findpeaks(-tracet, 'MinPeakHeight', max_pred);
Frame_Init(diff(Frame_Init) < min_frame_time*0.8) = [];
Frame_Init = Frame_Init';

% if first time init is missed, add one
if Frame_Init(1) > min_frame_time*0.8
    fprintf(['missing first frame start time: ', num2str(Frame_Init(1)), ' '])
	Frame_Init = [Frame_Init(1) - round(mean(Frame_Init(2:2:10^2) ...
        - Frame_Init(1:2:10^2))), Frame_Init];
end

if Frame_Init(1) <= 0; Frame_Init(1) = 1; end

Frame_Init = Frame_Init(1:numel(Frame_End));

fprintf([' FrameTimePoints ( ', num2str(numel(Frame_Init)), ' )'])

if frame_num == numel(Frame_Init)
    fprintf('\n');
else
    fprintf([' # of frames ~= from image (', ...
        num2str(frame_num), ')\n']);
end

% get the mode
mode_frame_width = mode(Frame_End - Frame_Init)*0.9;

% recalculate frame times

tracet = -diff(diff(smooth(mirror_y_trace(1, :), 10)));
% tracet_detrend = prctfilt(mirror_y_trace, 10, 10^3, 10, 0);
% think more about which detrending would be better

% amplitude threshold
max_pred = prctile(tracet, 99.8)*0.5;

% find frame ends & delete peaks that are closer than expected frame interval
[~, Frame_End] = findpeaks(tracet, 'MinPeakHeight', max_pred, ...
    'MinPeakDistance', mode_frame_width);
Frame_End(diff(Frame_End) < min_frame_time*0.8) = [];
Frame_End = Frame_End';

% correct for artifact peak at the beginning
if Frame_End(1) > min_frame_time*0.8
    Frame_End = Frame_End(1:end);
else
    Frame_End = Frame_End(2:end);
end

% amplitude threshold
max_pred = prctile(-tracet, 99.8)*0.5;

% find frame inits delete peaks that are closer than expected frame interval
[~, Frame_Init] = findpeaks(-tracet, 'MinPeakHeight', max_pred, ...
    'MinPeakDistance', mode_frame_width);
Frame_Init(diff(Frame_Init) < min_frame_time*0.8) = [];
Frame_Init = Frame_Init';

% manual shift of onset/offset
if ~contains(stim_file2load, 'fict')
    Frame_Init = Frame_Init + 15;
else
    Frame_End = Frame_End - 2;
    Frame_Init = Frame_Init + 4;
end

if Frame_Init(1) > min_frame_time*0.8
	Frame_Init = [Frame_Init(1) - round(mean(Frame_Init(2:2:10^2) ...
        - Frame_Init(1:2:10^2))), Frame_Init];
end

if Frame_Init(1) <= 0; Frame_Init(1) = 1; end

Frame_Init = Frame_Init(1:numel(Frame_End));

fprintf([' FrameTimePoints ( ', num2str(numel(Frame_Init)), ' )'])

if frame_num == numel(Frame_Init)
    fprintf('\n');
else
    fprintf([' # of frames ~= from image (', ...
        num2str(frame_num), ')\n']);
end

end

function stimuli_onset_offset = ...
    get_ttl_start_end(volt_trace, ...
    volt_thrs, time_init, minwidth, ...
    sampling_rate)
% get_ttl_start_end: get start and end of ttl pulses
%
% Usage:
%   [Frame_Init, Frame_End] = ...
%       colecttimestampNew(mirror_y_trace, frame_num, min_frame_time)
%
% Args:
%   volt_trace: input voltage trace
%   volt_thrs: voltage threshold
%       (4, default)
%   time_init: time to ignore
%       (100, default)
%   minwidth: minimun width of stimuli timepoints
%       (20, default)
%   sampling_rate: sampling rate
%       (10^4, default)

% voltage threshold
if ~exist('volt_thrs', 'var') || isempty(volt_thrs)
    volt_thrs = 4;
end

% time to ignore
if ~exist('time_init', 'var') || isempty(time_init)
    time_init = 100;
end

% time to ignore
if ~exist('minwidth', 'var') || isempty(minwidth)
    minwidth = 20;
end

if size(volt_trace, 1) > 1
    volt_trace = volt_trace';
end

% binarize vector:
mean_signal_start = ...
    mean(volt_trace(time_init:0.6*sampling_rate));
% sd_signal_start = ...
%     std(volt_trace(time_init:0.6*sampling_rate));

% remove artifact at the beginning of recording
volt_trace(1:time_init) = mean_signal_start;

% scale volt_trace
% volt_trace = (volt_trace - mean_signal_start)...
%     /sd_signal_start;

volt_trace = (volt_trace - mean_signal_start);

stimuli_onset_offset = ...
    find_stim_int(volt_trace, minwidth, volt_thrs);

end

function Y = local_binread(stimrel_var, cNum)
% local_binread: read '.bin' files
%
% Usage:
%   data = local_binread(stimrel_var, cNum)
%
% Args:
%   stimrel_var: stimuli related variable generated by LEDcontroler
%   cNum: number of channels to load
%
% Returns:
%   Y: vector or matrix contained in bin file

fID = fopen([strrep(stimrel_var.fName, '.mat', '') '.bin'], 'r');

if ~exist('cNum', 'var')
    
    Y = fread(fID, [length(stimrel_var.channels) inf], 'double');
    
else
    
    status = fseek(fID, 8*(cNum-1), -1);
    
    if status ~=0
        error('Could not seek for channel');
    end
    
    Y = fread(fID, 'double', (length(stimrel_var.channels)-1)*8)';
    
end

fclose(fID);

end
