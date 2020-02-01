function batch_plot_df_dfof_maxf(Filename, oDir, iparams)
% batch_plot_df_dfof_maxf: Generate
%   raw data videos (DF, DF\Fo, and MaxF) for a designated directory
%
% Usage:
%   batch_plot_df_dfof_maxf(Filename, oDir, iparams)
%
% Args:
%   Filename: name pattern of files to use
%       (default, [])
%   oDir: output directory.
%   iparams: parameters to update
%       (metadata_suffix: suffix of files to load, metadata files)
%       (rawdata_suffix: suffix of files to load, raw data files)
%       %%%%%%%%%%%% video settings %%%%%%%%%%%%
%       (vgate: save gate)
%       (vquality: video quality)
%       (frate: frame rate)
%       (cmap: foreground image colormap)
%       %%%%%%%%%%%% DF calculation %%%%%%%%%%%%
%       (baseline_tp: baseline timepoints)
%       (sign2use: use positive or negative changes (by multipliying it by 1/-1))
%           (default, 0 (absolute))
%           (-1 (negative))
%           (0 (absolute))
%       (df_flag: flag to calculate DF (if 1), or DF/Fo otherwise)
%           (default, [0 1 2] (get all))
%           (0 (get DF/Fo))
%           (1 (get DF))
%           (2 (get Max F))
%       (chunk_size: size of chunks)
%       %%%%%%%%%%%% parpool & server related %%%%%%%%%%%%
%       (serId: server id)
%           (default, 'int')
%       (corenum: number of cores)
%           (default, 4)
%       %%%%%%%%%%%% extra %%%%%%%%%%%%
%       (dir_depth: depth of directory search)
%           (default, 0)
%       (vstr: input stimuli name)
%           (default, [])
%
% Notes:
%   See get_df_MIP

% parameters
ipars.metadata_suffix = '_metadata.mat';
ipars.rawdata_suffix = '_rawdata.mat';
ipars.baseline_tp = 1:50;
ipars.vfontsiz = 10;
ipars.vgate = 1;
ipars.vquality = 100;
ipars.frate = 10;
ipars.cmap = parula;
ipars.sign2use = 0;
ipars.plottype = [0 1 2];
ipars.chunk_size = 10;
ipars.serId = [];
ipars.corenumber = [];
ipars.dir_depth = 0;
ipars.vstr = [];
ipars.axisratio = 1;

if ~exist('Filename', 'var')
    Filename = [];
end

if ~exist('oDir', 'var') || isempty(oDir)
    oDir = [pwd, filesep, 'dfrel_vid'];
end

if ~exist('iparams', 'var'); iparams = []; end
ipars = loparam_updater(ipars, iparams);

if ~isempty(oDir) && exist('oDir', 'var')
    mkdir(oDir)
end

% define parpool settings
if isempty(ipars.corenumber)
    if ispc
        ipars.corenumber = 4;
    else
        ipars.corenumber = 6;
    end
end

if isempty(ipars.serId)
    if ispc
        ipars.serId = 'int';
    else
        ipars.serId = 'spock';    
    end
end

% define files to use
if ipars.dir_depth == 0
    f2run = rdir(['.', filesep, '*', ...
        ipars.metadata_suffix]);
elseif ipars.dir_depth == 1
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', ipars.metadata_suffix]);
else
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', filesep, '*', ...
        ipars.metadata_suffix]);
end

f2run = {f2run.name}';
f2run = str2rm({'Zstack'}, f2run);
[filename, iDir] = split_path(f2run);
filename = strrep(filename, ipars.metadata_suffix, '');

if ~isempty(Filename)
    f2run = find(contains(filename, Filename));
    filename = filename(f2run);
    iDir = iDir(f2run);
end

fprintf(['Generating videos for ', ...
    num2str(numel(filename)), ' files\n'])

for i = 1:numel(filename)
    
        tic
        fprintf(['plotting ', filename{i}, '\n'])
        
        load([iDir{i}, filesep, filename{i}, ...
            ipars.metadata_suffix], 'wDat')
        time2load = 1:numel(wDat.fTime);
        
        MIP_proj = get_df_MIP(...
            [iDir{i}, filesep, filename{i}, ipars.rawdata_suffix], ...
            [iDir{i}, filesep, filename{i}, ipars.metadata_suffix], ...
            ipars.baseline_tp, [], time2load, ...
            ipars.sign2use, ipars.plottype, ...
            ipars.chunk_size, ipars.corenumber, ...
            ipars.serId);

        % get stim name
        if isempty(ipars.vstr)
            stimuli_name = cellfun(@(x, y) [x '_AMP_' num2str(y)], ...
                wDat.sPars.name(1, :), chunk2cell(wDat.sPars.int(1, :), 1), ...
                'UniformOutput', false);
            [stimuli_name, ~, ~] = unique(stimuli_name, 'stable');
            ipars.vstr = stimuli_name;
        end

        % get stim indeces
        ipars.vstrt = zeros(size(getStimVect(wDat, 1, [])));

        for j = 1:numel(ipars.vstr)
            ipars.vstrt = ipars.vstrt + getStimVect(wDat, 1, [], j)*j;
        end
        ipars.vstrt = ipars.vstrt(time2load);

        ipars.vname = [oDir, filesep, filename{i}, '_MIP_DF'];
        ipars.range = [0 1];
        
        if sum(ismember(ipars.plottype, 0))
            fprintf('plot DF \n')
            slice3Dmatrix(flip(MIP_proj{1}, 2), ipars)
        end
        
        ipars.range = [0 1];
        ipars.vname = [oDir, filesep, filename{i}, '_MIP_DFoF'];

        if sum(ismember(ipars.plottype, 1))
            fprintf('plot DFoF \n')
            slice3Dmatrix(flip(MIP_proj{2}, 2), ipars)
        end
        
        ipars.range = [0 1];
        ipars.vname = [oDir, filesep, filename{i}, '_MIP_maxF'];

        if sum(ismember(ipars.plottype, 2))
            fprintf('plot maxF \n')
            slice3Dmatrix(flip(MIP_proj{3}, 2), ipars)
        end
        
        clear wDat
        toc

end

end
