function MIP_proj = get_df_MIP(...
    rawdata_name, metadata_name, baseline_tp, ...
    slices2project, time2load, sign2use, ...
    df_flag, chunk_size, corenumber, serId)
% get_df_MIP: generate maximun intensity projection
%
% Usage:
%   MIP_proj = get_df_MIP(...
%   	rawdata_name, metadata_name, baseline_tp, ...
%   	slices2project, time2load, sign2use, ...
%       df_flag, chunk_size, corenumber)
%
% Args:
%   rawdata_name: name of '_prosdata.mat' file to use.
%   metadata_name: name of '_prosmetadata.mat' file to use.
%   baseline_tp: time points use as baseline.
%       (default, 1:9)
%   slices2project: slices to use for MIP projection, this is a cell that
%       can have different set of slices to collect.
%       (default, 1:wDat.vSize(3))
%   time2load: timepoints to load (usually 1:234).
%       (default, 1:numel(wDat.fTime))
%   sign2use: use positive or negative changes (by multipliying it by 1/-1)
%       (default, 1 (positive))
%       (-1 (negative))
%       (0 (absolute))
%   df_flag: flag to calculate DF (if 1), or DF/Fo otherwise.
%       (default, 1 (get DF))
%       (0 (get DF/Fo))
%       (2 (get Max F))
%   chunk_size: size of chunks.
%   corenumber: number of cores to use.
%   serId: server ID.
%       ('int' | ispc or ismac, to run locally)
%
% Notes:
% it requires a metadata variable wDat
% 2019-09-22: uses paralel parpool per patches

if ~exist('rawdata_name', 'var') || isempty(rawdata_name)
    rawdata_name = [];
end

if ~exist('metadata_name', 'var') || isempty(metadata_name)
    metadata_name = [];
end

if ~exist('baseline_tp', 'var') || isempty(baseline_tp)
    baseline_tp = 1:9;
end

if isempty(rawdata_name) || isempty(metadata_name)
    return
end

if ~exist('sign2use', 'var') || isempty(sign2use)
    sign2use = 1;
end

if ~exist('df_flag', 'var') || isempty(df_flag)
    df_flag = 1;
end

if ~exist('chunk_size', 'var') || isempty(chunk_size)
    chunk_size = 20;
end

if ~exist('corenumber', 'var') || isempty(corenumber)
    corenumber = 4;
end

if ~exist('serId', 'var') || isempty(serId)
    serId = 'int';
end

% set parpool
setup_parpool(serId, corenumber);

% load and mat required variables
dataObj = matfile(rawdata_name);
load(metadata_name, 'wDat')

% initialize internal variables
if ~exist('time2load', 'var') || isempty(time2load)
    time2load = 1:numel(wDat.fTime);
end

if exist('slices2project', 'var') && ~isempty(slices2project)
    z = slices2project;
else
    z = {1:wDat.vSize(3)};
end

MIP_proj = cell(numel(z), 3);

% get Y dimensions
n_d = numel(dataObj.sizY);

fprintf('load and generate baseline\n')
tic
if n_d == 4
    bas = mean(double(dataObj.Y(:, :, :, baseline_tp)), 4);
else
    bas = mean(double(dataObj.Y(:, :, baseline_tp)), 3);
end
toc

fprintf('load and generate DF, DFoFo, or MaxF\n')

tic
for i = 1:numel(z)
    
    % initialize variables
    temp_o_i = zeros([wDat.fSize numel(time2load)]);
    temp_o_ii = zeros([wDat.fSize numel(time2load)]);
    temp_o_iii = zeros([wDat.fSize numel(time2load)]);
    
    [~, ~, chunk_idx] = ...
        ppool_makechunks(chunk_size, ....
        corenumber, numel(time2load), ...
        time2load(1));
    
    % time patches    
    for j = 1:numel(chunk_idx)

        if j == 1; t0 = stic; end

        batch2run = chunk_idx{j};
        t_idx = cell(corenumber, 1);
        
        temp_i = cell(corenumber, 1);
        temp_ii = cell(corenumber, 1);
        temp_iii = cell(corenumber, 1);

        parfor ii = 1:numel(batch2run)

            t_idx{ii, 1} = (batch2run(ii):min(batch2run(ii) + ...
                chunk_size - 1, numel(time2load)))';
            iv = 1;

            for iii = t_idx{ii, 1}'

                [temp_i{ii}(:, :, iv), temp_ii{ii}(:, :, iv), ...
                    temp_iii{ii}(:, :, iv)] = ...
                    load_process_and_project_Y(dataObj, ...
                    time2load(iii), z{i}, n_d, sign2use, ...
                    df_flag, bas);
                iv = iv + 1;

            end
            
            temp_i{ii} = temp_i{ii}(:);
            temp_ii{ii} = temp_ii{ii}(:);
            temp_iii{ii} = temp_iii{ii}(:);
            
        end

        t_idx = cell2mat(t_idx);

        temp_o_i(:, :, t_idx) = reshape(cell2mat(temp_i), ...
            [wDat.fSize numel(t_idx)]);
        temp_o_ii(:, :, t_idx) = reshape(cell2mat(temp_ii), ...
            [wDat.fSize numel(t_idx)]);
        temp_o_iii(:, :, t_idx) = reshape(cell2mat(temp_iii), ...
            [wDat.fSize numel(t_idx)]);
        
        if j == 1 
            stocf(t0, 'Time per chunk ');
            fprintf(['Estimated time ', ...
                num2str(stoc(t0)*numel(chunk_idx)/60), ...
                ' minutes\n'])
        end

        if mod(j, 10) == 0
            fprintf('%2.1f%% of chunks completed \n', ...
                j*100/numel(chunk_idx));
        end

    end
    
    % collect projections
    MIP_proj{i, 1} = temp_o_i/max(temp_o_i(:));
    MIP_proj{i, 2} = temp_o_ii;
    MIP_proj{i, 3} = temp_o_iii/max(temp_o_iii(:));
    clear temp_o_i temp_o_ii temp_o_iii

end

toc

end

function [F_stat_i, F_stat_ii, F_stat_iii] = ...
    load_process_and_project_Y(dataObj, ...
    time_idx, planes_idx, n_d, sign2use, ...
    df_flag, bas)
% load_process_and_project_Y: load, process 
%   and project volumetric or plane timeseries
%
% Usage:
%   [F_stat_i, F_stat_ii, F_stat_iii] = ...
%       load_process_and_project_Y(dataObj, ...
%       time_idx, planes_idx, n_d, sign2use, ...
%       df_flag, bas)
%
% Args:
%   dataObj: data object, data is stored at Y.
%   time_idx: timepoints to load.
%   planes_idx: planes to load.
%   n_d: data dimensions.
%   sign2use: use positive or negative changes (by multipliying it by 1/-1).
%       (default, 1 (positive))
%       (-1 (negative))
%       (0 (absolute))
%   df_flag: flag to calculate DF (if 1), or DF/Fo otherwise.
%       (default, 1 (get DF))
%       (0 (get DF/Fo))
%       (2 (get Max F))
%   bas: baseline fluorescence.

if ~exist('planes_idx', 'var'); planes_idx = []; end

temp_i = [];
temp_ii = [];
temp_iii = [];

% load raw data
if n_d == 4
    temp_ = double(dataObj.Y(:, :, planes_idx, time_idx));
else
    temp_ = double(dataObj.Y(:, :, time_idx));
end

% 1) calculate DF
if sum(ismember(df_flag, 1))

    temp_i = imblur(temp_ - bas(:, :, planes_idx), 1, 3, n_d - 1);

    if sign2use == 0
        temp_i = abs(temp_i);
    else
        temp_i = temp_i*sign2use;
    end

end

% 2) calculate DF/Fo
if sum(ismember(df_flag, 0))

    temp_ii = imblur(temp_ - bas(:, :, planes_idx), 1, 3, n_d - 1);

    if sign2use == 0
        temp_ii = abs(temp_ii);
    else
        temp_ii = sign2use*temp_ii;
    end

    temp_ii = temp_ii./abs(imblur(bas(:, :, planes_idx), 1, 3, n_d - 1));

end

% 3) calculate Max
if sum(ismember(df_flag, 2))
    temp_iii = max(temp_, [], 3);
end

% collect MIP
F_stat_i = max(temp_i, [], 3);
F_stat_ii = max(temp_ii, [], 3);
F_stat_iii = max(temp_iii, [], 3);

end
