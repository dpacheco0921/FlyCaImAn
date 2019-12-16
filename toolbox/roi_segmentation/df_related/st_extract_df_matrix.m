function st_extract_df_matrix(obj, f2use, bs_range, ...
    sti_range, updgate, redo, stim2use)
% st_extract_df_matrix: generates a GreenChaDf and GreenChaDfof signal
% relative to baseline prior to selected stimuli
% Usage:
%   st_clus_perfile(obj, f2use, bs_range, ...
%       sti_range, updgate, redo, stim2use)
%
% Args:
%   obj: SpaTemp object
%   f2use: file to use
%   bs_range: time before stim and after stim to use for baseline and df
%       calculation in seconds (s)
%   sti_range: time before stim and after stim to use for sti in seconds (s)
%   updgate: gate to update prosmetadata file
%   redo: flag to redo
%   stim2use: stimuli to use
%
% Notes:
% this function requires obj fields:
%   (obj.wDat.fTime)
%   (obj.dirrel.cDir)
%   (obj.dirrel.foname)
%   (obj.dirrel.finame)
%   (obj.wDat.sVector)

if ~exist('f2use', 'var') || isempty(f2use)
    f2use = 1:numel(obj.dirrel.finame);
end
if ~exist('bs_range', 'var') || isempty(bs_range); bs_range = [-10 -0.5]; end
if ~exist('sti_range', 'var') || isempty(sti_range); sti_range = [0 5]; end
if ~exist('updgate', 'var') || isempty(updgate); updgate = 0; end
if ~exist('redo', 'var') || isempty(redo); redo = 0; end
if ~exist('stim2use', 'var') || isempty(stim2use); stim2use = []; end
if ~iscell(stim2use)
    stim2use = {stim2use};
end

if size(f2use, 1) > 1
    f2use = f2use';
end

dt = diff(obj.wDat.fTime{1}(1:2));
bs_range_f = bs_range/dt;
sti_range_f = sti_range/dt;

% get file's directory
fipath = [obj.dirrel.cDir, filesep, obj.dirrel.foname];

% generate wDat fields
if ~isfield(obj.wDat, 'GreenChaDf')
   obj.wDat.GreenChaDf = cell(numel(obj.dirrel.finame), ...
       size(sti_range, 1), numel(stim2use));
   obj.wDat.GreenChaDfof = cell(numel(obj.dirrel.finame), ...
       size(sti_range, 1), numel(stim2use));
end

for f_i = f2use
    if isempty(obj.wDat.GreenChaDf{f_i}) || redo
        fprintf(['Generating Df image from # ', ...
            obj.dirrel.finame{f_i}, ' files\n'])
        tic

        dataObj = matfile([fipath, filesep, ...
            obj.dirrel.finame{f_i}, '_prosdata.mat'], 'writable', false);
        
        for sti_i = 1:numel(stim2use)
            
            % get start and end of selected stimuli
            sti_En = getStim_InitEnd(obj.wDat.fTime{f_i}, ...
                obj.wDat.sTime{f_i}(ismember(obj.wDat.sIdx{f_i}, ...
                stim2use{sti_i}), :));

            siz = obj.wDat.bSiz{f_i};

            % iterate for each stimuli range

            for r_i = 1:size(sti_range, 1)

                fprintf('*')

                % initialize matrix size
                df_fo_t = single(zeros([siz(1:3), size(sti_En, 1), 2]));
                df_t_ = single(zeros([siz(1:3), size(sti_En, 1), 2]));

                for t_i = 1:size(sti_En, 1)

                    % load baseline fluorescence 'bf'
                    fo = dataObj.Y(:, :, :, (sti_En(t_i, 1) + bs_range_f(1)): ...
                        sti_En(t_i, 1) + bs_range_f(2));
                    fo = mean(fo, 4);
                    fo = imblur(fo, [2 2 2], [3 3 3], 3);

                    % load fluorescence during stimuli 'stif'
                    df_t = dataObj.Y(:, :, :, (sti_En(t_i, 1) + sti_range_f(r_i, 1)):  ...
                        (sti_En(t_i, 2) + sti_range_f(r_i, 2)));

                    % smooth (necessary otherwise it is pretty noisy per voxel)
                    df_t = imblur(df_t, [2 2 2], [3 3 3], 3);

                    % get 80th or 20th percentile
                    df_t = df_t - fo;
                    df_fo_h = prctile(df_t, 80, 4);
                    df_fo_l = prctile(df_t, 20, 4);

                    % concatenate both
                    df_t_(:, :, :, t_i, :) = cat(5, df_fo_h, df_fo_l);
                    df_fo_t(:, :, :, t_i, :) = df_t_(:, :, :, t_i, :)./abs(fo);

                end

                % get mean df or df/fo across trials
                obj.wDat.GreenChaDf{f_i, r_i, sti_i} = squeeze(mean(df_t_, 4));
                obj.wDat.GreenChaDfof{f_i, r_i, sti_i} = squeeze(mean(df_fo_t, 4));

                df_fo_t = []; df_t_ = [];

            end
            
            sti_En = []; siz = [];

        end
                
        toc
        fprintf('\n')
        
    end
end

% update wDat on matfile
for f_i = f2use
    
    if ~isempty(obj.wDat.GreenChaDf{f_i}) && updgate
        
        % load wDat
        load([fipath, filesep, obj.dirrel.finame{f_i}, '_prosmetadata.mat'], 'wDat')
        
        % update wDat
        wDat.GreenChaDf = obj.wDat.GreenChaDf(f_i, :, :);
        wDat.GreenChaDfof = obj.wDat.GreenChaDfof(f_i, :, :);
        
        save([fipath, filesep, obj.dirrel.finame{f_i}, '_prosmetadata.mat'], ...
            'wDat', '-append')
        
        clear wDat
        
    end
    
end

end
