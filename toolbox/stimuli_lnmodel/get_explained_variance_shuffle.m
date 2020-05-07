function [eVar_cs, corrcoef_cs] = ...
    get_explained_variance_shuffle(...
    Y, Y_pred, lgate, stim_bin, redo, ...
    filename, chunk_minsize, chunk_splitn, ...
    chunksiz, corenum, btn, shuffle2use, ...
    add_stim_lag, add_inverse)
% runridgeshuffle_chunks: Run Ridge regression to all 
%   the possible combinations of train and test data
%
% Usage:
%   eVar_cs = get_explained_variance_shuffle(...
%   	Y, Y_pred, lgate, stim_bin, redo, ...
%   	filename, chunk_minsize, chunk_splitn, ...
%   	chunksiz, corenum, btn, shuffle2use, ...
%       add_stim_lag, add_inverse)
%
% Args:
%   Y: time traces of variables [n, T]
%   Y_pred: time traces of predicted variables [n, T]
%   lgate: flag to generate permuted data
%   stim_bin: binary vector of stimuli
%   redo: redo flag
%   filename: name of temporary file to generate
%   chunk_minsize: minimun size of time chunk
%   chunk_splitn: number of times to split catrace_in in the time domain
%   chunksiz: number of chunks for parpool
%   corenum: number of cores
%   btn: number of permutations
%   shuffle2use: type of shuffling to use
%       (default, 0, random chunks, requires chunk_minsize & chunk_splitn & btn);
%       (1, random chunks + circular permutation);
%       (2, circular shuffling);
%   add_stim_lag: aditional lags of stimuli pattern 
%       to avoid (in timestamps) when doing circular permutation
%       (default, [-10 10]);
%   add_inverse: include reverse order when doing circular permutation
%       (default, 1);
%
% Outputs:
%   eVar_cs: explained variance per permutation for each ROI
%   corrcoef_cs: correlation coefficient per permutation for each ROI
%
% See also getsMod_ridgeperiter_ashuffle

if ~exist('stim_bin', 'var')
    stim_bin = [];
end

if ~exist('shuffle2use', 'var')
    shuffle2use = 0;
end

if ~exist('add_stim_lag', 'var')
    add_stim_lag = [-10 10];
end

if ~exist('filename', 'var') || isempty(filename)
    filename = [pwd, 'temp_file.mat'];
end

if ~exist('add_inverse', 'var') || isempty(add_inverse)
    add_inverse = 1;
end

siz1 = size(Y, 1);
dataObj = matfile(filename, 'Writable', true);

if lgate || redo
    % initialize parameters
    dataObj.eVar_cs = [];

    % generate shuffle in chunks over time
    if shuffle2use ~= 2
        [~, rperm_chunkshufIdx] = randchunkper(...
        Y(1, :), chunk_splitn, ...
        round(chunk_minsize), btn, stim_bin);
    end
    
    % add circular permutations
    if shuffle2use ~= 0
        rperm_circIdx = randcirshuffleper(...
            numel(Y(1, :)), 10, [], stim_bin, ...
            add_stim_lag, add_inverse);
    end

    if shuffle2use == 0
        rperm_chunkIdx = rperm_chunkshufIdx;
    elseif shuffle2use == 2
        rperm_chunkIdx = rperm_circIdx;
    else
        rperm_chunkIdx = [rperm_chunkshufIdx; rperm_circIdx];
    end
    
    dataObj.rperm_chunkIdx = rperm_chunkIdx;

    % make chunks to run
    [~, ~, chunk_idx] = ppool_makechunks(...
        chunksiz, corenum, siz1);
else
    % update initial index
    vect_init = size(dataObj.eVar_cs, 1) + 1;

    % Get shuffle time
    rperm_chunkIdx = dataObj.rperm_chunkIdx;

    % make chunks to run
    [~, ~, chunk_idx] = ppool_makechunks(...
        chunksiz, corenum, siz1, vect_init);       
end

btn = size(rperm_chunkIdx, 1);

for i = 1:numel(chunk_idx)
    
    if i == 1
        t0 = stic;
    end
    
    batch2run = chunk_idx{i};
    
    parfor ii = 1:numel(batch2run)
        
        idx_i{ii, 1} = (batch2run(ii):min(batch2run(ii)...
            + chunksiz - 1, siz1));
        k = 1;
        t_eVar = [];
        t_corrcoef = [];

        for iii = idx_i{ii, 1}
            
            if i == 1 && ii == 1 && k == 1
                t0_ = stic;
            end

            [t_eVar(k, :), t_corrcoef(k, :)] = ...
                get_explained_variance_per_row(...
                Y(iii, :), Y_pred(iii, :), rperm_chunkIdx);

            if i == 1 && ii == 1 && k == 1
                fprintf(['time per roi ', ...
                    num2str(stoc(t0_)), ' seconds\n']); 
                fprintf(['Estimated time ', ...
                    num2str(numel(chunk_idx)*chunksiz*stoc(t0_)/3600), ...
                    ' hours\n']); 
            end
        
            k = k + 1;
       end
        
        tb_eVar_cs_i{ii, 1} = t_eVar;
        tb_corrcoef_cs_i{ii, 1} = t_corrcoef;

    end
    
    dataObj.eVar_cs(cat(2, idx_i{:}), 1:btn) = ...
        cat(1, tb_eVar_cs_i{:});
    dataObj.corrcoef_cs(cat(2, idx_i{:}), 1:btn) = ...
        cat(1, tb_corrcoef_cs_i{:});
    
    clear tb_eVar_cs_i tb_corrcoef_cs_i idx_i
    
    if i == 1
        fprintf(['Estimated time ', ...
            num2str(numel(chunk_idx)*stoc(t0)/3600), ' hours\n']);
    end
    
    if mod(i, 1) == 0
        fprintf('%2.1f%% of chunks completed \n', ...
            i*100/numel(chunk_idx));
    end
    
end

eVar_cs = dataObj.eVar_cs;
corrcoef_cs = dataObj.corrcoef_cs;

end

function [eVar, corrcoef] = get_explained_variance_per_row(...
    Y, Y_pred, rperm_tIdx)
% get_explained_variance_per_row: get explained variance from many shuffle
%   iterations (as defined in rperm_tIdx)
%
% Usage:
%   eVar = get_explained_variance_per_row(...
%      Y, Y_pred, rperm_tIdx)
%
% Args:
%   Y: time traces [1, T]
%   Y_pred: predicted time traces [1, T]
%   rperm_tIdx: set of permutations to apply to Y [n, T],
%       n: permutations
%
% Outputs:
%   eVar: explained variance per permutation [n, T]
%   corrcoef: correlation coefficient [n, T]

T = size(Y, 2);

if size(rperm_tIdx, 1) == T
    rperm_tIdx = rperm_tIdx';
end

% apply shuffle to Y
Y_shuffle = zscorebigmem(Y(rperm_tIdx));
Y_pred = zscorebigmem(Y_pred);

% estimate explained variance per permutation
corrcoef = corr(Y_shuffle', Y_pred');
eVar = (corrcoef').^2;

end
