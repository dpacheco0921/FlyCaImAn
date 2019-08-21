function batch_LEDdenoise(FolderName, FileName, iparams)
% batch_LEDdenoise: Runs LED denoising in all folders and files within it
%
% Usage:
%   batch_LEDdenoise(FolderName, FileName, iparams)
%
% Args:
%   FolderName: Folder name to load
%   FileName: File name to load
%   iparams: parameters to update
%       (pgate: plot gate)
%           (default, 0)
%       (redo: redo gate)
%           (default, 0)
%       (debug: debug gate)
%           (default, 0)
%       (minput: min fluorescence input)
%           (default, 0)
%       (iter: number of denoising interations)
%           (default, 1)
%       (baseline: timepoints to use as baseline)
%           (default, [])
%       (delredcha: gate to delete red channel)
%           (default, 1)
% 
% Notes:
% it generates iDat.LED and iDat.PMT_fscore

% Default params
ledpars = [];
ledpars.pgate = 0;
ledpars.redo = 0;
ledpars.debug = 0;
ledpars.minput = 0;
ledpars.iter = 1; 
ledpars.baseline = [];
ledpars.delredcha = 1;

if ~exist('FolderName', 'var')|| isempty(FolderName)
    FolderName = [];
end
if ~exist('FileName', 'var') || isempty(FileName)
    FileName = [];
end
if ~exist('iparams', 'var'); iparams = []; end
ledpars = loparam_updater(ledpars, iparams);

% Selecting folders
cDir = pwd;
f2reject = {'.', '..', 'preprocessed', 'BData'};
f2run = dir;
f2run = str2match(FolderName, f2run);
f2run = str2rm(f2reject, f2run);
f2run = {f2run.name};

fprintf(['Running n-folders : ', num2str(numel(f2run)), '\n'])

for i = 1:numel(f2run)
    
    fprintf(['Running folder : ', f2run{i}, '\n']); 
    cd(f2run{i});
    runperfolder(FileName, ledpars);
    cd(cDir)
    
end

fprintf('... Done\n')

end

function runperfolder(FileName, ledpars)
% runperfolder: Getting rid of LED bleed through green channel,
% generates cDat which stores all LED correction parameters and results
%
% Usage:
%   runperfolder(FileName, ledpars)
%
% Args:
%   FileName: image related field
%   ledpars: internal parameters

% Denoising green channel

global GreenCha RedCha Data

if ~exist('FileName','var') || isempty(FileName)
    FileName = [];
end

% Files to load
sep = filesep;
File2Run = rdir(['.', sep, '*_rawdata.mat']);

if ~isempty(FileName)
    if iscell(FileName)
        for i = 1:numel(FileName)
            FileName{i} = [FileName{i}, '_rawdata.mat'];
        end
    else
        FileName = [FileName, '_rawdata.mat'];
    end
end

File2Run = str2match(FileName, File2Run);
File2Run = {File2Run.name};
fprintf('Correcting LED bleed to green channel\n')
qgate = 0;

for F2R_idx = 1:numel(File2Run)
    if qgate == 0
        
        % Default ledpars, cDat stores informative correction ledpars
        Data = [];
        GreenCha = [];
        RedCha = [];
        cDat.CorType = 'LTM';
        cDat.buffer = 500;
        
        % Loading files
        fprintf(['Correcting file: ', ...
            File2Run{F2R_idx}((max(strfind(File2Run{F2R_idx}, sep))+1):end), ' '])
        
        % Loading metadata
        load(strrep(File2Run{F2R_idx}, '_rawdata', '_metadata'), ...
            'iDat', 'lStim', 'fDat')
        load(strrep(File2Run{F2R_idx}, '_rawdata', ''), ...
            'sDat')
        
        % Loading Data (uint16) and keeping the format as double
        load(File2Run{F2R_idx})
        dDim = size(Data);
        Data = double(Data);
        
        if iDat.LEDCorr == 0 || ledpars.redo == 1
            
            if ledpars.redo == 1 && isfield(lStim, 'lStim2D')
               fprintf(' (reseting lStim2D & cDat) ');
               lStim = rmfield(lStim, 'lStim2D');
            end
            
            % Separating green / red channel
            ChannelSplitter(dDim)
            
            % Identifying baseline stacks / frames
            if isempty(ledpars.baseline) 
                ledpars.baseline = baselinegenerator(iDat, sDat);
            else
                ledpars.baseline
            end
            
            % Correcting green channel
            if ledpars.pgate
                axH = figuregenerator(ledpars.iter);
            else
                axH = zeros(1, ledpars.iter*3);
            end
            
            for iter_i = 1:ledpars.iter
                
                fprintf(['Running iteration # ', num2str(iter_i), '\n']);
                
                cDat.Iter = iter_i;
                
                if ledpars.iter == 1
                    AxIDx = iter_i:ledpars.iter:ledpars.iter*3;
                else
                    AxIDx = iter_i:3:ledpars.iter*3;
                end
                
                % it prunes pixels (reduces the heigth by 2 pixels)
                [GreenCha, cDat, lStim] = LEDdenoiser(GreenCha, ledpars.baseline, ...
                    iDat, lStim, axH(AxIDx), cDat, ledpars.pgate);
                
            end
            
            % Delete lStim.lStim2D field
            lStim = rmfield(lStim, 'lStim2D');
            
            % Do the same for RedCha (reduces the heigth by 2 pixels)
            RedCha = RedCha(1:end-2, :, :, :);
            lStim = rmfield(lStim, 'pstEn');
            
            % Calculate mean image of green and red channel
            [iDat, cDat.redchaths, hDat, qgate] = ...
                getMeanPerCha(iDat, sDat, fDat, ledpars);
            if qgate; return; end
            
            % Update fram size
            iDat.FrameSize(1) = size(GreenCha, 1);
            
            % Saving Data
            ChannelMerger(dDim);
            
            % Delete frame during flyback / update frameN / RedChaMean
            iDat = deleteflybackframe(iDat, fDat);
            iDat.LEDCorr = 1;
            
            % saving processed video to .mat file
            save(strrep(File2Run{F2R_idx}, '_rawdata', '_metadata'), ...
                'iDat', 'lStim', 'cDat', 'hDat', '-append')
            save(File2Run{F2R_idx}, 'Data', '-v7.3')
            
            GreenCha = [];
            RedCha = [];
            Data = [];
            
        else
            
            fprintf(' *already corrected*\n')
            
        end
        
        clear dDim iDat lStim sDat fDat cDat
        
    else
        
        fprintf('Stopped some files need manual editing\n')
        
    end
    
end

fprintf('****** Done ******\n')

end

function iDat = deleteflybackframe(iDat, fDat)
% deleteflybackframe: deleting z frame and related time stamps (Piezzo fly back)
%
% Usage:
%   iDat = deleteflybackframe(iDat, fDat)
%
% Args:
%   iDat: image related field
%   fDat: file related field

global Data

% choose plane to delete based on data input type
if ~isempty(strfind(fDat.DataType, '3DxT'))
    
    if ~isempty(strfind(fDat.DataType, 'old'))
        
        % old, Chop Data and update iDat
        z2del = size(Data, 3);
        iDat = rmfield(iDat, 'sstEn');
        
    else
        
        % new, Chop Data and update iDat
        z2del = 1;
        iDat = rmfield(iDat, 'sstEn');
        
    end
    
    % update frame time / Data / RedChaMean
    fprintf('Deleting flyback frame, final size: ')
    iDat.fstEn(z2del:iDat.FrameN:end, :) = [];
    
    Data(:, :, z2del, :, :) = [];
    
    if ~isempty(iDat.RedChaMean)
        iDat.RedChaMean(:, :, z2del) = [];
    end
    
    if ~isempty(iDat.GreenChaMean)
        iDat.GreenChaMean(:, :, z2del) = [];
    end
    
else
    
    fprintf('It is 2DxT data, so it does not require slice removal ')
    
end

if ~isempty(strfind(fDat.DataType, '3DxT'))
    iDat.FrameN = size(Data, 3);
elseif ~isempty(strfind(fDat.DataType, '2DxT'))
    iDat.FrameN = 1;
end

fprintf([num2str(iDat.FrameN), '\n'])

end

function [iDat, redchaths, hDat, qgate] = ...
    getMeanPerCha(iDat, sDat, fDat, ledpars)
% getMeanPerCha: Calculates mean image of green and red channel
%
% Usage:
%   [iDat, redchaths, hDat, qgate] = ...
%   getMeanPerCha(iDat, sDat, fDat, ledpars)
%
% Args:
%   iDat: image related field
%   sDat: stimuli related field
%   fDat: file related field
%   ledpars: stimuli related field
%
% Notes:
% For the red channel, it estimates a intensity threshold based on
% distribution of max F values per frame to tell apart volumes where the
% PMT was off or not.

global RedCha GreenCha
qgate = 0;
redchaths = [];
hDat = [];

% Get mean of GreenChannel
if ~isempty(strfind(sDat.stimFileName, 'Light'))
    if ~isempty(RedCha)
        
        iDat.PMT_fscore = squeeze(max(max(RedCha, [], 1), [], 2));
        
        if size(iDat.PMT_fscore, 1) > size(iDat.PMT_fscore, 2)
            % it assumes that the largest dimention should be time.
            iDat.PMT_fscore = iDat.PMT_fscore';
        end
        
        if sum(iDat.PMT_fscore(:) > 0) ~= 0
            
            % estimate redchaths
            FitType = 'gauss2';
            sigmaRed = 5;
            hDat = fitgauss1D(iDat.PMT_fscore(:), FitType, sigmaRed, 1);
            
            % get volumes / frames greater than ths
            if hDat.coef(5) - hDat.coef(2) >= 10 && ledpars.minput == 0
                
                if ledpars.pgate
                    fitgaussPlotter(hDat);
                end
                
                redchaths = hDat.xths(1, 2);
                
                if size(iDat.PMT_fscore, 1) == 1
                    iDat.PMT_fscore = iDat.PMT_fscore > redchaths;
                    iDat.RedChaMean = mean(RedCha(:, :, iDat.PMT_fscore == 1), 3);
                else
                    iDat.PMT_fscore = sum(iDat.PMT_fscore > redchaths, 1);
                    iDat.RedChaMean = mean(RedCha(:, :, :, iDat.PMT_fscore == iDat.FrameN), 4);
                end
                
            elseif ledpars.minput == 1
                
                % load Data from matfile
                load([fDat.FileName,'_redcha.mat'], 'Data')
                
                % delete last 2 pixels
                iDat.RedChaMean = double(Data(1:end-2, :, :));
                
            else
                
                if ledpars.debug == 1
                    
                    edit LEDDebug;
                    fitgaussPlotter(hDat);
                    keyboard;
                    
                else
                    
                    fprintf('File needs manual editing');
                    qgate = 1;
                    
                end
                
            end
        else
            
            iDat.RedChaMean = [];
            
        end
        
    else
        
        iDat.RedChaMean = [];
        
    end
    
    if ledpars.delredcha; RedCha = []; end
    
end

% Get mean of GreenChannel
if iDat.FrameN > 1
    iDat.GreenChaMean = mean(GreenCha, 4);
else
    iDat.GreenChaMean = mean(GreenCha, 3);
end

end

function ChannelSplitter(dDim)
% ChannelSplitter: splits channels to generate GreenCha RedCha
%
% Usage:
%   ChannelSplitter(dDim)
%
% Args:
%   dDim: original dimensions

global GreenCha RedCha Data

if length(dDim) == 4 && dDim(end) == 2
    
    % 2DxTxCh
    RedCha = Data(:, :, :, 1);
    GreenCha = Data(:, :, :, 2);
    
elseif length(dDim) == 5 && dDim(end) == 2
    
    % 3DxTxCh
    RedCha = Data(:, :, :, :, 1);
    GreenCha = Data(:, :, :, :, 2);
    
else
    
    % 2DxT
    GreenCha = Data;
    
end

Data = [];

end

function ChannelMerger(dDim)
% ChannelMerger: merges channels to generate Data
%
% Usage:
%   ChannelMerger(dDim)
%
% Args:
%   dDim: original dimensions

global GreenCha RedCha Data

if isempty(RedCha) 

    % Reducing Data to just green channel
    Data = GreenCha;
    
else    

    % Reconstituting Data from red and 
    %   corrected-green channel (1, 2, respectively) 
    Data = zeros([size(GreenCha), 2]);
    
    if length(dDim) == 4
        Data(:, :, :, 1) = RedCha; Data(:, :, :, 2) = GreenCha;
    elseif length(dDim) == 5
        Data(:, :, :, :, 1) = RedCha; Data(:, :, :, :, 2) = GreenCha;
    else
        Data = GreenCha;
    end
    
end

GreenCha = [];
RedCha = [];
Data = double(Data);

end

function Baseline = baselinegenerator(iDat, sDat)
% baselinegenerator: Generates baseline vector, depending on protocol details
%
% Usage:
%   Baseline = baselinegenerator(iDat, sDat)
%
% Args:
%   iDat: image related field
%   sDat: stimuli related field

baselineTime = (sDat.silencePre/10^3)*sDat.fs; % from ms to timepoints

if isfield(iDat, 'sstEn')
    % if this is a volume & low temporal resolution (0.5-1 Hz)
    Int = (iDat.sstEn(1, 2) - iDat.sstEn(1, 1));
    Baseline = 4:(ceil(baselineTime/Int)-2);
else
    % if this is a plane & high temporal resolution (8-10 Hz)
    Int = (iDat.fstEn(1, 2) - iDat.fstEn(1, 1));
    Baseline = 10:(ceil(baselineTime/Int)-10);
end

fprintf([' baseline (', num2str(Baseline(1)), ...
    ':', num2str(Baseline(end)), ') ']);

end

function axH = figuregenerator(Iterations)
% figuregenerator: generate axis for plotting results
%
% Usage:
%   axH = figuregenerator(Iterations)
%
% Args:
%   Iterations: number of iteration for LED denoising

if Iterations == 1
    figure('position', [248 752 1440 331])
elseif Iterations == 2
    figure('position', [248 460 1440 623])
elseif Iterations == 3
    figure('position', [248 138 1440 945])
end

k = 1;
for i = 1:Iterations
    for ii = 1:3
    	eval(['axH(', num2str(i),', ', num2str(ii), ...
            ') = subplot(', num2str(Iterations), ', 3, ', ...
            num2str(k), ');'])
        k = k + 1;
    end
end

end
