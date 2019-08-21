function wDat = getStimInfo(wDat, iDat, fDat, lStim, ...
    cDat, mcDat, fname, rep_i, reps, shifts_align, zlength)
% getStimInfo: compiles all relevant stimuli parameters into wDat structure
%
% Usage:
%   Dat = getStimInfo(wDat, iDat, fDat, lStim, ...
%       cDat, mcDat, fname, rep_i, reps, shifts_align, zlength)
%
% Args:
%   filename: file name
%   wDat: main metadata structure
%   iDat: image metadata structure
%   fDat: file name metadata structure
%   lStim: stimuli metadata structure
%   cDat: LED correction metadata structure
%   mcDat: motion correction metadata structure
%   fname: file name
%   rep_i: rep index
%   reps: total number of reps
%   shifts_align: stitch shifts
%   zlength: index of planes
%
% Notes:
% Input lStim.basPre and lStim.basPro are in ms, so it applies the right scaling

% initial unit is ms, final unit is seconds
unit_cor_factor = 10^3;

n_reps = numel(reps);

if rep_i == 1
    
    % general metadata
    if isfield(lStim, 'fStrain')
        wDat.fStrain = lStim.fStrain;
    else
        wDat.fStrain = [];
    end
    
    wDat.datatype = fDat.DataType;
    wDat.XYZres = iDat.MetaData(2:end); 
    wDat.fTime = iDat.Tres;
    wDat.Tn = length(wDat.fTime);
    wDat.sTime = iDat.lstEn;
    
    % preprocess steps
    wDat.prepros.LEDcor = iDat.LEDCorr;
    wDat.prepros.MotCor = iDat.MotCorr;
    wDat.prepros.sRes = iDat.XYresmatch;
    wDat.prepros.sSmooth = iDat.sSmooth;
    wDat.prepros.tRes = iDat.tResample;
    wDat.prepros.tSmooth = iDat.tSmooth;
    
end

% Collect song-related
if contains(fDat.DataType, 'song') || contains(fDat.DataType, 'prv')
    
    % it assumes that within a file name, the stimuli params are the same
    wDat.sTrace = lStim.trace;
    
    if isfield(lStim, 'sPars') && isfield(lStim.sPars, 'basPre')
        
        if rep_i == 1
            wDat.sPars.name = cell(n_reps, length(lStim.sPars.name));
            wDat.sPars.int = [];
            wDat.sPars.sr = [];
        end
        
        wDat.sPars.order(rep_i, :) = lStim.sPars.order;
        wDat.sPars.name(rep_i, :) = lStim.sPars.name;
        wDat.sPars.int(rep_i, :) = lStim.sPars.int';
        wDat.sPars.sr(rep_i, :) = lStim.sPars.sr';
        wDat.sPars.basPre(rep_i, :) = lStim.sPars.basPre/unit_cor_factor;
        wDat.sPars.basPost(rep_i, :) = lStim.sPars.basPost/unit_cor_factor;
        
    else
        
        load([fname, '_vDat.mat'], 'rDat');
        
        if rep_i == 1
            wDat.sPars.name = cell(n_reps, length(rDat.ctrl.stimFileName));
            wDat.sPars.int = [];
            wDat.sPars.sr = [];
        end
        
        wDat.sPars.order(rep_i, :) = rDat.stimOrder;
        wDat.sPars.name(rep_i, :) = rDat.ctrl.stimFileName;
        wDat.sPars.int(rep_i, :) = cell2mat(rDat.ctrl.intensity)';
        wDat.sPars.sr(rep_i, :) = rDat.ctrl.rate';
        wDat.sPars.basPre(rep_i, :) = rDat.ctrl.silencePre/unit_cor_factor;
        wDat.sPars.basPost(rep_i, :) = rDat.ctrl.silencePost/unit_cor_factor;
        clear rDat
        
    end
    
end

% Collect opto-related
if contains(fDat.DataType, 'opto') && ~contains(fDat.DataType, 'prv')
    
    wDat.sTrace = lStim.trace;
    
    if isfield(lStim, 'sPars') && isfield(lStim.sPars, 'basPre')
        wDat.sPars.int(rep_i, :) = lStim.sPars.int;
        wDat.sPars.freq(rep_i, :) = lStim.sPars.freq;
        wDat.sPars.width(rep_i, :) = lStim.sPars.width;
        wDat.sPars.sr(rep_i, :) = lStim.sPars.sr;
        wDat.sPars.basPre(rep_i, :) = lStim.sPars.basPre/unit_cor_factor;
        wDat.sPars.basPost(rep_i, :) = lStim.sPars.basPost/unit_cor_factor;
    else
        fprintf('Error lStim variable is missing fields\n')
    end
    
    wDat.sPars.led_mini(rep_i, 1) = cDat.minInit;
    wDat.sPars.led_mine(rep_i, 1) = cDat.minEnd;
    wDat.sPars.led_delta(rep_i, 1) = eval(['median(cDat.', cDat.CorType, ')']);
    
end

% Collect motion related
if ~isempty(mcDat)
    
    wDat.MotCor.sYshift(:, rep_i) = mcDat.rigid(1, :)';
    wDat.MotCor.sXshift(:, rep_i) = mcDat.rigid(2, :)';
    
    if size(mcDat.rigid, 1) == 3
        wDat.MotCor.Zshift(:, rep_i) = mcDat.rigid(3, :)';
    else
        wDat.MotCor.Zshift(:, rep_i) = zeros(size(wDat.MotCor.sYshift, 1), 1);
    end
    
end

% Stitch related
wDat.Zstitch.Zidx = cat(1, wDat.Zstitch.Zidx, ones(zlength, 1)*reps(rep_i));
if rep_i > 1
    wDat.Zstitch.Yshift(rep_i-1, 1) = shifts_align.shifts(1, 1, 1, 1);
    wDat.Zstitch.Xshift(rep_i-1, 1) = shifts_align.shifts(1, 1, 1, 2);
end

% save substracted background
if isfield(iDat, 'bs'); wDat.bF(rep_i, :) = iDat.bs(:); end 

% collect PMT_fscore (RedChannel)
if isfield(iDat, 'PMT_fscore')
    wDat.PMT_fscore(rep_i, :) = iDat.PMT_fscore;
end

wDat.power(rep_i, 1) = iDat.Power;

end