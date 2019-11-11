function [inferred, filtered, raw, filtered_, raw_] = ...
    signalExtraction_int(Y, Y_, A, C, b, f, d1, d2, d3, iparams)
% [inferred, filtered, raw, filtered_, raw_] = ...
%    signalExtraction_int(Y, Y_, A, C, b, f, d1, d2, d3, iparams)
%
% inputs: Y raw data (d X T matrix, d # number of pixels, T # of timesteps)
%         Y_ raw data from extra channel (d X T matrix, d # number of pixels, T # of timesteps)
%         A matrix of spatial components (d x K matrix, K # of components)
%         C matrix of temporal components (K x T matrix)
%         S matrix of deconvolved activity ((K-1) x T matrix) (optional)
%         b matrix of background spatial components (d x options.nb, options.nb # of components)
%         f matrix of background temporal components (options.nb x T matrix, options.nb # of components) 
%         iparams.b_prc: percentile of baseline of signals from CNMF to be substracted(default:0.25);
%           set to 0 to ignore this residue baseline 
%         iparams.thr: energy fraction included in the raw signal contour (default:0.99)
%
% outputs: inferred: spatial weighting on the ROI pixels, unmixing, background substraction, and denoising
%          filtered: spatial weighting on the ROI pixels, unmixing, background substraction
%          raw: uniform spatial weighting on the ROI pixels (with threshold to remove the very low energy pixels),
%           background substraction
% 
% Notes:
% df: df signal
% F: baseline
% dfof: dfof
% Ftrend: trend dusbtracted from C and added to F
% CR: information of the cell contour used in raw signal calculation
% Author: Eftychios A. Pnevmatikakis, and Weijian Yang
% edited from signalExtraction_int
% Notes:

% Default params
% detrend signal
p.dtype = 1; % 0 = lowpass filter, 1 = percentile filter
p.sr = 2;
p.lfreq = 0.002;
p.d_prct = 20;
p.d_window = 100;
p.d_shift = 10;
% A thresholding
p.medw = [3 3];
p.nrgthr = 0.99;
p.clos_op = strel('square', 3);
% background estimation
p.df_window = [];
p.df_prctile = [];
% gate to run on parallel
p.se_parallel = true; 

if ~exist('iparams', 'var'); iparams = []; end
p = loparam_updater(p, iparams);

% normalize spatial components to unit energy
fprintf('1) Normalize A and C\n')
A2 = [A, b];
C2 = [C; f];
K = size(C, 1);
K2 = size(C2, 1);
nA = sqrt(sum(A2.^2))';
% normalize spatial components to unit energy
A2 = A2/spdiags(nA, 0, K2, K2);
C2 = spdiags(nA, 0, K2, K2)*C2;
A = A2(:, 1:K);
C = C2(1:K, :);
b = A2(:, K+1:end);
f = C2(K+1:end, :);

% get extra values
fprintf('2) get AY\n')
tic; AY = mm_fun(A, Y, [], p.se_parallel, []); toc

% background + residue
Yf = AY - (A'*A)*C;

% residue + C
Y_r = Yf - (A'*full(b))*f + C;

% detrend C or C + residue: low pass, and detrend or use runing percentile
fprintf('3) detrend Y_r (C + Yr) and C\n')
if p.dtype == 0
    
    [filt_b, filt_a] = butter(2, p.lfreq/(p.sr/2), 'low');
    C_Yr_detred = filtfilt(filt_b, filt_a, Y_r')';
    C_detrend = filtfilt(filt_b, filt_a, C')';
    
else
    
    C_Yr_detred = prctfilt(Y_r, p.d_prct, p.d_window, p.d_shift, 0);
    C_detrend = prctfilt(C, p.d_prct, p.d_window, p.d_shift, 0);
    
end

fprintf('4) generate inferred and filtered signals\n')

% inferred signal
inferred.F = getFo(Yf, p.df_prctile, p.df_window) + C_detrend;
inferred.df = C - C_detrend;
inferred.dfof = spdiags(inferred.F, 0, K, K)\inferred.df;
inferred.Ftrend = C_detrend;

% filtered signal
filtered.F = getFo(Yf - Y_r + C, p.df_prctile, p.df_window) + C_Yr_detred;
filtered.df = Y_r - C_Yr_detred;
filtered.dfof = spdiags(filtered.F, 0, K, K)\filtered.df;
filtered.Ftrend = C_Yr_detred;
clear C_Yr_detred C_detrend AY

% ********** red signal **********
if ~isempty(Y_) 
    
    % processed extra channel
    % Notes: for the red signal, the F is basically any trend find in that
    % signal, and the signal or df is the detrended signal.
    
    fprintf('5) do steps 2 and 3 for extra channel and calculate filtered signal\n')
    tic; AY_ = mm_fun(A, Y_, [], p.se_parallel, []); toc
    
    % detrend C or C + residue: low pass, and detrend or use runing percentile
    if p.dtype == 0
        [filt_b, filt_a] = butter(2, p.lfreq/(p.sr/2), 'low');
        Yr_detred_ = filtfilt(filt_b, filt_a, AY_')';
    else
        Yr_detred_ = prctfilt(AY_, p.d_prct, p.d_window, p.d_shift, 0);
    end
    
    filtered_.F = Yr_detred_;
    filtered_.df = AY_ - Yr_detred_;
    filtered_.dfof = spdiags(filtered_.F, 0, K, K)\filtered_.df; 
    clear AY_ Yr_detred_
    
else
    
    filtered_ = [];
    
end

% ********** raw signal **********
fprintf('6) threshold A and make it uniform (uA)\n')

% spatial contour for raw signal
rawA = A; 
tic; [rawA, CR] = thresholdA(rawA, d1, d2, d3, ...
    p.medw, p.nrgthr, p.clos_op, p.se_parallel); toc
nA = sqrt(sum(rawA.^2))';

% normalize spatial components to unit energy
rawA = rawA/spdiags(nA, 0, K, K);
raw_bprc = (rawA'*full(b))*f;

fprintf('7) get uAY\n')
% no weighting on the ROI pixels, background substraction
tic; raw.df = mm_fun(rawA, Y, [], p.se_parallel, []) - raw_bprc; toc

fprintf('8) detrend uAY\n')
% detrend df: low pass, and detrend then threshold
if p.dtype == 0
    raw_detrend = filtfilt(filt_b, filt_a, raw.df')';
else
    raw_detrend = prctfilt(raw.df, p.d_prct, p.d_window, p.d_shift, 0);
end

fprintf('9) generate raw signals\n')
raw.F = getFo(raw_bprc, p.df_prctile, p.df_window) + raw_detrend;
raw.df = raw.df - raw_detrend;
raw.dfof = spdiags(raw.F, 0, K, K)\raw.df;
raw.Ftrend = raw_detrend;
raw.CR = CR;

% ********** red signal **********
if ~isempty(Y_) 
    
    % processed extra channel
    fprintf('10) do steps 7 and 8 for extra channel and calculate raw signal\n')
    tic; raw_.df = mm_fun(rawA, Y_, [], p.se_parallel, []); toc
    
    if p.dtype == 0
        raw_detrend_ = filtfilt(filt_b, filt_a, raw_.df')';
    else
        raw_detrend_ = prctfilt(raw_.df, p.d_prct, p.d_window, p.d_shift, 0);
    end
    
    raw_.F = raw_detrend_;
    raw_.df = raw_.df - raw_detrend_;
    raw_.dfof = spdiags(raw_.F, 0, K, K)\raw_.df;
    
else
    
    raw_ = [];
    
end

end

function [Ath, CR] = thresholdA(A, d1, d2, d3, ...
    medw, nrgthr, clos_op, se_parallel)

[d, nr] = size(A);
Ath = spalloc(d, nr, nnz(A));
indf = cell(nr, 1);
valf = cell(nr, 1);
CR = cell(nr, 1);

if se_parallel
    
    parfor i = 1:nr
        
        A_temp = reshape(full(A(:, i)), d1, d2, d3);
        
        for z = 1:d3
            A_temp(:, :, z) = medfilt2(A_temp(:, :, z), medw);
        end
        
        A_temp = A_temp(:);
        [temp, ind] = sort(A_temp(:).^2, 'ascend'); 
        temp =  cumsum(temp);
        ff = find(temp > (1-nrgthr)*temp(end), 1, 'first');
        BW = zeros(d1, d2, d3);
        BW(ind(ff:d)) = 1;
        
        for z = 1:d3
            BW(:, :, z) = imclose(BW(:, :, z), clos_op);
        end
        
        [L, NUM] = bwlabeln(BW, 8*(d3 == 1) + 6*(d3 ~= 1));
        
        if NUM > 0
            
            nrg = zeros(NUM, 1);
            
            for l = 1:NUM
                ff = (L == l);
                nrg(l) = sum(A_temp(ff).^2);
            end
            
            [~, indm] = max(nrg);
            ff = find(L == indm);
            indf{i} = ff;
            valf{i} = A_temp(ff);
        
            % Collect xyz indeces
            if d3 > 1
                [ii, jj, kk] = ind2sub([d1, d2, d3], ff);
                CR{i, 1} = [ii, jj, kk]';
            else
                [ii, jj] = ind2sub([d1, d2], ff);
                CR{i, 1} = [ii, jj]';
            end
            
        else
            valf{i} = 0;
        end
        
        if mod(i, 100) == 0
            fprintf('%2.1f%% of ROIs completed \n', i*100/nr);
        end
        
    end
    
else
    
    for i = 1:nr
        
        A_temp = reshape(full(A(:, i)), d1, d2, d3);
        
        for z = 1:d3
            A_temp(:, :, z) = medfilt2(A_temp(:, :, z), medw);
        end
        
        A_temp = A_temp(:);
        [temp, ind] = sort(A_temp(:).^2, 'ascend'); 
        temp =  cumsum(temp);
        ff = find(temp > (1-nrgthr)*temp(end), 1, 'first');
        BW = zeros(d1, d2, d3);
        BW(ind(ff:d)) = 1;
        
        for z = 1:d3
            BW(:, :, z) = imclose(BW(:, :, z), clos_op);
        end
        
        [L, NUM] = bwlabeln(BW, 8*(d3 == 1) + 6*(d3 ~= 1));
        
        if NUM > 0
            
            nrg = zeros(NUM, 1);
            
            for l = 1:NUM
                ff = (L == l);
                nrg(l) = sum(A_temp(ff).^2);
            end
            
            [~, indm] = max(nrg);
            ff = find(L == indm);
            indf{i} = ff;
            valf{i} = A_temp(ff);
            
            % Collect xyz indeces
            if d3 > 1
                [ii, jj, kk] = ind2sub([d1, d2, d3], ff);
                CR{i, 1} = [ii, jj, kk]';
            else
                [ii, jj] = ind2sub([d1, d2], ff);
                CR{i, 1} = [ii, jj]';
            end
            
        else
            valf{i} = 0;
        end
        
        if mod(i, 100) == 0
            fprintf('%2.1f%% of ROIs completed \n', i*100/nr);
        end
        
    end
    
end

for i = 1:nr
    
    if ~isempty(indf{i})
        Ath(indf{i}, i) = valf{i};
    end
    
end

end

function Fo = getFo(signal, df_prctile, df_window)
% Basically computes Fo or background signal

% similar to extract_DF_F
if ~isempty(df_prctile)
    
    if isempty(df_window) || (df_window > size(signal,2))
        Fo = prctile(signal, df_prctile, 2);
    else
        
        Fo = zeros(size(signal));
        
        for i = 1:size(Fo, 1)
            df_temp = running_percentile(signal(i, :), df_window, df_prctile);
            Fo(i, :) = df_temp(:)';
        end
        
    end
    
else
    Fo = signal;
end

if size(Fo, 2) == 1
    Fo = repmat(Fo, 1, size(signal, 2));
end

end