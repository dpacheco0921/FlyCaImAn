function [P,Y] = preprocess_data_int(Y,p,options)

% data pre-processing for:
% (i)   identifying and interpolating missing entries (assumed to have the
%       value NaN). Interpolated entries are passed back to Y.
% (ii)  identifying saturated pixels
% (iii) estimating noise level for every pixel
% (iv)  estimating global discrete time constants (if needed)
% This function replaces arpfit, present in the previous versions of the code.

% Author: Eftychios A. Pnevmatikakis
%           Simons Foundation, 2015

defoptions = CNMFSetParms;

if nargin < 3 || isempty(options); options = defoptions; end
if nargin < 2 || isempty(p); p = 2; end
P.p = p;

if ~isfield(options,'noise_range'); options.noise_range = defoptions.noise_range; end
if ~isfield(options,'noise_method'); options.noise_method = defoptions.noise_method; end
if ~isfield(options,'block_size'); options.block_size = defoptions.block_size; end
if ~isfield(options,'flag_g'); options.flag_g = defoptions.flag_g; end
if ~isfield(options,'lags'); options.lags = defoptions.lags; end
if ~isfield(options,'include_noise'); options.include_noise = defoptions.include_noise; end; include_noise = options.include_noise;
if ~isfield(options,'split_data'); split_data = defoptions.split_data; else split_data = options.split_data; end
if ~isfield(options,'cluster_pixels'); cluster_pixels = defoptions.cluster_pixels; else cluster_pixels = options.cluster_pixels; end
if ~isfield(options,'extract_max'); extract_max = defoptions.extract_max; else extract_max = options.extract_max; end
if ~isfield(options,'max_nlocs'); options.max_nlocs = defoptions.max_nlocs; end
if ~isfield(options,'max_width'); options.max_width = defoptions.max_width; end
if ~isfield(options,'sr'); options.sr = 4; end % sampling rate of data Hz
if ~isfield(options,'freq_th'); options.freq_th = 0.5; end % Hz over which to log-avegare

%% interpolate missing data

if any(isnan(Y(:)))
    Y_interp = interp_missing_data(Y);      % interpolate missing data
    mis_data = find(Y_interp);
    Y(mis_data) = Y_interp(mis_data);       % introduce interpolated values for initialization
else
    Y_interp = sparse(size(Y));
    mis_data = [];
end
P.mis_values = full(Y_interp(mis_data));
P.mis_entries = mis_data;

%% indentify saturated pixels

P.pixels = find_unsaturatedPixels(Y);                % pixels that do not exhibit saturation

%% estimate noise levels

fprintf('Estimating the noise power for each pixel from a simple PSD estimate...');
[P.sn, P.psdx, ff] = get_noise_fft(Y,options);
P.sn = P.sn(:);
fprintf('  done \n');

%% Reducing psd to relevant frequencies: average all F above options.freq_th, 0.5Hz 
ff_l = ff*options.sr;
ff_idx = find(ff_l > options.freq_th, 1, 'first') + 1;
lofreq = sqrt(exp(mean(log(P.psdx(:, ff_idx:end)/2), 2)));
P.psdx(:, (ff_idx+1):end) = [];
P.psdx(:, ff_idx) = lofreq;
P.psdx = single(P.psdx);
% For future use of psdx consider doing the following:
%P.psdx = sqrtbigmem(P.psdx(:,3:end-1));
%P.psdx = centerbigmem(P.psdx);
%P.psdx = sdnormbigmem(P.psdx);

%% extract maximum activity for each pixel
if extract_max
    [LOCS,Ym] = extract_max_activity(Y,options.max_nlocs,options.max_width);
    P.max_locs = LOCS;
    P.max_data = Ym;
end

fprintf('  done \n');
end