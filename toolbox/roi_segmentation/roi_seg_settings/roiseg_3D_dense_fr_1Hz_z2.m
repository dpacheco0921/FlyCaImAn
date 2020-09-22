%% Parameters for 3D data at 1Hz
% default setting for
% voxel size: [1.2 1.2 2] um
% relatively short recording 9-15 min
% sampling rate 1Hz

roiparams = [];
roiparams.init_method = 'greedy';               % 'sparse_NMF'; % 'greedy'
roiparams.merge_thr = 0.90;                     % merging threshold
roiparams.p = 0;                                % order of autoregressive system
                                                %(p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
roiparams.K = 600;                              % number of components to be found
roiparams.patch_size = [300, 300, 11];          % size of each patch along each dimension (optional, default: [32,32])
roiparams.overlap = [10, 10, 3];                % amount of overlap in each dimension (optional, default: [4,4])
roiparams.patch_minsiz = [4, 4, 1];             % minimun size of patch
roiparams.cropseg = 0;                          % lines to prune from border
roiparams.ssub = 1;                             % spatial subsampling   
roiparams.tsub = 1;                             % temporal subsampling
roiparams.croptime = [];
roiparams.cluster_pixels = 1;
roiparams.noise_range = [0.0903 0.5];           % noise is anything from 0.368-2Hz
roiparams.se = strel(ones(6, 6, 5));            % or even 1, try
roiparams.nrgthr = 0.99;                        % threshold for A energy
roiparams.thr_method = 'nrg';
roiparams.noise_norm = 1;
roiparams.fudge_factor = 0.98;
roiparams.temporal_iter = 2;
roiparams.dist = 2;
roiparams.search_method = 'dilate';             % default, 'ellipse'
roiparams.tau = [4 4 2];                        % std of gaussian kernel
roiparams.gSiz = [9 9 5];                       % size of gaussian kernel
roiparams.sr = 1;                               % sampling rate of data Hz
roiparams.spatial_method = 'constrained';       % regularized vs constrained

% parameters related to detrending
roiparams.freq_th = 0.25;                       % Hz over which to log-avegare
roiparams.rem_prct = 10;                        % percentile to be removed
roiparams.df_prctile = [];                      % percentile for dfof for background stimation
roiparams.dtype = 1;                            % 0 = lowpass filter, 1 = percentile filter
roiparams.d_prct = 20;                          % running percentile
roiparams.d_window = 200;                       % window
roiparams.d_shift = 10;                         % shift
roiparams.se_parallel = 1;                      % gate for parallel stuff in extractsignal

% size threshold for ROI
roiparams.Asize_ths = 27;                       % Min size of A to consider it a real ROI

%For more parameter to play with check:
%edit CNMFSetParms
