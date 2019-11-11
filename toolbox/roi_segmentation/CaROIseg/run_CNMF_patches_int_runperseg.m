function run_CNMF_patches_int_runperseg(data, K, patches, tau, p, options, pacthIdx, filename)
% run each segment and save RESULTS
defoptions = CNMFSetParms;

if nargin < 6 || isempty(options)
    options = defoptions;
end

if ~isfield(options, 'cluster_pixels') || isempty(options.cluster_pixels)
    options.cluster_pixels = defoptions.cluster_pixels;
end
if ~isfield(options, 'create_memmap') || isempty(options.create_memmap)
    options.create_memmap = defoptions.create_memmap;
end
if ~isfield(options, 'gnb') || isempty(options.gnb)
    options.gnb = defoptions.gnb;
end
if ~isfield(options, 'classify_comp') || isempty(options.classify_comp)
    options.classify_comp = defoptions.classify_comp;
end

memmaped = isobject(data);
if memmaped
    sizY = data.sizY;
else    % create a memory mapped object named data_file.mat    
    Y = data;
    clear data;
    sizY = size(Y);
    Yr = reshape(Y, prod(sizY(1:end-1)), []);
    nY = min(Yr(:));
    %Yr = Yr - nY;
    if options.create_memmap
        save('data_file.mat', 'Yr', 'Y', 'nY', 'sizY', '-v7.3');
        data = matfile('data_file.mat', 'Writable',true);
        memmaped = true;
    else
        data = Yr;
    end
end

if ~isfield(options, 'd1') || isempty(options.d1); options.d1 = sizY(1); end
if ~isfield(options, 'd2') || isempty(options.d2); options.d2 = sizY(2); end
if ~isfield(options, 'd3') || isempty(options.d3)
    if length(sizY) == 3
        options.d3 = 1;
    else
        options.d3 = sizY(3);
    end
end

if nargin < 5 || isempty(p)
    p = 0;
end

if nargin < 4 || isempty(tau)
    tau = 5;
end

if nargin < 3 || isempty(patches)
    patches = construct_patches(sizY(1:end-1), [50, 50]);
end

Yc = cell(length(patches), 1);
if ~memmaped
    for j = 1:length(patches)
        if length(sizY) == 3
            Yc{j} = Y(patches{j}(1):patches{j}(2), ...
                patches{j}(3):patches{j}(4), :);
        else
            Yc{j} = Y(patches{j}(1):patches{j}(2), ...
                patches{j}(3):patches{j}(4), ...
                patches{j}(5):patches{j}(6), :);
        end
    end
end        

if nargin < 2 || isempty(K)
    K = 10;
end

RESULTS = struct();
cDir = pwd;

% Process each patch
fprintf('Running all patches individually\n')

if length(sizY) == 3
    if memmaped
        Y = data.Y(patches{pacthIdx}(1):patches{pacthIdx}(2), ...
            patches{pacthIdx}(3):patches{pacthIdx}(4), :);
    else
        Y = Yc{pacthIdx};
    end
    [d1, d2, T] = size(Y);
    d3 = 1;
else
    if memmaped
        Y = data.Y(patches{pacthIdx}(1):patches{pacthIdx}(2), ...
            patches{pacthIdx}(3):patches{pacthIdx}(4), ...
            patches{pacthIdx}(5):patches{pacthIdx}(6), :);
    else
        Y = Yc{pacthIdx};
    end
    [d1, d2, d3, T] = size(Y);
end

Y = double(Y);
d = d1*d2*d3;
options_temp = options;

% chop the brain mask if it exist
if isfield(options, 'brainmask') && ~isempty(options.brainmask)
    if length(sizY) == 3
        options_temp.brainmask = ...
            options.brainmask(patches{pacthIdx}(1):patches{pacthIdx}(2), ...
            patches{pacthIdx}(3):patches{pacthIdx}(4));
    else
        options_temp.brainmask = ...
            options.brainmask(patches{pacthIdx}(1):patches{pacthIdx}(2), ...
            patches{pacthIdx}(3):patches{pacthIdx}(4), ...
            patches{pacthIdx}(5):patches{pacthIdx}(6));
    end
    options_temp.brainmask = options_temp.brainmask(:); 
end

options_temp.d1 = d1;
options_temp.d2 = d2;
options_temp.d3 = d3;
options_temp.nb = 1;

[P, Y] = preprocess_data_int(Y,p,options_temp);
[Ain, Cin, bin, fin] = initialize_components_int(Y,K,tau,options_temp,P);  % initialize
Yr = reshape(Y, d, T);

%clear Y;
options_temp.spatial_parallel = 0; % turn off parallel updating for spatial components
[A,b,Cin,P] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options_temp);

P.p = 0;
options_temp.temporal_parallel = 0; % turn off parallel updating for temporal components
[C,f,P,S] = update_temporal_components(Yr,A,b,Cin,fin,P,options_temp); 

fprintf(['Initial number of components : ', num2str(size(C, 1)), '\n'])
[Am,Cm,~,mC,P] = merge_components(Yr,A,b,C,f,P,S,options_temp);

% Collect all merged components
[RESULTS.mC, RESULTS.mCi] = getCfrommC(mC, C, Cm);
fprintf(['New number of components : ', num2str(size(Cm, 1)), '\n'])

% threshold components
fprintf(['Deleting components with a delta C lower than ', num2str(options_temp.l_df), '\n'])
idx2keep = threshold_C(Cm, Am, [], options_temp);
Am = Am(:, idx2keep);
Cm = Cm(idx2keep, :);

% update mC and mCi
[~,i2del] = setdiff(RESULTS.mCi, idx2keep);
RESULTS.mC(i2del) = []; RESULTS.mCi(i2del) = [];
fprintf(['New number of components : ', num2str(size(Cm, 1)), '\n'])

if ~isempty(Am)
    
    % final spatial and temporal update
    [A2,b2,Cm,P] = update_spatial_components(Yr,Cm,f,[Am,b],P,options_temp);
    P.p = p;
    [C2,f2,P2,S2,YrA] = update_temporal_components(Yr,A2,b2,Cm,f,P,options_temp);
    
else
    
    % reruning initialization so we estimate b and f accurately
    K = 1;
    [P,Y] = preprocess_data_int(Y,p,options_temp);
    [Am,Cm,b,f] = initialize_components_int(Y,K,tau,options_temp,P);  % initialize
    Yr = reshape(Y,d,T);
    options_temp.spatial_parallel = 0; % turn off parallel updating for spatial components
    [A2,b2,~,P] = update_spatial_components(Yr,Cm,f,[Am,b],P,options_temp);
    P.p = p;
    options_temp.temporal_parallel = 0; % turn off parallel updating for temporal components
    [C2,f2,P2,S2,YrA] = update_temporal_components(Yr,A2,b2,Cm,f,P,options_temp);
    
end

RESULTS.A = A2;
RESULTS.C = C2;
RESULTS.b = b2;
RESULTS.f = f2;
RESULTS.S = S2;
RESULTS.P = P2;
RESULTS.YrA = YrA;

fprintf(['Finished processing patch # ', num2str(pacthIdx), ...
    ' out of ', num2str(length(patches)), '.\n']);

% save tempfile
save([cDir, filesep, filename, '_t_', num2str(pacthIdx), '.mat'], 'RESULTS')

end