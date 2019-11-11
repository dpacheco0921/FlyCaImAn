function [A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches_int_compile(data,K,patches,tau,p,options,filename)
% For detail description look at original code run_CNMF_patches
%% Loading data and setting some parameters
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
        data = matfile('data_file.mat','Writable',true);
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

Yc = cell(length(patches),1);
if ~memmaped
    for i = 1:length(patches)
        if length(sizY) == 3
            Yc{i} = Y(patches{i}(1):patches{i}(2), ...
                patches{i}(3):patches{i}(4),:);
        else
            Yc{i} = Y(patches{i}(1):patches{i}(2), ...
                patches{i}(3):patches{i}(4), ...
                patches{i}(5):patches{i}(6),:);
        end
    end
end        

if nargin < 2 || isempty(K)
    K = 10;
end

RESULTS(length(patches)) = struct('mC', [], ...
    'mCi', [], 'A', [], 'C', [], 'b', [], ...
    'f', [], 'S', [], 'P', []);
cDir = pwd;

%% Load all RESULTS and compile them
fprintf('Loading generated data...');
f2load = rdir(['.', filesep, filename, '_t_*.mat']);
f2load = {f2load.name};
f2load = strrep(f2load, ['.', filesep], '');
f2load = strrep(f2load, '.mat', '');
repnum = zeros(1, numel(f2load));
for ii = 1:numel(f2load)
    TempS = strsplit2(f2load{ii}, '_t_');
    repnum(1, ii) = str2double(TempS{2});
end
repnum = sort(repnum);

for ij = repnum
    fprintf('*')
    lRESULTS = load([cDir, filesep, filename, '_t_', num2str(ij), '.mat'], 'RESULTS');
    %% dF thresholding
    idx2keep = threshold_C(lRESULTS(1).RESULTS.C, lRESULTS(1).RESULTS.A, lRESULTS(1).RESULTS.YrA, options);
    [~, i2del] = setdiff(lRESULTS(1).RESULTS.mCi, idx2keep);
    lRESULTS(1).RESULTS.mC(i2del) = [];
    lRESULTS(1).RESULTS.mCi(i2del) = [];
    lRESULTS(1).RESULTS.A = lRESULTS(1).RESULTS.A(:, idx2keep);
    lRESULTS(1).RESULTS.C = lRESULTS(1).RESULTS.C(idx2keep, :);
    lRESULTS(1).RESULTS.S = lRESULTS(1).RESULTS.S(idx2keep, :);
    if isfield(lRESULTS(1).RESULTS, 'YrA')
        lRESULTS(1).RESULTS = rmfield(lRESULTS(1).RESULTS, 'YrA');
    end
    RESULTS(ij) = lRESULTS(1).RESULTS;
    clear lRESULTS
end
fprintf(' done. \n');

%% Combine results into one structure
fprintf('Combining results from different patches...');
tic
d = prod(sizY(1:end-1));
A = sparse(d,length(patches)*K);
P.sn = zeros(sizY(1:end-1));
if isfield(RESULTS(1).P,'sn_ds'); P.sn_ds = zeros(sizY(1:end-1)); end
P.b = {};
P.c1 = {};
P.gn = {};
P.neuron_sn = {};

if length(sizY) == 3
    P.psdx = zeros(patches{end}(2),patches{end}(4),size(RESULTS(1).P.psdx,2));
else
    P.psdx = []; % do not collect anything if 3DxT
    % P.psdx = zeros(patches{end}(2),patches{end}(4),patches{end}(6),size(RESULTS(1).P.psdx,2));
end

cnt = 0;
B = sparse(prod(sizY(1:end-1)),length(patches));
MASK = zeros(sizY(1:end-1));
F = zeros(length(patches),sizY(end));
% get number of components per RESULTS
mC = []; mCi = [];
numC = cellfun(@(x) size(x,1), {RESULTS(:).C}); numC = [0 cumsum(numC)];

for i = 1:length(patches)
    for k = 1:K
        if k <= size(RESULTS(i).A,2)
            cnt = cnt + 1;
            Atemp = zeros(sizY(1:end-1));
            if length(sizY) == 3
                Atemp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = ...
                    reshape(RESULTS(i).A(:,k),patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);            
            else
                Atemp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) = ...
                    reshape(full(RESULTS(i).A(:,k)),patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,patches{i}(6)-patches{i}(5)+1);
            end
            A(:, cnt) = sparse(Atemp(:));
        end
    end
    if length(sizY) == 3
        b_temp = sparse(sizY(1),sizY(2));
        b_temp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = ...
            reshape(RESULTS(i).b,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);
        MASK(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = ...
            MASK(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) + 1;
        P.sn(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = ...
            reshape(RESULTS(i).P.sn,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);
        if isfield(RESULTS(i).P,'sn_ds')
            P.sn_ds(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = ...
                reshape(RESULTS(i).P.sn_ds,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);
        end
        P.psdx(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),:) = ...
            reshape(RESULTS(i).P.psdx,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,[]);
    else
        b_temp = zeros(sizY(1),sizY(2),sizY(3));
        b_temp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) = ...
            reshape(full(RESULTS(i).b),patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,patches{i}(6)-patches{i}(5)+1);
        MASK(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) = ...
            MASK(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) + 1;
        P.sn(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) = ...
            reshape(RESULTS(i).P.sn,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,patches{i}(6)-patches{i}(5)+1);
        if isfield(RESULTS(i).P,'sn_ds')
            P.sn_ds(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) = ...
                reshape(RESULTS(i).P.sn_ds,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,patches{i}(6)-patches{i}(5)+1);
        end      
        %P.psdx(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6), :) = ...
        %    reshape(RESULTS(i).P.psdx,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1, patches{i}(6)-patches{i}(5)+1,[]);
    end
    P.b = [P.b;RESULTS(i).P.b];
    P.c1 = [P.c1;RESULTS(i).P.c1];
    P.gn = [P.gn;RESULTS(i).P.gn];
    P.neuron_sn = [P.neuron_sn;RESULTS(i).P.neuron_sn];
    B(:,i) = sparse(b_temp(:));
    F(i,:) = RESULTS(i).f;
    % dealing with traces of merged components
    % i)collect all the cells and ii) get the global index of the merged component
    mC = [mC; RESULTS(i).mC];
    mCi = [mCi, RESULTS(i).mCi + numC(i)];
end
A(:,cnt+1:end) = [];
A = spdiags(1./MASK(:),0,prod(sizY(1:end-1)),prod(sizY(1:end-1)))*A;
B = spdiags(1./MASK(:),0,prod(sizY(1:end-1)),prod(sizY(1:end-1)))*B;
C = cell2mat({RESULTS(:).C}');
S = cell2mat({RESULTS(:).S}');
ff = find(sum(A,1)==0);
A(:,ff) = [];
C(ff,:) = [];
S(ff,:) = [];
[mCi, mCi2del] = update_mCi(mCi, ff);
mC(mCi2del) = [];
mCi(mCi2del) = [];
fprintf(' done. \n');
toc

clear MASK b_temp Atemp ff defoptions k K

RESULTS = [];
%% Merge results
fprintf('Merging overlaping components...\n')
tic
Km = 0;
Kn = size(A,2);
while Km < Kn
    Kn = size(A,2);
    fprintf(['Initial number of components : ', num2str(Kn), '\n'])
    Cinit = C;
    [A,C,~,imCi,P,S] = merge_components([],A,[],C,[],P,S,options);
    % update the cells that stores merged C, and their indexes to the C they generate
    [mC, mCi] = update_mC(mC, mCi, Cinit, C, imCi);
    Km = size(A,2);
    fprintf(['New number of components : ', num2str(Km), '\n'])
end
clear Km Kn S
fprintf(' done. \n');
toc

%% compute spatial and temporal background using a rank-1 fit
fprintf('Computing background components...')
fin = [mean(F, 1);rand(options.gnb-1,length(F))];
for iter = 1:150
    fin = diag(sqrt(sum(fin.^2,2)))\fin;
    bin = max(B*(F*fin')/(fin*fin'),0);
    fin = max((bin'*bin)\(bin'*B)*F,0);
end
clear B
fprintf(' done. \n');

%% update spatial components
fprintf('Updating spatial components...');
tic
options.nb = options.gnb;
options.d1 = sizY(1);
options.d2 = sizY(2);
if length(sizY) == 4; options.d3 = sizY(3); end
if ~isfield(P, 'mis_values'); P.mis_values = []; end
if ~isfield(P, 'mis_entries'); P.mis_entries = []; end
options.spatial_method = 'constrained';
options.spatial_parallel = 1;
[A,b,C,P] = update_spatial_components(data,C,fin,[A,bin],P,options);
fprintf(' done. \n');
toc
%% update temporal components
fprintf('Updating temporal components... ')
tic
P.p = p;
options.temporal_iter = 2;
options.temporal_parallel = 1;
[C,f,P,S,YrA] = update_temporal_components(data,A,b,C,fin,P,options);
fprintf(' done. \n');
toc

%% remove small components
% 3by3 xy matrix, it could even be harder (18 or 27, given that rois usually spans many planes
ff = find(sum(A~=0) < options.Asize_ths);
fprintf(['Deleting small ROIs: ', num2str(numel(ff)), '\n'])
A(:,ff) = [];
C(ff,:) = [];
YrA(ff, :) = [];
S(ff, :) = [];
[mCi, mCi2del] = update_mCi(mCi, ff);
mC(mCi2del) = [];
mCi(mCi2del) = [];

%% Save additional info
P.f_perseg = F; % the temporal f
P.mC = mC;
P.mCi = mCi; 
% clear variable to reduce memory requirements, consider not collecting it to begin with
P.psdx = [];

fprintf(' done. \n');
end