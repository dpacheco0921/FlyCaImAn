function [P] = freqtransform(data, patches, options)
%% Gets the fft profile per pixel
memmaped = isobject(data);
if memmaped
    sizY = data.sizY;
else    % create a memory mapped object named data_file.mat
    Y = data;
    clear data;
    sizY = size(Y);
    nY = min(Y(:));
    Y = Y - nY; % always zero min
    save('data_file.mat','Y','nY','sizY','-v7.3');
    data = matfile('data_file.mat','Writable',true);
end

defoptions = CNMFSetParms;
if nargin < 6 || isempty(options)
    options = defoptions;
end

if nargin < 3 || isempty(patches)
    patches = construct_patches(sizY(1:end-1),[150, 150, 20]);
end

%% Generate a brain mask from pixel-psd profile
% get psx per segment
fprintf(['Getting a brain mask based on pixel-psd profile from ', num2str(numel(patches)), ' patches\n'])
RESULTS(length(patches)) = struct();
tic
for i = 1:length(patches)
    if length(sizY) == 3
        Y = data.Y(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),:);
        [d1,d2,T] = size(Y);
        d3 = 1;
    else
        Y = data.Y(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6),:);
        [d1,d2,d3,T] = size(Y);
    end
    if ~(isa(Y,'single') || isa(Y,'double')); Y = single(Y);  end
    options_temp = options;
    options_temp.d1 = d1; options_temp.d2 = d2; options_temp.d3 = d3;
    RESULTS(i).P = getPSDandSN(Y, options);
    fprintf(['Finished processing patch # ', num2str(i), ' out of ', num2str(length(patches)), '.\n']);
end

fprintf('Compiling outputs')

P.sn = zeros(sizY(1:end-1));
if length(sizY) == 3
    P.psdx = zeros(patches{end}(2),patches{end}(4),size(RESULTS(1).P.psd,2));
else
    P.psdx = zeros(patches{end}(2),patches{end}(4),patches{end}(6),size(RESULTS(1).P.psd,2));
end

% compile P
for i = 1:length(patches)
    if length(sizY) == 3
        P.sn(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),:) = ...
            reshape(RESULTS(i).P.sn,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);
        P.psdx(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),:) = ...
            reshape(RESULTS(i).P.psd,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,[]);
    else
        P.sn(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6), :) = ...
            reshape(RESULTS(i).P.sn,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1, patches{i}(6)-patches{i}(5)+1);
        P.psdx(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6), :) = ...
            reshape(RESULTS(i).P.psd,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1, patches{i}(6)-patches{i}(5)+1,[]);
    end
end
toc
clear RESULTS
%% Downsample
% posibly downsample 2-->1 in xyz
% P.psdx = single(interp3DxT(P.psdx, [1 1 1], [2 2 2], 3));
% update size
% sizY = size(P.psdx);

%% Reduce dimentionality
% linearize P.psdx
% P.psdx = reshape(P.psdx, prod(sizY(1:end-1)), []);

%% do pca
% fprintf('Getting pca\n')
% tic
% if size(P.psdx, 1) < size(P.psdx, 2)
%     [P.PCA.coeff, P.PCA.score, P.PCA.latent, P.PCA.tsquared, P.PCA.explained, P.PCA.mu] = ...
%         pca(P.psdx' , 'NumComponents', 40, 'Algorithm', 'als');
% else
%     [P.PCA.coeff, P.PCA.score, P.PCA.latent, P.PCA.tsquared, P.PCA.explained, P.PCA.mu] = ...
%         pca(P.psdx , 'NumComponents', 40, 'Algorithm', 'als');
% end
% toc
% fprintf('Done\n')

%% do tsne
% fprintf('Getting pca\n')
% tic
% finalDims = 2;
% initialDims = 20;
% perplexity = 10;
% P.tsne = fast_tsne(P.psdx, finalDims, initialDims, perplexity, 0.7);
% toc
% fprintf('Done\n')

%% reshape P
P.psdx = reshape(P.psdx, sizY(1:end-1), []);

fprintf('done\n')
end

function [P] = getPSDandSN(Y, options)
%% Calculate psx and sn
defoptions = CNMFSetParms;
if nargin < 3 || isempty(options); options = defoptions; end
if ~isfield(options,'noise_range'); options.noise_range = defoptions.noise_range; end
if ~isfield(options,'noise_method'); options.noise_method = defoptions.noise_method; end
if ~isfield(options,'block_size'); options.block_size = defoptions.block_size; end
if ~isfield(options,'flag_g'); options.flag_g = defoptions.flag_g; end
if ~isfield(options,'lags'); options.lags = defoptions.lags; end
if ~isfield(options,'sr'); options.sr = 4; end
if ~isfield(options,'freq_th'); options.freq_th = 0.5; end

%% estimate noise levels
fprintf('Estimating the noise power for each pixel from a simple PSD estimate...');
[P.sn, P.psd, ff] = get_noise_fft(Y, options);

%% Reducing psd to relevant frequencies: average all F above options.freq_th, 0.5Hz 
ff_i = ff*options.sr;
idx = find(ff_i > options.freq_th, 1, 'first') + 1;
lofreq = sqrt(exp(mean(log(P.psd(:, idx:end)/2), 2)));
P.psd(:, (idx+1):end) = [];
P.psd(:, idx) = lofreq;

%% preprocess psx
% if ndims(P.psd) > 2
%     sizX = size(P.psd);
%     P.psd = reshape(P.psd, prod(sizX(1:end-1)), sizX(end));
% end
% P.psd = sqrtbigmem(P.psd);
% P.psd = centerbigmem(P.psd);
% P.psd = sdnormbigmem(P.psd);
fprintf('  done \n');
end