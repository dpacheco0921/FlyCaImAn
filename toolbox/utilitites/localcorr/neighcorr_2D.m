function rho = neighcorr_2D(Y, iparams)
% neighcorr_2D: Function that uses lc2Dgen to generate a correlation image
%   based on neighboring pixels, but in this case it deals better with big
%   matrices (memmap).
%
% Usage:
%   rho = neighcorr_2D(Y, iparams)
%
% Args:
%   Y: M x N x T movie
%   iparams: parameter variable
%
% Returns:
%   rho: M x N matrix, cross-correlation with adjacent pixel
%
% Notes: 
% splits volume in n parts and runs lc2Dgen on each

neigh_2d.overlap = [2 2];
neigh_2d.chunksize = [50 50];
neigh_2d.timestamps = [];

neigh_2d = loparam_updater(neigh_2d, iparams);

memmaped = isobject(Y);

if memmaped
    sizY = Y.sizY;
    patches = construct_patches(sizY(1:end-1), ...
        neigh_2d.chunksize, neigh_2d.overlap);
else
    sizY = size(Y);
end

if memmaped
    
    fprintf('Running neighcorr_2D for big matrices ')
    
    for p_i = 1:length(patches)
        xi = patches{p_i}(1):patches{p_i}(2);
        yi = patches{p_i}(3):patches{p_i}(4);
        temp_rho{p_i} = lc2Dgen(Y.Y(xi, yi, :), neigh_2d.timestamps);
        fprintf('*')
    end
    
    fprintf(' ')
    % Correct interseccions
    rho = seampatch(patches, temp_rho, sizY);
    fprintf('\n')
    
else
    
    rho = lc2Dgen(Y);
    
end

end

function rho = lc2Dgen(data, timestamps)
% lc2Dgen: get local correlation of the calcium imaging data (8 neighboring pixels)
%
%
% Usage:
%   rho = lc2Dgen(data, timestamps)
%
% Args:
%   data: M x N x T movie
%   timestamps: parameter variable
%
% Returns:
%   rho: M x N matrix, cross-correlation with adjacent pixel
%
% Author: Yuanjun Gao

if exist('timestamps', 'var') && ~isempty(timestamps)
    data = data(:, :, timestamps);
end

[M, N, ~] = size(data);
rho = zeros(M, N);

%normalizing data
data = bsxfun(@minus, data, mean(data, 3));
dataStd = sqrt(mean(data.^2, 3));
data = bsxfun(@rdivide, data, dataStd);

%calculate cross-correlation
rho1 = mean(data(1:(end - 1), :, :) .* data(2:end, :, :), 3);
rho2 = mean(data(:, 1:(end - 1), :) .* data(:, 2:end, :), 3);
rho3 = mean(data(2:end, 2:end, :) .* data(1:(end - 1), 1:(end - 1), :), 3);
rho4 = mean(data(1:(end - 1), 2:end, :) .* data(2:end, 1:(end - 1), :), 3);

%add all cross-correlation
rho(1:(end - 1), :) = rho(1:(end - 1), :) + rho1;
rho(2:end, :) = rho(2:end, :) + rho1;
rho(:, 1:(end - 1)) = rho(:, 1:(end - 1)) + rho2;
rho(:, 2:end) = rho(:, 2:end) + rho2;
rho(1:(end - 1), 1:(end - 1)) = rho(1:(end - 1), 1:(end - 1)) + rho3;
rho(2:end, 2:end) = rho(2:end, 2:end) + rho3;
rho(1:(end - 1), 2:end) = rho(1:(end - 1), 2:end) + rho4;
rho(2:end, 1:(end - 1)) = rho(2:end, 1:(end - 1)) + rho4;

%normalize by number of adjacent pixels
nNgbr = 8 + zeros(M, N);
nNgbr([1,end], :) = nNgbr([1,end], :) - 3;
nNgbr(:, [1,end]) = nNgbr(:, [1,end]) - 3;
nNgbr(1, 1) = nNgbr(1, 1) + 1;
nNgbr(1, end) = nNgbr(1, end) + 1;
nNgbr(end, 1) = nNgbr(end, 1) + 1;
nNgbr(end, end) = nNgbr(end, end) + 1;

rho = rho ./ nNgbr;

rho(isnan(rho)) = 0;

end

