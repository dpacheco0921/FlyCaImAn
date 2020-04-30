function [Y_, idx] = build_matrix_with_lags(...
    Y, filter_onset_offset, pad_val)
% build_matrix_with_lags: using filter onset and offset timepoints
%   it builds the stimulus matrix to be used for model fitting
%
% Usage:
%   build_matrix_with_lags(Y, filter_onset_offset, pad_val)
%
% Args:
%   Y: stimuli matrix [stimuli_n x time]
%   filter_onset_offset: timepoints before and after stimuli to use for
%       matrix building
%   pad_val: value to pad stimulus ends with before circshit

Y_ = cell(size(Y, 1), 1);

for i = 1:size(Y, 1)
    
    padsize = sum(abs(filter_onset_offset(i, :)));
    Y_in = Y(i, :);
    Y_in = padarray(Y_in, [0 padsize], nan, 'both');
    
    % do pre-stim
    k = 1;
    if filter_onset_offset(i, 1) ~= 0
        for j = filter_onset_offset(i, 1):1:-1
            Y_{i}(k, :) = circshift(Y_in, j);
            k = k + 1;
        end
    end

    % do lag 0
    Y_{i}(k, :) = Y_in;
    k = k + 1;
    
    % do post stim
    if filter_onset_offset(i, 2) ~= 0
        for j = 1:1:(filter_onset_offset(i, 2)-1)
            Y_{i}(k, :) = circshift(Y_in, j);
            k = k + 1;
        end
    end
    
    Y_{i, 1} = flip(Y_{i}(:, padsize + 1:end - padsize), 1);
    
    idx{i, 1} = ones(k - 1, 1)*i;
    
end

Y_ = cell2mat(Y_);
idx = cell2mat(idx);

Y_(isnan(Y_)) = pad_val;
Y_ = Y_';

end
