function [Y_output] = addNoise2Signal(SNR, Y_input, reps)
% addNoise2Signal: adds random noise to Y_input
%
% Usage:
%   [Y_output] = addNoise2Signal(SNR, Y_input, reps)
%
% Args:
%   SNR: target signal to noise ratio
%   Y_input: input vector
%   reps: number of copies to generate
%
% Returns:
%   Y_output: cell array with one entry per SNR value, each
%       [length(Y_input), reps] (one noisy copy of Y_input per column)

% amplitude of the signal, used to scale the noise to the requested SNR
signalA = std(Y_input);
noise2use = randn([length(Y_input), reps]);
% scale each noise column by its own population std (std(..., 1))
noise2use = bsxfun(@times, noise2use, std(noise2use, 1));

for i = 1:numel(SNR)
    Y_output{i} = repmat(Y_input', [1, reps]) ...
        + noise2use*signalA/SNR(i);
end

end
