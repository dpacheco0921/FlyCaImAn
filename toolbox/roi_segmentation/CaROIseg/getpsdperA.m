function mean_psdx = getpsdperA(A, P)
% get mean PSD per spatial component
% flatten psd
ndim = size(P.psdx);
P.psdx = reshape(P.psdx, prod(ndim(1:end-1)), ndim(end));
for i = 1:size(A,2)
    pixIdx = find(A(:, i)~= 0);
    mean_psdx(i, :) = mean(P.psdx(pixIdx, :), 1);
    clear pixIdx
end