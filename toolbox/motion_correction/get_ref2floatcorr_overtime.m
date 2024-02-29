function CM = get_ref2floatcorr_overtime(templateIm, floatIm)
% get_ref2floatcorr_overtime: calculate correlation of template to each frame of floating image
%
% Usage:
%   CM = get_ref2floatcorr_overtime(templateIm, floatIm)
%
% Args:
%   templateIm: reference image to use
%   floatIm: floating image, its last dimension is time.
%       The size of the matrix is [size(templateIm), time].
%
% Returns:
%   CM: pearson correlation overtime

siz_ = size(floatIm);
templateIm = reshape(templateIm, [prod(siz_(1:end-1)) 1]);
floatIm = reshape(floatIm, [prod(siz_(1:end-1)) siz_(end)]);

CM = corr(templateIm, floatIm)';

end