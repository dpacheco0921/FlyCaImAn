function Data1DxT = interp1DxTixTj(Data1D, iniT, endT, intmethod)
% interp1DxTixTj: This function interpolates the calcium matrix M (n, t) 
%   using Time matrix iniT
%
% Usage:
%   Data1DxT = interp1DxTixTj(Data1D, iniT, endT, method)
%
% Args:
%   Data1D: traces per ROI (roi_i, t)
%   iniT: input time
%   endT: output time
%   intmethod: interpolation method
%       ('linear', 'nearest', 'cubic', 'spline')
%       (default, 'linear')
%
% Returns:
%   Data1DxT: interpolated traces [roi_i, numel(endT)] (same units as
%       Data1D); points outside the range of iniT are NaN
%
% Notes
% iniT is [1, t] (shared by all ROIs) or [roi_i, t] (one row per ROI).
% The NaN re-masking at the end assumes numel(endT) == size(Data1D, 2).

if ~exist('intmethod', 'var'); intmethod = 'linear'; end

roiNum = size(Data1D, 1);

Data1DxT = zeros(roiNum, length(endT));

for roi_i = 1:roiNum
    if size(iniT, 1) > 1
        Data1DxT(roi_i, :) = interp1(iniT(roi_i, :), ...
            Data1D(roi_i, :), endT, intmethod);
    else
        Data1DxT(roi_i, :) = interp1(iniT, ...
            Data1D(roi_i, :), endT, intmethod);
    end
end

% do not interpolate NaN values
Data1DxT(isnan(Data1D)) = nan;

end
