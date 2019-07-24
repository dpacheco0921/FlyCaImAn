function Data1DxT = interp1DxTixTj(Data1D, iniT, endT, intmethod)
% interp1DxTixTj: This function interpolates the calcium matrix M (n, t) 
%   usint Time matrix iniT
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
% Notes

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
