function Data2DxT = interp2DxT(Data2DxT, ...
    XYinit, XYend, intmethod)
% interp2DxT: resampling 2D x T volumes or 3D
%
% Usage:
%   Data2DxT = interp2DxT(Data2DxT, ...
%       XYinit, XYend, method)
%
% Args:
%   Data2DxT: 2DxT matrix
%   XYinit: resolution of input data
%   XYend: target resolution of data
%   intmethod: interpolation method
%       ('linear', 'nearest', 'cubic', 'spline')
%       (default, 'linear')
%
% Notes:
% units in microns order: [X Y T] which is:
%   (height, width, time), XYT and XYT are (1, 3) vectors

XYinit = XYinit(1:2); XYend = XYend(1:2);
fprintf('Resampling data ... ')
Data2DxT = single(Data2DxT);
Dims = double(size(Data2DxT));
nDim = length(Dims);
DimEnd = (Dims(1:2)-1).*XYinit;
% last dimention to interpolate
if ~exist('intmethod', 'var'); intmethod = 'linear'; end

% Order of dimentions to resample
AxesOrder = [1 2 3; 2 1 3; 3 1 2];  %(height, width)
AxesOrder2 = [1 2 3; 2 1 3; 2 3 1];  %(height, width)
AxesOrder = AxesOrder(:, 1:nDim);
AxesOrder2 = AxesOrder2(:, 1:nDim);

% resample
for Dim_i = 1:2
    
    % do linear interpolation of each dimention separately from 1 to 3 (3rd is time)
    fprintf('*')
    Axes_i = 0:XYinit(Dim_i):DimEnd(Dim_i); 
    Axes_e = 0:XYend(Dim_i):DimEnd(Dim_i);
    
    if XYinit(Dim_i) ~= XYend(Dim_i)
        fprintf('*')
        Data2DxT = permute(interp1(Axes_i, permute(Data2DxT, AxesOrder(Dim_i, :)), ...
            Axes_e, intmethod), AxesOrder2(Dim_i,:));
    end
    
    clear Init_Axes New_Axes
    
end

fprintf(' ')

end


