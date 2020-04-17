function Data3DxT = interp3DxT(Data3DxT, ...
    XYZinit, XYZend, LastDim, intmethod)
% interp3DxT: resampling 3DxT volumes (or 3DxTxCh),
%   but just up to the 3rd dimention
%
% Usage:
%   Data3DxT = interp3DxT(Data3DxT, ...
%       XYZinit, XYZfinal, LastDim, method)
%
% Args:
%   Data3DxT: 3DxT matrix
%   XYZinit: resolution of input data
%   XYZend: target resolution of data
%   LastDim: last dimension
%   intmethod: interpolation method
%       ('linear', 'nearest', 'cubic', 'spline')
%       (default, 'linear')
%
% Notes:
% units in microns order: [X Y Z] which is:
%   (height, width, depth), XYZinit and XYZfinal are (1, 3) vectors

fprintf('Resampling data ... ')
Data3DxT = single(Data3DxT);
Dims = double(size(Data3DxT));
nDim = length(Dims);
DimEnd = (Dims(1:3)-1).*XYZinit;

% last dimention to interpolate
if ~exist('LastDim', 'var')
    LastDim = 3;
end

if ~exist('intmethod', 'var')
    intmethod = 'linear';
end

% Order of dimentions to resample
AxesOrder = [1 2 3 4 5; 2 1 3 4 5; 3 1 2 4 5]; %(height, width, depth)
AxesOrder2 = [1 2 3 4 5; 2 1 3 4 5; 2 3 1 4 5]; %(height, width, depth)

% chop dimentions in case data is just 3D
AxesOrder = AxesOrder(:, 1:nDim);
AxesOrder2 = AxesOrder2(:, 1:nDim);

% resample
for Dim_i = 1:LastDim
    
    % do linear interpolation of each dimention separately from 1 to 3
    Init_Axes = 0:XYZinit(Dim_i):DimEnd(Dim_i);
    New_Axes = 0:XYZend(Dim_i):DimEnd(Dim_i);
    
    if XYZinit(Dim_i) ~= XYZend(Dim_i)
        fprintf('*')
        Data3DxT = permute(interp1(Init_Axes, permute(Data3DxT, AxesOrder(Dim_i, :)), ...
            New_Axes, intmethod), AxesOrder2(Dim_i, :));
    end
    
    clear Init_Axes New_Axes
    
end

fprintf(' done')

end
