function [DataNew] = interp3DxTixTj(Data, ...
    initXY, endXY, iniT, endT, intmethod)
% interp3DxTixTj: resampling 2DxT or 3DxT volumes (or 3DxTxCh), resamples XYZ, and
% re-slices XY planes per XYZ volume using time per plane to align all
% planes per volume
%
% Usage:
%   Data3DxT = interp3DxT(Data, ...
%       XYZinit, XYZfinal, LastDim, method)
%
% Args:
%   Data: 2DxT or 3DxT matrix
%   initXY: resolution of input data
%   endXY: target resolution of data
%   iniT: input time per frame [z, n]
%   endT: output time per frame [z, n]
%   intmethod: interpolation method
%       ('linear', 'nearest', 'cubic', 'spline')
%       (default, 'linear')
%
% Notes:
% assumes Dim 3 is Z, and Dim 4 is time

fprintf('Resampling data ... ')
if ~exist('intmethod', 'var'); intmethod = 'linear'; end

dDim = size(Data);
tDim = numel(endT);

if length(dDim) > 3
    
    DataNew = zeros([dDim(1:3), tDim]);
    
    for z_i = 1:size(Data, 3)
        
        fprintf('*')
        if size(iniT, 1) == 1
            
            % assumes all planes were collected a the same time
            DataNew(:, :, z_i, :) = interp2DxTixTj(squeeze(Data(:, :, z_i, :)), ...
                initXY, endXY, iniT, endT, intmethod);
            
        else
            DataNew(:, :, z_i, :) = interp2DxTixTj(squeeze(Data(:, :, z_i, :)), ...
                initXY, endXY, iniT(z_i, :), endT, intmethod);
        end
        
    end
    
else
    
    DataNew(:, :, :) = interp2DxTixTj(squeeze(Data), ...
        initXY, endXY, iniT, endT, intmethod);
    
end

fprintf(' done\n')

end

function Data2DxT = interp2DxTixTj(Data2DxT, ...
    XYinit, XYend, Ti, Te, intmethod)
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
%   Ti: input time resolution
%   Te: target time resolution
%   method: interpolation method
%       ('linear', 'nearest', 'cubic', 'spline')
%       (default, 'linear')

Data2DxT = double(Data2DxT);
Dims = double(size(Data2DxT));
DimEnd = (Dims(1:2)-1).*XYinit;

% last dimention to interpolate
LastDim = 3;

% Order of dimentions to resample
AxesOrder = [1 2 3; 2 1 3; 3 1 2];
AxesOrder2 = [1 2 3; 2 1 3; 2 3 1];

% resample
for Dim_i = 1:LastDim
    
    Init_Axes = [];
    New_Axes = [];
    
    % do linear interpolation of each dimention separately from 1 to 3 (3rd is time)
    if Dim_i <= 2 && XYinit(Dim_i) ~= XYend(Dim_i)
        
        Init_Axes = 0:XYinit(Dim_i):DimEnd(Dim_i);
        New_Axes = 0:XYend(Dim_i):DimEnd(Dim_i);
        Data2DxT = permute(interp1(Init_Axes, permute(Data2DxT, AxesOrder(Dim_i, :)), ...
            New_Axes, intmethod), AxesOrder2(Dim_i,:));
        
    elseif Dim_i == 3
        
        Init_Axes = Ti;
        New_Axes = Te;
        Data2DxT = permute(interp1(Init_Axes, permute(Data2DxT, AxesOrder(Dim_i, :)), ...
            New_Axes, intmethod), AxesOrder2(Dim_i,:));
        
    end
    
    clear Init_Axes New_Axes
    
end

end
