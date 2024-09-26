function [iObject, err_status] = ...
    to_atlas2atlas_xform(iObject, ...
    refi, refo, xDir, esuffix)
% to_atlas2atlas_xform: function that performs reformatting
%   of trees from atlas to atlas coordinates using CMTK or CPD
%
% Usage:
%   [iObject, err_status] = ...
%       to_atlas2atlas_xform(iObject, ...
%       refi, refo, xDir, esuffix)
%
% Args:
%   iObject: input object (tree or xyz matrix)
%   refi: input reference atlas
%   refo: output reference atlas
%   xDir: directory of transformation
%   esuffix: extra suffix

if ~exist('esuffix', 'var') || ...
        isempty(esuffix)
    esuffix = '';
end

if ~exist('xDir', 'var') || ...
        isempty(xDir)
    xDir = '.';
end

if contains(refi, 'nsybIVAi') || contains(refo, 'nsybIVAi')

    % load transformation
    if contains(refo, 'nsybIVAi')

        transform_ = 'reg_1to2_woOL_cbrain_lambda_1_beta_5_outliers_02_smp_1.mat';
        load([xDir, filesep, transform_], 'reg_struct')
        reg_struct = reg_struct.reg_1to2;

    elseif contains(refi, 'nsybIVAi')

        transform_ = 'reg_2to1_woOL_cbrain_lambda_1_beta_10_outliers_02_smp_1.mat';
        load([xDir, filesep, transform_], 'reg_struct')                 
        reg_struct = reg_struct.reg_2to1;

    end

    xyz = getObj_xyz(iObject);
    
    % invert axis (trasformation takes 3D matrices were first axis is y instead of x)
    xyz = xyz(:, [2 1 3]);
    
    if size(xyz, 1) > 10^4
        
        xyz = chunk2cell(xyz, 10^4);
        
        for i = 1:numel(xyz)
            stic
            xyz_temp{i, 1} = apply_tform(xyz{i}, reg_struct);
            stoc
        end
        
        xyz = cell2mat(xyz_temp);
        clear xyz_temp i
        
    else
        
        xyz = apply_tform(xyz, reg_struct);
        
    end
    
    clear reg_struct
    xyz = xyz(:, [2 1 3]);
    iObject = updObj_xyz(iObject, xyz);
    err_status = 0;

else

    ibpars.xDir = xDir;
    ibpars.esuffix = ['*', esuffix];
    [~, ~, repoDirs] = customdirs_deafult;
    ibpars.xDir = repoDirs{3};

    [iObject, err_status] = cmtk_xform(iObject, refi, refo, ibpars);
    % edit cmtk_streamxform
    % edit cmtk_xform

end

end
