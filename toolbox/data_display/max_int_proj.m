function [figH, axH] = max_int_proj(Y, iparams)
% max_int_proj: basic maximun intensity projection function
% 
% Usage:
%   max_int_proj(Y, iparams)
%
% Args:
%   Y: 3D matrix or a cell array of 3D matrices
%   iparams: basic figure and axis parameters
% 
% Returns:
%   figH: figure handle
%   axH: axes handle

mippars.figH = [];
mippars.axH = [];
mippars.figpos = [];
mippars.irange = [];

if ~exist('iparams', 'var'); iparams = []; end
mippars = loparam_updater(mippars, iparams);

if ~iscell(Y); Y = {Y}; end

for i = 1:numel(Y)
    
    if isempty(mippars.figpos)
        mippars.figpos = genfigpos(1, 'center', [1306 443]);
    end
    
    if isempty(mippars.figH)
        mippars.figH = figure('position', mippars.figpos);
    end
    
    if isempty(mippars.axH)
        for j = 1:3
            if j == 1
                mippars.axH = subplot(1, 3, j);
            else
                mippars.axH(j) = subplot(1, 3, j);
            end
        end
    end
    
    if isempty(mippars.irange)
        
        mippars.irange(i, :) = [min(Y{i}(:)) max(Y{i}(:))];
        
    else
        
        if size(mippars.irange, 1) == 1
            mippars.irange = repmat(mippars.irange, [3, 1]);
        end
        
    end
    
    % plot projections in different subplots
    imagesc(max(Y{i}, [], 3), 'Parent', mippars.axH(1));
    caxis(mippars.axH(1), mippars.irange(1, :));
    mippars.axH(1).XLabel.String = 'X';
    mippars.axH(1).YLabel.String = 'Y';
    
    imagesc(flip(squeeze(max(Y{i}, [], 2))', 2), ...
        'Parent', mippars.axH(2));
    caxis(mippars.axH(2), mippars.irange(1, :));
    mippars.axH(2).XLabel.String = 'Y';
    mippars.axH(2).YLabel.String = 'Z';
    
    imagesc(flip(squeeze(max(Y{i}, [], 1))', 2), ...
        'Parent', mippars.axH(3));
    caxis(mippars.axH(3), mippars.irange(1, :));
    mippars.axH(3).XLabel.String = 'X';
    mippars.axH(3).YLabel.String = 'Z';
    
    figEdit(mippars.axH, mippars.figH)
    
end

axH = mippars.axH;
figH = mippars.figH;

end