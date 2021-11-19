function plot_plane(handles)
% function that plots each plane

global p

% get variables
p.z = round(get(handles.slider1, 'Value'));
p.cha = round(get(handles.slider2, 'Value'));
p.min_max(p.cha, :) = [str2num(get(handles.edit1, 'String')), ...
    str2num(get(handles.edit2, 'String'))];

% plot plane
imagesc(p.cha_matrix(:, :, p.z, p.cha), 'parent', handles.axes1)
handles.axes1.XTick = [];
handles.axes1.YTick = [];

caxis(handles.axes1, p.min_max(p.cha, :))

% show channel name
handles.axes1.Title.String = [p.chaname{p.cha} ' z-' num2str(p.z)];

% plot ROIs in this plane
if ~isempty(p.roi_center)
    
    roi2use = find(p.roi_center(:, 3) == p.z);
    
    roi_center = p.roi_center(roi2use, :);
    roi_idx = p.roi_idx(roi2use, :);
    roi_radious = p.roi_radious(roi2use, :);
    
    for i = 1:numel(roi_idx)
        
        text(roi_center(i, 1)-3, roi_center(i, 2)-3, ...
            num2str(roi_idx(i)), 'color', [1 0 0],...
            'FontSize', max(roi_radious(i)/4, 14))

    end
    
end

end
