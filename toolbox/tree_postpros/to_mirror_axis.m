function iObject = to_mirror_axis(iObject, axis2flip, boundingbox)
% to_mirror_axis: mirror iObject (xyz coordinates within it)
%
% Usage:
%   to_mirror_axis(iObject, flip_axis, boundingbox)
%
% Args:
%   iObject: tree obj
%   axis2flip: axis to flip
%   boundingbox: bounding box

if ~iscell(axis2flip)
    axis2flip = {axis2flip};
end

boundingbox_int = zeros(2, 3);
flip_axis = ones(1, 3);

axis_ = {'X', 'Y', 'Z'};
idx = contains(axis_, axis2flip);
flip_axis(idx) = -1;

boundingbox_int(:, idx) = boundingbox(:, idx);

if iscell(iObject)
    
    for i = 1:numel(iObject)
        fprintf('*')
        iObject{i} = mirror_funct(iObject{i}, ...
            flip_axis, boundingbox_int);
    end
    
    fprintf('\n')
    
elseif isobject(iObject)
    % it is an array of traces
    
    for i = 1:numel(iObject.traces)
        fprintf('*')
        iObject.traces(i) = mirror_funct(iObject.traces(i), ...
            flip_axis, boundingbox_int);
    end
    
    fprintf('\n')    
    
elseif isstruct(iObject)

    for i = 1:numel(iObject)
        fprintf('*')
        iObject(i) = mirror_funct(iObject(i), ...
            flip_axis, boundingbox_int);
    end
    
    fprintf('\n')    
else
    
    iObject = mirror_funct(iObject, ...
        flip_axis, boundingbox_int);
    
end

end

function iObject = mirror_funct(...
    iObject, axis2flip, boundingbox)
% mirror_funct: mirror iObject (xyz coordinates within it)
%
% Usage:
%   iObject = mirror_funct(...
%       iObject, flip_axis, boundingbox)
%
% Args:
%   iObject: tree obj
%   axis2flip: axis to flip
%   boundingbox: bounding box

% make it compatible with any input type (need to update tree)
xyz = getObj_xyz(iObject);

% mirror
xyz = xyz.*(axis2flip) - boundingbox(1, :) ...
    + boundingbox(2, :);

% update xyz
iObject = updObj_xyz(iObject, xyz);

end
