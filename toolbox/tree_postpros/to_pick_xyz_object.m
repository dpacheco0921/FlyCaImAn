function iObject = to_pick_xyz_object(obj, idx2sel, imatrix)
% to_pick_xyz_object: function to parse xyz field to use
%
% Usage:
%   to_pick_xyz_object(iObject, idx2sel, imatrix)
%
% Args:
%   obj: traceObj obj
%   idx2sel: files/indeces to use
%   imatrix: field to use

if ~exist('imatrix', 'var') || ...
        isempty(imatrix)
    imatrix = 0;
end

if ~exist('idx2sel', 'var') || ...
        isempty(idx2sel)

    if imatrix == 0
        iObject = obj.xyz;
    elseif imatrix == 1
        iObject = obj.xyz_s;
    elseif imatrix == 2
        iObject = obj.xyz_pix;
    elseif imatrix == 3
        iObject = obj.traces;
    end
    
else
    
    if imatrix == 0
        iObject = obj.xyz(idx2sel, 1);
    elseif imatrix == 1
        iObject = obj.xyz_s(idx2sel, 1);
    elseif imatrix == 2
        iObject = obj.xyz_pix(idx2sel, 1);
    elseif imatrix == 3
        iObject = obj.traces(idx2sel);
    end
    
end

end
