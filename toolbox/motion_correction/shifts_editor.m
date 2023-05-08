function Out = shifts_editor(shifts_i, edit2use, ishifts, idiff)
% shifts_editor: read and write shifts from NoRMCorre
%
% Usage:
%   Out = shifts_editor(shifts_i, edit2use, ishifts)
%
% Args:
%   shifts_i: input shifts structure
%   edit2use: command to read or write
%   ishifts: input shifts as vector to replace structure
%   idiff: 
% 
% Notes:
%   beware shifts is a struct with fields: 'shifts', 'shifts_up' and 'diff'
%       for a linear or cubic interpolation (using apply_shifts) you only need 'shifts'
%       but for fft interpolation you also need 'shifts_up' and 'diff', so far not saved
%   see normcorre_batch, apply_shifts

if ~exist('idiff', 'var')
    idiff = [];
end

switch edit2use
    case 'read'
        
        % works when only the last dimension has the information
        shift_type = size(shifts_i(1).shifts);
        
        if sum(shift_type(1:3) == 1) == 3
            
            shift_format = numel(shift_type);
            Out = horzcat(shifts_i(:).shifts);
            
            if shift_format == 4
                Out = squeeze(Out);
                if numel(shifts_i) == 1; Out = Out'; end
            else
                Out = Out';
            end
            
        else
            
            ldDim = size(shifts_i(1).shifts, 4);
            shifts_r = cat(ndims(shifts_i(1).shifts)+1, shifts_i(:).shifts);
            shifts_r = reshape(shifts_r, [], ldDim, numel(shifts_i));
            Out{1} = squeeze(shifts_r(:, 1, :))';
            Out{2} = squeeze(shifts_r(:, 2, :))';
            
            if ldDim == 3
                Out{3} = squeeze(shifts_r(:, 3, :))';
            end
            
        end
        
    case 'write_r'
        
        siz = size(ishifts, 2);
        
        % only edits the shifts and shifts_up fields
        for i = 1:numel(shifts_i)
            shifts_i(i).shifts(1, 1, 1, 1:siz) = ishifts(i, :)';
            shifts_i(i).shifts_up(1, 1, 1, 1:siz) = ishifts(i, :)';
            if ~isempty(idiff)
                shifts_i(i).diff(1) = idiff(i);
            end
        end
        
        Out = shifts_i;
        
end

end
