function [trace2keep_idx, xyz] = to_get_trace_with_nans(iObject, imatrix)
% to_get_trace_with_nans: identify traces woth consecutive nans
% Usage:
%   to_get_trace_with_nans(iObject)
%
% Args:
%   iObject: tree obj
%   imatrix: field to use
%
% Note:
% see getObj_xyz

if ~exist('imatrix', 'var'); imatrix = 0; end

% get xyz matrix
iObject = to_pick_xyz_object(iObject, [], imatrix);
xyz = getObj_xyz(iObject);

% by default works with groups of xyz
if ~iscell(xyz)
    xyz = {xyz};
end

xyz_n = numel(xyz);

trace2keep_idx = true(xyz_n, 1);

for i = 1:xyz_n
    
    fprintf('*')
    
    % find NaNs
    i_nan = isnan(xyz{i, 1});
    i_nan = sum(i_nan, 2) > 0;
    
    if sum(i_nan) > 0
        
        csc_nan = find(i_nan == 1);
        
         % more than one consecutive NaN
         if sum(diff(csc_nan) == 1) == 0
             
            for j = csc_nan'
                if j ~= 1 && j ~= numel(i_nan)
                    xyz{i, 1}(j, :) = ....
                        (xyz{i, 1}(j - 1, :) + xyz{i, 1}(j + 1, :))/2;
                end
            end
            
            if ~isempty(find(csc_nan == 1, 1))
               xyz{i, 1}(1, :) = []; 
            end
            
            if ~isempty(find(csc_nan == numel(i_nan), 1))
               xyz{i, 1}(end, :) = []; 
            end 
            
            trace2keep_idx(i) = true;
            
        else
            
            trace2keep_idx(i) = false;
            
        end
        
    end
    
end

end

