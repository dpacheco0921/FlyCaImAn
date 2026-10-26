function [trace2keep_idx, xyz] = to_get_trace_with_nans(iObject, imatrix)
% to_get_trace_with_nans: identify traces woth consecutive nans
% Usage:
%   [trace2keep_idx, xyz] = to_get_trace_with_nans(iObject, imatrix)
%
% Args:
%   iObject: tree obj
%   imatrix: field to use (0 xyz, 1 xyz_s, 2 xyz_pix, 3 traces) (default, 0)
%
% Returns:
%   trace2keep_idx: [trees by 1] logical, false for traces that have
%       consecutive NaNs (cannot be repaired)
%   xyz: cell with one [points by 3] matrix per tree; isolated NaN points
%       are replaced by the mean of their two neighbours (NaNs at the first
%       or last point are removed)
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

