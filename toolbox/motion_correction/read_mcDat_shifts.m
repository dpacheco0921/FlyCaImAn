function motion_out = read_mcDat_shifts(motion_cell, type)
% read_mcDat_shifts: function that reads mcDat shifts provided by NoRMCOrre
% 
% Usage:
%   read_mcDat_shifts(motion_cell, type)
%
% Args:
%   motion_cell: input shifts (cell)
%   type: output type
%       (cumulative (1) or all as 3D matrix (2))

if ~exist('type', 'var'); type = 1; end

if type == 1

    % provide cumulative motion across iterations
    motion_out = zeros(size(motion_cell{1}));
    
    for i = 1:numel(motion_cell)
        motion_out = motion_cell{i} + motion_out;
    end 
    
elseif type == 2
    
    % provide motion of each iteration in the 3rd dim
    motion_out = zeros([size(motion_cell{1}), ...
        numel(motion_cell)]);
    
    for i = 1:numel(motion_cell)
        motion_out(:, :, i) = motion_cell{i};
    end
    
end

end