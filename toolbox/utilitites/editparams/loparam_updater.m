function input_struct = loparam_updater(...
    input_struct, input_var, sel_field)
% loparam_updater: function that updates fields of input structure, or
%   sub-fields of selected field from the same structure
% 
% Usage:
%   input_struct = loparam_updater(...
%       input_struct, input_var, sel_field)
%
% Args:
%   input_struct: input structure
%   input_var: structure with input_struct fields to be updated
%   sel_field: selected field name from input_struct to update sub-fields
% 
% Returns:
%   input_struct: updated structure
%
% Notes:

if exist('input_var', 'var') && ~isempty(input_var)
    if exist('sel_field', 'var') && ~isempty(sel_field)
        
        allfields = fields(input_struct.(sel_field));

        for ii = 1:length(allfields)
            if isfield(input_var, allfields{ii})
                input_struct.(sel_field).(allfields{ii}) = ...
                    input_var(1).(allfields{ii}); 
            end
        end
    else
        
        allfields = fields(input_struct);
        
        for ii = 1:length(allfields)
            if isfield(input_var, allfields{ii})
                input_struct.(allfields{ii}) = ...
                    input_var(1).(allfields{ii}); 
            end
        end        
    end   
end

end
