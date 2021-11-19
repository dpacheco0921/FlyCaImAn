function [sheet_1, sheet_2, overlapping_values, ...
    overlapping_idx_1, overlapping_idx_2] = ...
    intersect_gsheets(DOCID_1, GID_1, column_n1, ...
    DOCID_2, GID_2, column_n2, column2display)
% intersect_gsheets: intersect a particular column of two different sheets
%
% Usage:
%   [sheet_1, sheet_2, overlapping_values, ..
%       overlapping_idx_1, overlapping_idx_2] = ...
%       intersect_gsheets(DOCID_1, GID_1, column_n1, ...
%       DOCID_2, GID_2, column_n2)
%
% Args:
%   DOCID_1: see the value after 'key=' in your spreadsheet's url
%   GID_1: value after 'gid=' in url - optional
%   column_n1: column to compare from sheet 1
%   DOCID_2: see the value after 'key=' in your spreadsheet's url
%   GID_2: value after 'gid=' in url - optional
%   column_n2: column to compare from sheet 2
%   column2display: columnns to display for selecting duplicated
%
% Returns:
%   sheet_1, sheet_2: sheets downloaded from DOCID_1 and DOCID_2
%   overlapping_values: overlapping IDs/values between selected columns
%   overlapping_idx_1, overlapping_idx_2: indeces of overlapping IDs/values
%       for each sheet.
% 
% Notes:
% might be worth making DOCID_2 and GID_2 compatible with multiple sheets,
% so you compare a single one to many which would be the most frequent case
% with flywire sheets.

% download sheets and assign them to variables
sheet_1 = pull_gsheet(DOCID_1, GID_1);
column_1 = sheet_1(1:end, column_n1);

sheet_2 = cell(1, 1);
column_2 = cell(1, 1);

if iscell(DOCID_2) || iscell(GID_2)
    
    if iscell(DOCID_2) && iscell(GID_2)
        
        n_sheets = numel(DOCID_2);
        
        if size(column_n2, 2) ~= n_sheets
            column_n2 = repmat(column_n2, [1, n_sheets]);
        end
        
        for i = 1:n_sheets
        	sheet_2{i} = pull_gsheet(DOCID_2{i}, GID_2{i});
            column_2{i} = sheet_2{i}(1:end, column_n2(i));
        end
        
    elseif iscell(GID_2)
        
        n_sheets = numel(GID_2);
        if size(column_n2, 2) ~= n_sheets
            column_n2 = repmat(column_n2, [1, n_sheets]);
        end
        
        for i = 1:n_sheets
        	sheet_2{i} = pull_gsheet(DOCID_2, GID_2{i});
            column_2{i} = sheet_2{i}(1:end, column_n2(i));
        end
        
    else
        
        fprintf('Error')
        return
        
    end
    
else
    sheet_2{1} = pull_gsheet(DOCID_2, GID_2);
    column_2{1} = sheet_2{1}(1:end, column_n2);
end

% find intersection between selected columns
for i = 1:numel(column_2)
    [overlapping_values{i}, ~, overlapping_idx_2{i}] = intersect(column_1, column_2{i});
    overlapping_idx_1{i} = zeros(size(overlapping_values{i}));
end

% check if there is only one match
for j = 1:numel(column_2)
    for i = 1:numel(overlapping_values{j})

        if sum(contains(column_1, overlapping_values{j}(i))) ~= 1

            % get all indeces
            idx2use = find(contains(column_1, overlapping_values{j}(i)));

            display([column_1(idx2use), sheet_1(idx2use, column2display)])

            prompt = 'Select index';
            idx = input(prompt);

            overlapping_idx_1{j}(i, 1) = idx2use(idx);

        else

            overlapping_idx_1{j}(i, 1) = find(contains(column_1, overlapping_values{j}(i)));

        end

    end
end

end
