function status = pushrows_gsheets(DOCID, GID, ...
    sheet_template, sheet_edits, columns2update, ...
    row_idx, header_str)
% pushrows_gsheets: push rows to be updated in DOCID
%
% Usage:
%   status = pushrows_gsheets(DOCID, GID, ...
%       sheet_template, sheet_edits, columns2update, ...
%       row_idx, header_str)
%
% Args:
%   DOCID: see the value after 'key=' in your spreadsheet's url
%   GID: value after 'gid=' in url - optional
%   sheet_template: a nxm cell with a template
%   sheet_edits: a nxm cell with edits to make to template
%   columns2update: columnns to update
%   row_idx: indeces of rows in DOCID to update
%   header_str: optional header update

% download sheets and assign them to variables

% make sure sheet_templatehas the right size
if size(sheet_template, 2) < max(columns2update(2, :))
    sheet_template(:, max(columns2update(2, :))) = cell(size(sheet_template, 1), 1);
end

% update header
if exist('header_str', 'var') && ~isempty(header_str)
    mat2sheets(DOCID, GID, [1 1], header_str);
end

fprintf(['Updating ', num2str(size(sheet_template, 1)), ' rows'])

% make replacements
sheet_template(:, columns2update(2, :)) = sheet_edits(:, columns2update(1, :));
display(sheet_template)

for i = 1:numel(row_idx)
    fprintf('*')
    status(i, :) = mat2sheets(DOCID, GID, [row_idx(i) 1], sheet_template(i, :));
end
fprintf('... done\n')

end
