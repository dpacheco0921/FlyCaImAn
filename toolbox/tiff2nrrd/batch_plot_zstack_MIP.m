function batch_plot_zstack_MIP(...
    Filename, oDir, iparams)
% batch_plot_zstack_MIP: plot MIP of zstack
%
% Usage:
%   batch_plot_zstack_MIP(...
%       Filename, oDir, iparams)
%
% Args:
%   Filename: name pattern of files to use
%       (default, [])
%   oDir: output directory.
%   iparams: parameters to update
%       (fsuffix: suffix of metadata file)
%           (default, '_metadata.mat')
%       (dir_depth: depth of directory search)
%           (default, 0)
% Notes

% default params
metpars = [];
metpars.fsuffix = '_Zstack.nrrd';
metpars.dir_depth = 0;
metpars.nchannels = 2;

if ~exist('Filename', 'var')
    Filename = [];
end

if ~exist('oDir', 'var') || isempty(oDir)
    oDir = [pwd, filesep, 'Zstack_MIP'];
end

% update variables
if ~exist('iparams', 'var'); iparams = []; end
metpars = loparam_updater(metpars, iparams);

if ~isempty(oDir) && ~exist(oDir, 'dir')
    mkdir(oDir)
end

% define files to use
if metpars.dir_depth == 0
    f2run = rdir(['.', filesep, '*', ...
        metpars.fsuffix]);
elseif metpars.dir_depth == 1
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', metpars.fsuffix]);
else
    f2run = rdir(['.', filesep, '*', ...
        filesep, '*', filesep, '*', ...
        metpars.fsuffix]);
end

f2run = {f2run.name}';
[filename, iDir] = split_path(f2run);

if ~isempty(Filename)
    f2run = find(contains(filename, Filename));
    filename = filename(f2run);
    iDir = iDir(f2run);
end

fprintf(['Generating plots for ', ...
    num2str(numel(filename)), ' files\n'])

% plot sMod results
for i = 1:numel(filename)

    % load nrrd
    [Data, ~] = nrrdread([iDir{i}, filesep, filename{i}]);
    siz = size(Data);
    Data = double(reshape(Data, ...
        [siz(1:2), metpars.nchannels, prod(siz(3:end))/(metpars.nchannels)]));
    Data = permute(Data, [1 2 4 3]);
    filename_ = strrep(filename{i}, '_Zstack.nrrd', '');
    
    % plot MIP
    [figH, axH] = max_int_proj(Data(:, :, :, 1));
    [figH_ , axH_ ] = max_int_proj(Data(:, :, :, 2));
    
    figEdit(axH, figH, [], [], [], 'off')
    savefig_int(figH, oDir, [filename_, '_r_MIP'], ...
        [0 0 0 0 0 0 0 0 1])
    close(figH)
    
    figEdit(axH_, figH_, [], [], [], 'off')
    savefig_int(figH_, oDir, [filename_, '_g_MIP'], ...
        [0 0 0 0 0 0 0 0 1])
    close(figH_)

end

fprintf('... Done\n')

end
