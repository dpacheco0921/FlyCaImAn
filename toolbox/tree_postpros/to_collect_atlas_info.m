function to_collect_atlas_info(obj, atlasname, aDir, deviceId)
% to_collect_atlas_info: code to collect atlas information and save it in a mat file
% Usage:
%   to_collect_atlas_info(obj, atlases, aDir, deviceId)
%
% Args:
%   obj: trace obj
%   atlases: name of atlases
%   aDir: atlas directory
%   deviceId: device id

if ~exist('deviceId', 'var') || isempty(deviceId); deviceId = 1; end
if ~exist('aDir', 'var') || isempty(aDir); aDir = customdirs_deafult(deviceId); end
if ~exist('atlasname', 'var') || isempty(atlasname)
    atlasname = {'JFRC2', 'VNCIS1', 'IS2', 'FCWB', 'IBNWB', 'nsybIVAi'};
end

atlas = [];

for i = 1:numel(atlasname)
    atlas(i).name = atlasname{i};
    [im, meta] = nrrdread([aDir, filesep, atlasname{i}, '.nrrd']);
    atlas(i).voxsiz = nrrdread_res(meta);
    atlas(i).siz = size(permute(im, [2 1 3]));
    atlas(i).boundingbox = [0 0 0; atlas(i).voxsiz.*(atlas(i).siz - 1)];
end

save([strrep(which('traceObj'), 'traceObj.m', ''), 'atlas_meta.mat'], 'atlas')
end