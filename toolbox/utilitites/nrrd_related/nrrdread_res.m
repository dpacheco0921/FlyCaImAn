function iXYZres = nrrdread_res(meta)
% nrrdread_res: read pixel image resolution from nrrd metadata
%
% Usage:
%   iXYZres = nrrdread_res(meta)
%
% Args:
%   meta: nrrd metadata as output from nrrdread.m
%
% Returns:
%   iXYZres: 1 x 3 pixel size (x y z, i.e. width height depth), rounded to 1e-4 of the
%       metadata units (microns)

fres = 10^4; % round to 1e-4 (avoids floating-point noise in the metadata)
init = strfind(meta.spacedirections, '(');
endt = strfind(meta.spacedirections, ')');

for i_axes = 1:numel(init)
    string_i = strsplit2(meta.spacedirections(init(i_axes)+1:endt(i_axes)-1), ',');
    iXYZres(1, i_axes) = round(str2double(string_i{i_axes})*fres)/fres;
end

end