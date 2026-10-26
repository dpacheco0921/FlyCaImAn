function meta = nrrd_homogenize_units(meta)
% nrrd_homogenize_units: homogenize spaceunit nomenclature
%
% Usage:
%   meta = nrrd_homogenize_units(meta)
%
% Args:
%   meta: nrrd metadata as output from nrrdread.m
%
% Returns:
%   meta: updated meta (spaceunits set to "microns" for all 3 axes; other units are
%       reported on screen and left unchanged)

if contains(meta.spaceunits, 'µm') || ...
        contains(meta.spaceunits, 'microns')
    meta.spaceunits = '"microns" "microns" "microns"';
else
    fprintf('found other units: %s\n', meta.spaceunits); % print as data, not as format
end

end