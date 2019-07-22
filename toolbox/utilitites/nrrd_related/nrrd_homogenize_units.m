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
%   meta: updated meta

if contains(meta.spaceunits, 'µm') || ...
        contains(meta.spaceunits, 'microns')
    meta.spaceunits = '"microns" "microns" "microns"';
else
    fprintf(['found other units', meta.spaceunits, '\n'])
end

end