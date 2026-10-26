% copy, edit and save this function as 'fiji_fullpath.m'
% provide the right fiji directory per system (PC, mac)

function fijiDirs = fiji_fullpath
% fiji_fullpath: user defined directory of FIJI executable (full path) 
%
% Usage:
%   fijiDirs = fiji_fullpath
%
% Args: none (the path is hardcoded below, per computer)
%
% Returns:
%   fijiDirs: directory of FIJI executable (full path)

% fijiDir

% returns [] on other systems (e.g. Linux), so callers must check for empty
fijiDirs = [];

if ispc
    fijiDirs = 'C:\Users\*\Fiji.app\ImageJ-win64.exe';
elseif ismac
    fijiDirs = '/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx';
end

end
