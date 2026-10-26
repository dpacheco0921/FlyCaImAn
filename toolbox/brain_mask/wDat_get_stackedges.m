function lnbis = wDat_get_stackedges(wDat, ...
    planes, planeorder)
% wDat_get_stackedges: uses input 'planes' to 
%   define stack edges along the z axis of wDat.vSize
%
% Usage:
%   lnbis = wDat_get_stackedges(wDat, ...
%       planes, planeorder)
%
% Args:
%   wDat: main metadata structure
%   planes: planes that define stack edges
%       (default, [], one stack covering the whole volume)
%   planeorder: direction bottom2up or up2bottom
%       (default, 'bottom2up')
%
% Returns:
%   lnbis: 1 x (number of stacks + 1) row vector of z plane indices that
%       delimit the stacks (stack i spans planes lnbis(i)+1 to lnbis(i+1))

if ~exist('planeorder', 'var')
    planeorder = 'bottom2up';
end

if ~exist('planes', 'var')
    planes = [];
end

switch planeorder
    
    case 'up2bottom'
        
        lnbis(1, 1) = 0;
        
        for i = 1:numel(planes)
            lnbis(1, i + 1) = ...
                find(wDat.Zstitch.Zidx > planes(i), 1 , 'first');
        end
        
        lnbis(1, end + 1) = wDat.vSize(3);
        
    case 'bottom2up'
        
        lnbis(1, 1) = wDat.vSize(3);
        
        for i = 1:numel(planes)
            lnbis(1, i + 1) = ...
                find(wDat.Zstitch.Zidx > planes(i), 1 , 'last');
        end
        
        lnbis(1, end + 1) = 0;
        lnbis = flip(lnbis, 2);
        
end

end
