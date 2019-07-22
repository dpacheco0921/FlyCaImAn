function Y = framegapfill(tp2fill, Y, ftype, fside)
% framegapfill: fills frames or volumes from input timepoints
%   with the mean of the surrounding frames or volumes
%
% Usage:
%   Y = framegapfill(tp2fill, Y, ftype, fside)
%
% Args:
%   tp2fill: indexes of frames or volumes to fill
%   Y: original 3D or 4D matrix
%   ftype: type of replacing ('mean', default, or other (direct replacement))
%   fside: side to use (right or left to the timepoint to replace)

if ~exist('ftype', 'var'); ftype = 'mean'; end
if ~exist('fside', 'var'); fside = 'r'; end

if ~isempty(tp2fill)
    
    siz = size(Y); 
    tEnd = 4; % default
    
    if length(siz) == 4
        
        str2run{1} = 'Y(:, :, :, t_idx)';
        str2run{2} = 'Y(:, :, :, tNext)';
        dstr = num2str(tEnd);
        
    elseif length(siz) == 5
        
        str2run{1} = 'Y(:, :, :, t_idx, :)';
        str2run{2} = 'Y(:, :, :, tNext, :)';
        dstr = num2str(tEnd);
        
    else
        
        str2run{1} = 'Y(:, :, t_idx)'; 
        str2run{2} = 'Y(:, :, tNext)';
        dstr = num2str(length(siz));
        tEnd = 3;
        
    end
    
    fprintf(' filling empty frames:\n')
    
    if size(tp2fill, 1) > 1; tp2fill = tp2fill'; end
    
    for t_idx = tp2fill
        
        fprintf([' ', num2str(t_idx)])
        
        if t_idx == 1 || t_idx == size(Y, tEnd)
            
            % using preceding timepoint
            if t_idx == 1
                tNext = volsearchForward(tp2fill, t_idx);
            else
                tNext = volsearchBackward(tp2fill, t_idx);
            end
            
            eval([str2run{1},' = ', str2run{2}, ';'])
            
        else 
            
            % using following timepoint
            tNext(1, 1) = volsearchBackward(tp2fill, t_idx);
            tNext(1, 2) = volsearchForward(tp2fill, t_idx);
            
            if tNext(1, 1) < 1 || ...
                    (~contains(ftype, 'mean') && contains(fside, 'l'))
                
                fprintf(' l');
                tNext = tNext(1, 2);
                eval([str2run{1},' = ', str2run{2}, ';'])  
                
            elseif tNext(1, 2) > size(Y, tEnd) || ...
                    (~contains(ftype, 'mean') && contains(fside, 'r'))
                
                fprintf(' r');
                tNext = tNext(1, 1);
                eval([str2run{1},' = ', str2run{2}, ';'])
                
            else
                
                fprintf(' b');
                eval([str2run{1},' = mean(', str2run{2}, ', ', dstr,');'])
                
            end
            
        end
        
        fprintf(['(', num2str(tNext), ')\n'])
        
    end
    
    fprintf(')\n')
    
end

end

function t_idx = volsearchForward(vol2fill, t_idx)
% volsearchForward: Search for closest non-empty 
%   frame to average to (forward)
%
% Usage:
%   t_idx = volsearchForward(vol2fill, t_idx)
%
% Args:
%   vol2fill: 2D or 3D matrix
%   t_idx: starting time point

emptyf = 1;

while emptyf == 1
    
    t_idx = t_idx + 1;
    if sum(vol2fill == t_idx) == 0; emptyf = 0; end
    
end

end

function t_idx = volsearchBackward(vol2fill, t_idx)
% volsearchBackward: Search for closest non-empty
%   frame to average to (backwards)
%
% Usage:
%   t_idx = volsearchBackward(vol2fill, t_idx)
%
% Args:
%   vol2fill: 2D or 3D matrix
%   t_idx: starting time point

emptyf = 1;

while emptyf == 1
    
    t_idx = t_idx - 1;
    if sum(vol2fill == t_idx) == 0; emptyf = 0; end
    
end

end
