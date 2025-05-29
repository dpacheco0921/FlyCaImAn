function rho = neighcorr_3D(Y, iparams)
% neighcorr_3D: Function that uses correlation_image_3D
%   to generate a correlation image based on neighboring 
%   pixels, but in this case it deals better with big
%   matrices (memmap).
%
% Usage:
%   rho = neighcorr_3D(Y, iparams)
%
% Args:
%   Y: M x N x O x T movie or memmap object (looks for Y.Y)
%   iparams: parameter variable
%
% Returns:
%   rho: M x N x O matrix, cross-correlation with adjacent pixel
%
% Notes: 
% splits volume in n parts and runs correlation_image_3D on each

% deafult params
p.overlap = [2 2 2];
p.timestamps = [];
p.chunksize = [50 50 10];

if ~exist('iparams', 'var'); iparams = []; end
p = loparam_updater(p, iparams);

% get file info
memmaped = isobject(Y);

if memmaped
    
    sizY = Y.sizY; 
    if ~isfield(p, 'chunksize') || isempty(p.chunksize)
        p.chunksize = [sizY(1:2) 10]; 
    end
    patches = construct_patches(sizY(1:end-1), p.chunksize, p.overlap);
    
else
    sizY = size(Y);
end

if memmaped
    
    % Do partial volumes
    fprintf('Running neighcorr_3D for big volumes, patch number : ')
    fprintf([num2str(length(patches)), ' '])
    
    for p_i = 1:length(patches)
        
        temp_rho{p_i} = lc3Dgen(Y.Y(patches{p_i}(1):patches{p_i}(2), ...
            patches{p_i}(3):patches{p_i}(4), ...
            patches{p_i}(5):patches{p_i}(6), :), p.timestamps);
        planeidx(p_i, :) = [patches{p_i}(5) patches{p_i}(6)];
        yidx(p_i, :) = [patches{p_i}(1) patches{p_i}(2)];
        xidx(p_i, :) = [patches{p_i}(3) patches{p_i}(4)];
        fprintf('*')
        
    end
    
    fprintf(' ')
    
    % Correct interseccions
    planenum = sort(unique(planeidx(:, 1))); 
    ynum = sort(unique(yidx(:, 1))); 
    xnum = sort(unique(xidx(:, 1))); 
    rho = zeros(sizY(1:end-1));
    
    for p_ii = 1:length(patches)
        
        % Prunning edges (depth)
        if planeidx(p_ii, 1) == planenum(1)
            
            if numel(planenum) == 1
                zi = planeidx(p_ii, 1):planeidx(p_ii, 2);
                zii = 1:planeidx(p_ii, 2);
            else
                zi = planeidx(p_ii, 1):(planeidx(p_ii, 2) - 1);
                zii = 1:(planeidx(p_ii, 2)-planeidx(p_ii, 1));                
            end
            
        elseif planeidx(p_ii, 1) == planenum(end)
            
            zi = (planeidx(p_ii, 1)+1):(planeidx(p_ii, 2));
            zii = 2:(planeidx(p_ii, 2)-planeidx(p_ii, 1) + 1);
            
        else
            
            zi = (planeidx(p_ii, 1)+1):(planeidx(p_ii, 2)-1);
            zii = 2:(planeidx(p_ii, 2)-planeidx(p_ii, 1));
            
        end
        
        if yidx(p_ii, 1) == ynum(1)
            
            if numel(ynum) == 1
                yi = yidx(p_ii, 1):yidx(p_ii, 2);
                yii = 1:yidx(p_ii, 2);
            else
                yi = yidx(p_ii, 1):(yidx(p_ii, 2) - 1);
                yii = 1:(yidx(p_ii, 2)-yidx(p_ii, 1));
            end
            
        elseif yidx(p_ii, 1) == ynum(end)
            
            yi = (yidx(p_ii, 1)+1):(yidx(p_ii, 2));
            yii = 2:(yidx(p_ii, 2)-yidx(p_ii, 1) + 1);
            
        else
            
            yi = (yidx(p_ii, 1)+1):(yidx(p_ii, 2)-1);
            yii = 2:(yidx(p_ii, 2)-yidx(p_ii, 1));
            
        end

        if xidx(p_ii, 1) == xnum(1)
            
            if numel(xnum) == 1
                xi = xidx(p_ii, 1):xidx(p_ii, 2);
                xii = 1:xidx(p_ii, 2);
            else            
                xi = xidx(p_ii, 1):(xidx(p_ii, 2) - 1);
                xii = 1:(xidx(p_ii, 2)-xidx(p_ii, 1));
            end
            
        elseif xidx(p_ii, 1) == xnum(end)
            
            xi = (xidx(p_ii, 1)+1):(xidx(p_ii, 2));
            xii = 2:(xidx(p_ii, 2)-xidx(p_ii, 1) + 1);
            
        else
            
            xi = (xidx(p_ii, 1)+1):(xidx(p_ii, 2)-1);
            xii = 2:(xidx(p_ii, 2)-xidx(p_ii, 1));
            
        end
        
        rho(yi, xi, zi) = temp_rho{p_ii}(yii, xii, zii);
        fprintf('+')
        
    end
    
else
    
    rho = lc3Dgen(Y);
    
end

fprintf('\n')

end

function rho = lc3Dgen(data, timestamps)
% lc3Dgen: get local correlation (26 neighboring pixels)
%
%
% Usage:
%   rho = lc3Dgen(data, timestamps)
%
% Args:
%   data: M x N x T movie
%   timestamps: timestamps to use

xdims = size(data);

if exist('timestamps', 'var') && ~isempty(timestamps)
    data = data(:, :, :, timestamps);
end

rho = correlation_image_3D(data, 26, xdims(1:3));

end
