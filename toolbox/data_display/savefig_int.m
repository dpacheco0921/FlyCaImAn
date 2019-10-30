function savefig_int(figH, tDir, ...
    figname, figformat, resolution_)
% savefig_int: save figure in various formats
% 
% Usage:
%   savefig_int(figH, tDir, ...
%       figname, figformat, resolution_)
%
% Args:
%   figH: figure handle
%   tDir: target directory
%   figname: figure name
%   figformat: which formats to save (R = [1x12])
%       (1, fig, 2, tif)
%       (3, eps (build-in), 4, eps-ogl (build-in))
%       (5, expfig_eps_nocrop, 6 expfig_eps_crop)
%       (7, expfig_eps-ogl_nocrop, 8 expfig_eps-ogl_crop)
%       (9, png (build-in), 10, expfig_png)
%       (11, svg (build-in), 12 svg)
%   resolution_: resolution (default -r300), in dots per pixel
% 
% Notes:

if strcmp(filesep, tDir(end))
    tDir = tDir(1:end-1);
end

if ~exist('figformat', 'var') || isempty(figformat)
    figformat = zeros(1, 12);
end

if length(figformat) < 12
    figformat(12) = 0;
end

if ~exist('resolution_', 'var') || isempty(resolution_)
    resolution_ = '-r300';
end

if sum(figformat)
    
    if figformat(1)
        saveas(figH, [tDir, filesep, figname, '.fig']);
    end
    
    if figformat(2)
        saveas(figH, [tDir, filesep, figname, '.tif']);
    end
    
    % ********** EPS-build-in **********
    %   this is necessary for overlays or figures with transparency
    if figformat(3)
        print(figH, [tDir, filesep, figname, '_pt'], ...
            '-depsc2', resolution_);
    end
    
    if figformat(4)
        print(figH, [tDir, filesep, figname, '_ogl'], ...
            '-opengl', '-depsc2', resolution_);
    end
    
    % ********** EPS-export_fig **********
    if figformat(5)
        export_fig(figH, [tDir, filesep, figname, '_fept'], ...
            '-eps', '-painters', '-q100', '-transparent', '-nocrop');
    end
    
    if figformat(6)
        export_fig(figH, [tDir, filesep, figname, '_fept'], ...
            '-eps', '-painters', '-q100', '-transparent');
    end
    
    if figformat(7)
        export_fig(figH, [tDir, filesep, figname, '_feogl'], ...
            '-eps', '-opengl', '-q100', '-transparent', '-nocrop');
    end
    
    if figformat(8)
        export_fig(figH, [tDir, filesep, figname, '_feogl'], ...
            '-eps', '-opengl', '-q100', '-transparent');
    end
    
    % ********** PNG-build-in **********
    if figformat(9)
        print(figH, [tDir, filesep, figname, '.png'], ...
            '-dpng', '-opengl', resolution_);
    end
    
    % ********** PNG-export_fig **********
    if figformat(10)
        export_fig(figH, [tDir, filesep, figname, '_fept'], ...
            '-png', '-painters', '-q100', '-transparent', '-nocrop');
    end
    
    % ********** SVG-build-in **********
    if figformat(11)
        print(figH, [tDir, filesep, figname, '.svg'], ...
            '-dsvg', '-painters');
    end
    
    % ********** using plot2svg package: deprecated **********
    % if figformat(12)
    %     plot2svg([tDir, filesep, figname, '_ec.svg'], figH);
    % end
    
else
    
    fprintf('not saving figure\n')
    
end

end
