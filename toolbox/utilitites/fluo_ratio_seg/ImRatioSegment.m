function [med_mask, expr] = ...
    ImRatioSegment(x_, y_, iparams)

% matlab version of ben's python code: make_mask
% c1 == x_, c2 == y_
% main parameters

pmm = []; 
pmm.thrs_f = 25;
pmm.thrs_std = 1;
pmm.slope_thresh_1 = 4;
pmm.slope_thresh_2 = 10;
pmm.outlier_std_thresh = 10^3; % for my data it seems like this cutoff is too low, so increase this ths
pmm.step = 10;
%pmm.initial_window = 0.05; pmm.slide_step = 0.01;
pmm.initial_window = 0.7;
pmm.slide_step = 0.02;
pmm.type = 1;
pmm.maxn = 8*10^6;
pmm.minpixsiz = 10^2;
% blur sigma / kernel size
pmm.bgate = 0;
pmm.siz = 3;
pmm.sig = 2;
% median filtering
pmm.sd_sig = 3;
pmm.med_siz = 3;
% for plotting
pmm.verbose = false;
pmm.debug = false;
pmm.plotr = true;
pmm.deltaidx = 10;
pmm.deltabnoridx = 0.01;
% for saving
pmm.fname = [];
pmm.oDir = '.';

if ~exist('iparams', 'var'); iparams = []; end
pmm = loparam_updater(pmm, iparams);

% 1) preprocess images
% get input size
sizY = size(x_);

% homogenize class
x_ = single(x_);
y_ = single(y_);

% smooth image
if pmm.bgate
    x_ = imblur(x_, ...
        pmm.siz, pmm.sig, numel(sizY));
    y_ = imblur(y_, ...
        pmm.siz, pmm.sig, numel(sizY));
end

% flatten
x_ = x_(:);
y_ = y_(:);

% threshold intensity in absolute values
idx2del = x_ <= pmm.thrs_f;
x_(idx2del) = nan;
y_(idx2del) = nan;
clear idx2del

% normalize intensity
x_ = (x_ - min(x_(:)))/max(x_(:));
y_ = (y_ - min(y_(:)))/max(y_(:));

% hard threshold intensity in relative values
idx2del = x_ <= 0.005 | y_ <= 0.005;
x_(idx2del) = nan;
y_(idx2del) = nan;
clear idx2del

% remove major outliers
outliers = x_ > (mean(x_(:)) + pmm.outlier_std_thresh*std(x_(:)));
figH = [];

if pmm.verbose || pmm.plotr
    
    fprintf(['Outliers ', num2str(sum(outliers)), '\n'])
    figH = figure();
    axH(1) = subplot(1, 1, 1);
    
    % image params
    % density resolution
    wpar.width = 1000;
    wpar.height = 1000;
    
    % smoothing params
    wpar.sig = [2 2];
    wpar.siz = [6 6];
    
    densityplot(double(x_(:))', double(y_(:))', ...
        wpar.width, wpar.height, [], ...
        wpar.sig, wpar.siz, axH(1), 1)
    
    hold(axH(1), 'on'); 
    plot(x_(outliers), y_(outliers), 'r.', 'Parent', axH(1))
    
end

y_(outliers) = nan;
x_(outliers) = nan;

% fit line
if pmm.type
    
    % explore all y axes (come from top to bottom), aim is to cut the very
    % positive values
    i = 0;
    m = 1;
    m_all = [];
    i_all = [];
    
    mask = ~isnan(x_);
    mask_t = mask;
    win = pmm.initial_window;
    sstep = pmm.slide_step;
    
    while (1 - win - i*sstep) >= 0 || i == 100
        
        idx = (y_ >= 1 - win - i*sstep) & mask;
        [m, yint] = regression_int(x_, y_, idx, 0);
        
        % reset mask
        mask = mask_t;
        
        % update mask
        resid = (m*x_+ yint - y_);
        resid_sd = [std(resid(resid > 0)) std(resid(resid < 0))];
        mask(resid > resid_sd(2)*pmm.sd_sig) = 0;
        
        fprintf(['Iteration # ', num2str(i), ...
            ' slope ', num2str(m), ...
            ' std ', num2str(resid_sd), '\n']);
        
        if pmm.verbose
            [~, AxH] = densityplot(double(x_(idx)'), ...
                double(y_(idx)'), wpar.width, wpar.height, ...
                [], wpar.sig, wpar.siz, [], 1); 
            hold(AxH, 'on');
            plot(0:pmm.deltabnoridx:1, ...
                (0:pmm.deltabnoridx:1)*m + yint, ...
                'g', 'Linewidth', 2,  'Parent', AxH)
            drawnow;
        end
        
        if pmm.plotr
            hold(axH(1), 'on');
            plot(0:pmm.deltabnoridx:1, ...
                (0:pmm.deltabnoridx:1)*m + yint, ...
                'r', 'Linewidth', 2,  'Parent', axH(1))
            plot([0 1], [1 - win - i*sstep, 1 - win - i*sstep], ...
                'k', 'Linewidth', 2,  'Parent', axH(1))
            xlim(axH(1), [0 1]);
            ylim(axH(1), [0 1]);
            drawnow;
        end
        
        m_all(i + 1) = m;
        i_all(i + 1) = i;
        i = i + 1;
        
    end
    
    clear mask_t
    
    % conditions: toss the first segment and the ones below 1    
    [m, yint] = regression_int(x_, y_, mask, 0);
    
    fprintf(['Iteration # ', num2str(i), ...
        ' slope ', num2str(m), ...
        ' std ', num2str(resid_sd), '\n']);
    
    if pmm.plotr
        hold(axH(1), 'on');
        plot(0:pmm.deltabnoridx:1, ...
            (0:pmm.deltabnoridx:1)*m + yint, ...
            'g', 'Linewidth', 2,  'Parent', axH(1))
    end
    
else
    
    % generate initial mask and regression
    i = 0;
    m = 0;
    win = pmm.initial_window;
    
    while m < pmm.slope_thresh_1
        
        [~, m, yint] = ...
            regression(x_(x_ >= i*win & x_ < (i+1)*win)', ...
            y_(x_ >= i*win & x_ < (i+1)*win)');
        i = i + 1;
        
        if pmm.verbose
            fprintf(['Iteration # ', num2str(i), ...
                ' slope ', num2str(m), '\n']);
            hold(axH(1), 'on');
            
            plot(0:pmm.deltabnoridx:1, ...
                (0:pmm.deltabnoridx:1)*m + yint, ...
                'k', 'Linewidth', 2,  'Parent', axH(1))
            xlim(axH(1), [0 1]);
            ylim(axH(1), [0 1]);
            drawnow;
        end
        
    end
    
    mask = (x_ >= i*win);
    [m, yint, dif, pos, neg, ...
        res, pc_pos, pc_neg] = ...
        regression_int(x_, y_, mask, pmm.step);
    
    if pmm.verbose
        plot(0:pmm.deltabnoridx:1, ...
            (0:pmm.deltabnoridx:1)*m + yint, ...
            'g', 'Linewidth', 2,  'Parent', axH(1))
    end
    
    % regression iterations
    i = 0;
    k = true;
    
    while k
        
        last_mask = mask;
        last_par = [m, yint];
        
        if dif < 0 && m > pmm.slope_thresh_2
            % it mainly should depend on dif to pass 0, but if slope_thresh_2
            % is too low it will stop earlier
            k = False;
        end
        
        mask(res' > pc_pos) = false;
        [m, yint, dif, pos, neg, ...
            res, pc_pos, pc_neg] = ...
            regression_int(x_, y_, mask, pmm.step);
        i = i + 1;
        
        if pmm.verbose
            %plot(x_(~mask), y_(~mask), 'b.', 'Parent', Ax(2));
            plot(0:pmm.deltabnoridx:1, ...
                (0:pmm.deltabnoridx:1)*m + yint, ...
                'color', [1-(i/100) 0 0], 'Linewidth', 2,  'Parent', axH(2))
            fprintf(['Iteration # ', num2str(i), ...
                ' slope ', num2str(m), ...
                ' dif ', num2str(dif), '\n'])
            drawnow;
        end
        
    end
    
    mask = last_mask;
    m = last_par(1);
    yint = last_par(2);
    
end

% generate sd index
maskresid = m*x_(mask) + yint - y_(mask);
residstd = std(maskresid(maskresid < 0));
resid = (m*x_+ yint - y_);

% expression normalized to std of residuals in this brain
expr = resid./residstd;
expr(resid < 0 & isnan(expr)) = 0;
expr = reshape(expr, sizY);

% delete small rois
[expr_m] = split_im_into_cc(expr > 0, pmm.minpixsiz);
expr = expr.*expr_m;

% generate median mask
expr_l = expr;
expr_l = expr_l > pmm.thrs_std;
med_mask = zeros(size(expr));

if numel(sizY) == 3
    med_mask = medfilt3(expr_l, ...
        [pmm.med_siz pmm.med_siz pmm.med_siz]);
else
    med_mask = medfilt2(expr_l, ...
        [pmm.med_siz pmm.med_siz]);
end

% delete small rois
[med_mask] = split_im_into_cc(med_mask > 0, pmm.minpixsiz);

% save results
if pmm.plotr && ~isempty(pmm.fname)
    
    % save line fitting
    [filename, ~] = split_path(pmm.fname);
    saveas(figH, [pmm.oDir, filesep, filename, '.tif']);
    close(figH)
    
    % save video of segmentation
    pi.lag = 0;
    pi.vgate = 1;
    pi.vname = [pmm.oDir, filesep, filename, '_seg'];
    pi.range = [0 min([6, max(expr(:))])];
    pi.sizY = size(med_mask);
    pi.Y1 = med_mask; 
    slice3Dmatrix(expr, pi);
    close(gcf);
    clear pi
    
end

end

function [m, yint, pc_dif, pos, neg, res, ...
    pc_pos, pc_neg] = ...
    regression_int(x_, y_, mask, step, maxn)

if ~exist('maxn', 'var')
    maxn = 8*10^6;
end

x_ = x_(mask);
y_ = y_(mask);

if size(x_, 1) > 1
    x_ = x_';
end

if size(y_, 1) > 1
    y_ = y_';
end

idx2pick = randperm(numel(x_));
idx2pick = idx2pick(1:min(maxn, numel(x_)));
[~, m, yint] = regression(x_(idx2pick), y_(idx2pick));

res = m*x_ + yint - y_;
pos = res>0;
neg = res<0;
pc_pos = prctile(abs(res(pos)), 100-step);
pc_neg = prctile(abs(res(neg)), 100-step);
pc_dif = pc_pos - pc_neg;

end
