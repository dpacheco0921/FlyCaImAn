%% Debugging LED denoising
%% 1) Difference between the two distributions of max fluorescence is not clear
iDat.RedChaMean = median(RedCha(:, :, :), 3);

%% 2) Difference between the two distributions of max fluorescence was not properly calculated
%% 2.1) Trying other gauss-number or selecting the right gauss-2use
gaussnumber = 'gauss3'; % default gauss2, so try lower or higher
sigmaRed = 5;
gauss2use = 2; % default 1, so if using more gaussians, try higher values
hDat = fitgauss1D(iDat.PMT_fscore(:), gaussnumber, sigmaRed, gauss2use);
fitgaussPlotter(hDat)

remotegate = 1;
if remotegate
    fhandle = gcf;
    fhandle.Position = [680 558 560 420];
end

%% 2.2) Select volumes to average (based on hDat.xths and iDat.PMT_fscore)
if params.pgate; fitgaussPlotter(hDat); end
redchaths = hDat.xths(1, 2);
if size(iDat.PMT_fscore, 1) == 1
    iDat.PMT_fscore = iDat.PMT_fscore > redchaths;
    iDat.RedChaMean = mean(RedCha(:, :, iDat.PMT_fscore == 1), 3);
else
    iDat.PMT_fscore = sum(iDat.PMT_fscore > redchaths, 1);
    iDat.RedChaMean = mean(RedCha(:, :, :, iDat.PMT_fscore == iDat.FrameN), 4);
end

figure('Position', [560 528 801 184])
plot(iDat.PMT_fscore, 'linewidth', 5)