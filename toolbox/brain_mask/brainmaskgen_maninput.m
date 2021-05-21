%% Generating brain mask
%% Manual curation/editing of mask
% 1.1) test if generated 'wDat.bF_ths' is good
display(wDat.bF_ths) % display current threshold
display([min(wDat.GreenChaMean(:)) ...
    max(wDat.GreenChaMean(:))])
im_range = [min(wDat.GreenChaMean(:)) ...
    max(wDat.GreenChaMean(:))*0.7];
% generate mask using this threshold
brainmask = wDat_generatemask(wDat, 5^3, bGreen, 0); 

% 1.2) slide through planes and contour of mask per plane
pi = [];
pi.iter = 0;
pi.range = [im_range(1) im_range(1) + 100];
pi.maskout = 0;
pi.lag = 0.1;
pi.Y1 = brainmask;
pi.figpos = genfigpos(1, 'nw', [879 659]);

% 1.3 optional save this movie to local directory (uncomment the code below)
% pi.vgate = 1;
% pi.vname = [pwd, filesep, f2run, '_brainmask_temp'];
% pi.frate = 1;

slice3Dmatrix(wDat.GreenChaMean, pi)
% slice3Dmatrix(wDat.RedChaMean, pi)
% 1.3) if it looks good then just continue (F5 or Run)
% 1.4) if not happy, edit wDat.bF_ths, then repeat 1.1), 1.2) and 1.3)

%% save stitched image
pi.vgate = 1;
pi.vname = ['C:\Users\diego\Dropbox\tempfiles\stitch', filesep, f2run];
pi.frate = 1;
slice3Dmatrix(wDat.GreenChaMean, pi)

%% Use Red channel for mask
wDat.bF_ths = 21;
display(wDat.bF_ths) % display current threshold
bRed = imblur(wDat.RedChaMean, [2 2 2], [5 5 3], 3);

% generate mask using this threshold
brainmask = wDat_generatemask(wDat, 5^3, bRed);

% 1.2) slide through planes and contour of mask per plane
pi = [];
pi.iter = 0;
pi.range = [-10 200]; 
pi.maskout = 0;
pi.lag = 0.001;
pi.Y1 = brainmask;
pi.figpos = genfigpos(1, 'nw', [879 659]);
slice3Dmatrix(wDat.GreenChaMean, pi)

%% Plot mean image and generated/saved mask
% 1.2) slide through planes and contour of mask per plane
im_range = [min(wDat.GreenChaMean(:)) ...
    max(wDat.GreenChaMean(:))*0.7];

pi = [];
pi.iter = 0;
pi.range = [im_range(1) im_range(1) + 100];
pi.maskout = 0;
pi.lag = 0.1;
pi.Y1 = wDat.bMask;
pi.figpos = genfigpos(1, 'nw', [879 659]);

slice3Dmatrix(wDat.GreenChaMean, pi)
