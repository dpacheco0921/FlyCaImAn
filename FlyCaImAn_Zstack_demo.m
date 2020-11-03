%% demo for analysis of ball video only
%% 1) add paths
% it assumes you have already add the repository folders to your path
addpath(genpath(pwd))

% add paths of all dependencies
% CaImAn, NoRMCorre, CMTK_matlab_wrapper

%% 2) Move to folder and Download demo data
tDir = strrep(which('FlyCaImAn_demo'), 'FlyCaImAn_demo.m', '');
cd(tDir)

url = 'https://www.dropbox.com/s/purvtala87l1hpe/20200322.zip?dl=1';
filename = '20200322.zip';

if ~exist('demodata', 'dir')
    mkdir('demodata')
end
cd demodata

outfilename = websave(filename, url);
unzip(outfilename);
clear url outfilename
    
%% 3) process z-stack tiffs
%% 3.1) Collect structural Zstacks

FolderName = '20200322'; FileName = [];
zstackpar = [];
zstackpar.ths = 50;
zstackpar.redo = 0;
zstackpar.integrate_flag = 1;
zstackpar.shift_f = [50 50];

batch_zstacktiff2mat(FolderName, FileName, zstackpar)
% batch_zstacktiff2mat(FolderName, FileName, params)

%% 3.2) Stitch structural Zstack

FolderName = '20200322'; FileName = [];
stitchpar.refcha = 1;
stitchpar.peaknum = 5;
stitchpar.init_xyz = [0 0 0 0 0 0];
stitchpar.redo = 1;
stitchpar.debug_flag = 0;
stitchpar.ijscript = 'refstitcher_1to2_3to4.ijm';

batch_stitcher_nrrd(FolderName, FileName, stitchpar)
%batch_stitcher_nrrd(FolderName, FileName, stitchpar)
