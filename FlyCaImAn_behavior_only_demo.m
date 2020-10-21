%% demo for analysis of ball video only
%% 1) add paths
% it assumes you have already add the repository folders to your path
addpath(genpath(pwd))

% add paths of all dependencies
% CaImAn, NoRMCorre, CMTK_matlab_wrapper

%% 2) Move to folder and Download demo data
tDir = strrep(which('FlyCaImAn_demo'), 'FlyCaImAn_demo.m', '');
cd(tDir)

url = 'https://www.dropbox.com/s/n65nq41dn0lcsoe/20200911.zip?dl=1';
filename = '20200911.zip';

if ~exist('demodata', 'dir')
    mkdir('demodata')
end
cd demodata

outfilename = websave(filename, url);
unzip(outfilename);
clear url outfilename
    
%% 3) process behavior videos with fictract
% 1.2.1) Generate fictrac input files
FolderName = '20200911'; FileName = [];
setupID = 'Alevel_2p1';
% readout a tiff for mask or extract from video mp4
im_mode = 1;

batch_gen_fictrac_input_files(FolderName, FileName, setupID, im_mode)
% this function generates:
% 20200911_1_1_maskim.tiff
% 20200911_1_1_calibration-transform.dat
% which are used by fictrac

% 1.2.2) Run fictrac on each file per current directory
fictracpars = [];
fictracpars.redo = 1;
fictracpars.do_search_flag = 1; % to do global search
%fictracpars.load_template_flag = 1; % to load pre-generated template
%fictracpars.template_fn = 'template_im.jpg'; % pre-generated template image

batch_fictrac_perfile(FolderName, FileName, fictracpars)
% this function generates:
% 20200911_1_1_fictrac.txt
% 20200911_1_1_fictrac-debug.avi
% 20200911_1_1_fictrac-raw_frames.avi
% *.avi are for visualization and debugging
% 20200911_1_1_vDat.mat
% fictrac variables are saved as 'fictrac'
