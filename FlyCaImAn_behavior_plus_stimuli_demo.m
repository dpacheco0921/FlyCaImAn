%% demo for analysis of ball video + stimuli
%% 1) add paths
% it assumes you have already add the repository folders to your path
addpath(genpath(pwd))

% add paths of all dependencies
% CaImAn, NoRMCorre, CMTK_matlab_wrapper

%% 2) Move to folder and Download demo data
tDir = strrep(which('FlyCaImAn_demo'), 'FlyCaImAn_demo.m', '');
cd(tDir)

url = 'https://www.dropbox.com/s/5jql2d2tptjdp6e/20200910.zip?dl=1';
filename = '20200910.zip';

if ~exist('demodata', 'dir')
    mkdir('demodata')
end
cd demodata

outfilename = websave(filename, url);
unzip(outfilename);
clear url outfilename

%% 3) process behavior videos with fictract
% 1.2.1) Generate fictrac input files
FolderName = '20200910'; FileName = [];
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

%% 4) add stimuli metadata

FolderName = '20200910'; FileName = [];
m2vpar = [];
m2vpar.SpMode = '3DxT_song_prv';
m2vpar.ch2save = [1 2];
m2vpar.Zres = [];
m2vpar.ch2plot = [];

batch_tiff2mat(FolderName, FileName, m2vpar)
% generates '20200910_1_1_metadata.mat'
% generates '20200910_1_1_rawdata.mat' with empty frames

cmpar = [];
cmpar.pgate = 1;
cmpar.pgates = 1;
cmpar.mode = 1;

batch_collectmetada(FolderName, FileName, cmpar)
% passes stimuli info to '20200910_1_1_metadata.mat'

stpar = [];
stpar.sigma = [];
stpar.size = [];
stpar.newtimeres = 0.1;
stpar.time = [];
stpar.debug = 1;
stpar.direction = 'invert';
stpar.fshift = [6 6];

batch_SpaTemp_ResFilt(FolderName, FileName, stpar)
% passes both fictrac and stimuli info to '20200910_1_1_metadata.mat'
%   stored in wDat

%% 5) plot stim vs fictrac

FolderName = '20200910'; FileName = [];
fictracpar = [];
fictracpar.dir_depth = 1;
fictracpar.timeres = 0.05;

batch_plot_fictrac_results(FolderName, FileName, fictracpar)
