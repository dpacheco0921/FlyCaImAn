%% demo
%% Move to folder and Download demo data
tDir = 'C:\Users\Diego\Dropbox\LabScripts\myRepos\CaImProPi\demodata';
cd(tDir)

url = 'https://www.dropbox.com/s/1s2h6yigfmdhodf/20161129.zip?dl=1';
filename = '20161129.zip';
outfilename = websave(filename, url);
unzip(outfilename);
clear url outfilename
    
%% Test batch_tiff2mat
FolderName = {'20161129'}; FileName = [];
m2vpar = [];
m2vpar.SpMode = '3DxT_song';
m2vpar.ch2save = [1 2];
m2vpar.Zres = 2; % Z resolution in um
m2vpar.pixelsym = 0;
batch_tiff2mat(FolderName, FileName, m2vpar)

FolderName = {'20190712'}; FileName = [];
m2vpar.SpMode = '3DxT';
batch_tiff2mat(FolderName, FileName, m2vpar)

%% Test batch_collectmetada
FolderName = {'20161129'}; FileName = [];
cmpar = [];
cmpar.pgate = 1;
cmpar.pgates = 1;
cmpar.mode = 1;
cmpar.minframet = 900;
batch_collectmetada(FolderName, FileName, cmpar)

FolderName = {'20190712'}; FileName = [];
cmpar.minframet = 900;
batch_collectmetada(FolderName, FileName, cmpar)

%% Test batch_NoRMCorre
FolderName = []; FileName = [];
mcpar = [];
mcpar.debug = 0;
mcpar.rigidg = 1;
mcpar.nrigidg = 0;
mcpar.parallel = 1;
mcpar.stack2del = 1:3;
mcpar.sgate = 1; %(1 = smooth and zero, 2 = smooth, 3 = zeroing, 0 = raw)
batch_NoRMCorre(FolderName, FileName, mcpar)

FolderName = {'20190712'}; FileName = [];
batch_NoRMCorre(FolderName, FileName, mcpar)

batch_plotShiftPerStack([], [], mcpar)

%% Test batch_SpaTemp_ResFilt_3D
FolderName = []; FileName = [];
stpar = [];
stpar.sigma = [];
stpar.size = [];
stpar.newtimeres = 0.5;
stpar.time = [];
stpar.debug = 1;

batch_SpaTemp_ResFilt_3D(FolderName, FileName, stpar)
