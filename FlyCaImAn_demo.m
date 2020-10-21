%% Demo for preprocessing imaging and/or behavioral data
%% 1) add paths
% it assumes you have already add the repository folders to your path
addpath(genpath(pwd))

% add paths of all dependencies
% CaImAn, NoRMCorre, CMTK_matlab_wrapper

%% 1) preprocess behavior only
edit FlyCaImAn_behavior_only_demo

%% 2) preprocess behavior + stimuli (song)
edit FlyCaImAn_behavior_plus_stimuli_demo.m

%% 3) preprocess imaging + behavior + stimuli (song)
edit FlyCaImAn_behavior_plus_stimuli_plus_imaging_demo.m

%% 4) preprocess imaging only (song) (with stitching across sequentially imaged segments)
edit FlyCaImAn_imaging_only_demo.m
