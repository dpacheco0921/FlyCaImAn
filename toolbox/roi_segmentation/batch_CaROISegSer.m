function batch_CaROISegSer(fname, inputparams, ...
    serverid, jobpart, memreq, patchtype, ...
    corenum, roi_n_init, stitch_flag)
% batch_CaROISegSer: batch function to ROI segment files within
%   subdirectory
%
% Usage:
%   batch_CaROISegSer(fname, inputparams, ...
%       serverid, jobpart, memreq, patchtype, ...
%       corenum, roi_n_init, stitch_flag)
%
% Args:
%   fname: name of files to run
%   inputparams: string defining script to generate 
%       ROI segmentation parameters
%       (default, roiseg_3D_dense_fr_2Hz_z2)
%   serverid: server ID 'int', 'spock', 'della'
%       (deafult, 'spock')
%   jobpart: type of job to run
%       1) run roi segmentation for patches of whole volume
%       2) compile patches
%       (deafult, 1)
%       3) extract processed signal
%   memreq: RAM memory to request
%       (deafult, 48)
%   patchtype: flag on how to define patches
%       (default, 1, it splits it into nXmXt segments)
%       (alternatively, for stitched volumes, it splits them into segments
%       of the same stack)
%   corenum: number of cores to request
%       (deafult, 4)
%   roi_n_init: number of ROIs to segment
%       (deafult, [])
%   stitch_flag: flag indicating that file is stitched (used to then find input files)
%       (deafult, 1)
% 
% Notes:
% it requires that the script is run in the data folder
% internal variable with important parameters CaRp
%   CaRp (controls which files to load and details of server usage)
%         %%%%%% input files %%%%%%
%        (fisuffix: input data suffix)
%        (fimetsuffix: metadata suffix)
%        (rfisuffix: refdata suffix)
%         %%%%%% function related %%%%%%
%        (functype: functions available to run)
%        (maxjobs: maximun number of segments per submition)
%         %%%%%% directory related %%%%%%
%        (tDir: directory where the output files are saved)
%        (fDir: function directory)
%        (cDir: default data folder)

% running ROI segmentation
CaRp = [];
CaRp.fisuffix = [];
CaRp.fimetsuffix = [];
CaRp.rfisuffix = [];
CaRp.functype = {'CaROISegSer.m'};
CaRp.cDir = pwd;
maxjobs = 4;

if exist('fname', 'var') && ~isempty(fname)
   fname = [];
end

if exist('serverid', 'var') && ~isempty(serverid)
    serverid = 'spock';
end

if exist('jobpart', 'var') && ~isempty(jobpart)
    jobpart = 1;
end

if exist('memreq', 'var') && ~isempty(memreq)
    memreq = 48;
end

if exist('corenum', 'var') && ~isempty(corenum)
    corenum = 4;
end

if exist('roi_n_init', 'var') && ~isempty(roi_n_init)
    roi_n_init = [];
end

if ~exist('inputparams', 'var') || isempty(inputparams)
    inputparams = 'roiseg_3D_dense_fr_2Hz_z2';
end

if ~exist('patchtype', 'var') || isempty(patchtype)
    patchtype = 0;
end

if ~exist('stitch_flag', 'var') || isempty(stitch_flag)
    stitch_flag = 1;    
end

if jobpart == 2
    corenum = 4;
end

% 1) define string patter to find input files
if stitch_flag
    % input files for stitched data
    CaRp.fimetsuffix = '_prosmetadata';
    CaRp.fisuffix = '_prosdata';
    CaRp.rfisuffix = '_prosref';
else
    % input files for non-stitched data
    CaRp.fimetsuffix = '_metadata';
    CaRp.fisuffix = '_rawdata';
    CaRp.rfisuffix = '_refdata';
end

% 2) set deault roi seg parameters
eval(inputparams)
% existing settings:
% roiseg_2D_fr_10Hz (2D data)
% roiseg_3D_dense_fr_2Hz_z2
%   used for whole brain (song/opto) / dense labelling (9-15 min recordings)
% roiseg_3D_dense_fr_1Hz_z1
%   used for whole brain (opto) / dense labelling (9-15 min recordings)
% roiseg_3D_dense_fr_01Hz_z2
%   used for whole brain (opto) / dense labelling (>70 min recordings)

% 3) set a custom number of ROIs
%   this overwrites roiparams.K defined earlier at 2) 
if ~isempty(roi_n_init)
    roiparams.K = roi_n_init;
end

if ispc || ismac
    % running locally
    corenum = 4;
    maxjobs = 4;
end

% get scratch (temporary) and bucket (permanent) directories
[~, username, ~, scratchdir, ~] = ...
    user_defined_directories(serverid);
if ~exist([scratchdir, sep, 'jobsub', sep, 'roirel'], 'dir')
    mkdir([scratchdir, sep, 'jobsub', sep, 'roirel']);
end
tDir = [fastscratchdir, filesep, 'jobsub', filesep, 'roirel'];

% Find input files withing data folder CaRp.cDir
[filename, patchidx] = getinputfiles(...
    CaRp.fisuffix, CaRp.fimetsuffix, ...
   fname, maxjobs, jobpart, ....
    roiparams, patchtype);

% Submitting job or running it
if ~isempty(filename)
    
    fprintf(['Files to preprocess: ', num2str(numel(filename)), '\n'])
    fprintf('roiparams used: \n')
    roiparams
    
    fprintf('Files to run \n')
    for f_Num = 1:numel(filename)
        fprintf([filename{f_Num}, '\n'])
    end
    
    % saving params (in outputdir)
    param_file = datestr(now, 'yymmddHHMMSS');
    
    % pass some input variables
    CaRp.serId = serverid;
    CaRp.corenum = corenum;
    CaRp.jobpart = jobpart;
    p = CaRp;
    roiparams.patchtype = patchtype;

    % save in the current folder
    save([tDir, filesep, param_file, '_impre.mat'], ...
        'filename', 'patchidx', 'p', 'roiparams')
    save([tDir, filesep, param_file, '_status.mat'], ...
        'filename')
    
    % Executing File
    numT = numel(filename); 
    cd(tDir)
    
    CaRp
    
    submitjob(param_file, tDir, ...
        username, corenum, ...
        jobpart, serverid, ...
        numT, memreq)
    
else
    
    fprintf('No files selected\n')
    
end

end

function [filename, patchidx] = ...
    getinputfiles(fisuffix, fimetsuffix, ...
    fname, maxjobs, jobpart, ....
    roiparams, patchtype)
% getinputfiles: get all files to run
%
% Usage:
%   [filename, patchidx] = ...
%       getinputfiles(fisuffix, fimetsuffix, ...
%       fname, maxjobs, jobpart, ....
%       roiparams, patchtype)
%
% Args:
%   fisuffix: input data suffix
%   fimetsuffix: input metadata suffix
%   fname: pattern to match
%   maxjobs: maximun number of segments per submition
%   jobpart: type of job to run
%       1) run roi segmentation for patches of whole volume (generate RESULTS)
%       2) compile patches
%       (deafult, 1)
%       3) extract processed signal
%   roiparams: roi parameters
%   patchtype: how to define patches
%
% Notes

% look for input files in the local dir
patchidx = [];
filename = [];

f2run = rdir(['.', filesep, '*', fisuffix, '.mat']);
f2run = str2match(fname, f2run);
f2run = {f2run.name};
f2run = strrep(f2run, ['.', filesep], '');
f2run = strrep(f2run, '.mat', '');
f2run = strrep(f2run, fisuffix, '');

if jobpart == 1 || jobpart == 3
    
    for ii = 1:numel(f2run)
        
        load([f2run{ii}, fimetsuffix, '.mat'], 'wDat')
        
        % get number of patches
        if patchtype == 0
            
            % split into a pre-defined size patches
            patchnum = length(construct_patches(wDat.vSize, ...
                roiparams.patch_size, roiparams.overlap));
            
        else
            
            % split into patches recorded simultaneously
            %   (uses wDat.Zstitch.Zidx)
            patchnum = numel(construct_patches_seg(wDat.vSize, ...
                    wDat.Zstitch.Zidx));
                
        end
        
        % patches = construct_patches(sizY,siz,overlap,min_size)
        tempidx = arrayfun(@(i) i:min((i+maxjobs-1), patchnum), ...
            1:maxjobs:patchnum, 'UniformOutput', false)';
        tempname = arrayfun(@(i) f2run{ii}, 1:length(tempidx), ...
            'UniformOutput', false)';
        clear wDat
        
        patchidx = cat(1, patchidx, tempidx);
        filename = cat(1, filename, tempname);
        
    end
    
elseif jobpart == 2
    
    filename = f2run;
    
end

end

function submitjob(name, tDir, ...
    username, corenum, jobpart, ...
    serverid, numT, memreq)
% submitjob: submit jobs to spock/della
%
% Usage:
%   submitjobsubmitjob(name, tDir, ...
%       username, corenum, jobpart, ...
%       serverid, numT)
%
% Args:
%   name: name of matfile with parameters to use
%   tDir: directory where the output files are saved
%   username: used to update directories to use
%   corenum: maximun number of cores to use per task
%   jobpart: type of job to run
%       1) run roi segmentation for patches of whole volume (generate RESULTS)
%       2) compile patches
%       (deafult, 1)
%       3) extract processed signal
%   serverid: server ID 'int', 'spock', 'della'
%       (deafult, 'spock')
%   numT: number of jobs
%   memreq: RAM memory to request
%       (deafult, 48)

functype = 'CaROISegSer.m';

switch serverid
    case 'spock'
        
        LogFileName = fullfile([name, '.slurm']);
        if exist(LogFileName, 'file')
            delete(LogFileName)
        end
        
        % open/create log file
        fid = fopen(LogFileName, 'a+');
        fprintf(fid, '#!/bin/bash\n\n');
        fprintf(fid, '#SBATCH -N 1\n');
        fprintf(fid, ['#SBATCH --cpus-per-task=', num2str(corenum), '\n']);
        if jobpart == 1 || jobpart == 3
            fprintf(fid, '#SBATCH --time=4:00:00\n');
            fprintf(fid, ['#SBATCH --mem=', num2str(memreq), '000\n']);            
        elseif jobpart == 2
            fprintf(fid, '#SBATCH --time=4:00:00\n');
            fprintf(fid, ['#SBATCH --mem=', num2str(memreq), '000\n']);            
        end
        fprintf(fid, '#SBATCH --mail-type=END\n');
        fprintf(fid, '#SBATCH --mail-user=dpacheco@princeton.edu\n');
        fprintf(fid, ['#SBATCH --array=1-', num2str(numT), '\n\n']);
        
        fprintf(fid, 'module load matlab/R2016b\n');
        fprintf(fid, '# Create a local work directory\n');
        fprintf(fid, 'mkdir -p /tmp/$USER-$SLURM_JOB_ID\n');
        fprintf(fid, ['matlab -nodesktop -nodisplay -nosplash -r "', ...
            functype(1:end-2),'(''', name, ''',''', serverid, ''')"\n']);
        fprintf(fid, '# Cleanup local work directory\n');
        fprintf(fid, 'rm -rf /tmp/$USER-$SLURM_JOB_ID\n');
        
        % close log file
        fclose(fid);
        
    case 'della'
        
        LogFileName = fullfile([name, '.slurm']);
        if exist(LogFileName, 'file')
            delete(LogFileName)
        end
        
        % open/create log file
        fid = fopen(LogFileName, 'a+');
        fprintf(fid, '#!/bin/bash\n\n');
        fprintf(fid, '#SBATCH -N 1\n');
        fprintf(fid, ['#SBATCH --cpus-per-task=', num2str(corenum), '\n']);
        if jobpart == 1 || jobpart == 3
            fprintf(fid, '#SBATCH --time=4:00:00\n');
            fprintf(fid, ['#SBATCH --mem=', num2str(memreq), '000\n']);
        elseif jobpart == 2
            fprintf(fid, '#SBATCH --time=4:00:00\n');
            fprintf(fid, ['#SBATCH --mem=', num2str(memreq), '000\n']);
        end
        fprintf(fid, '#SBATCH --mail-type=END\n');
        fprintf(fid, '#SBATCH --mail-user=dpacheco@princeton.edu\n');
        fprintf(fid, ['#SBATCH --array=1-', num2str(numT), '\n\n']);
        
        fprintf(fid, 'module load matlab/R2016b\n');
        fprintf(fid, '# Create a local work directory\n');
        fprintf(fid, 'mkdir -p /tmp/$USER-$SLURM_JOB_ID\n');
        fprintf(fid, ['matlab -nodesktop -nodisplay -nosplash -r "', ...
            functype(1:end-2),'(''', name, ''',''', serverid, ''')"\n']);
        fprintf(fid, '# Cleanup local work directory\n');
        fprintf(fid, 'rm -rf /tmp/$USER-$SLURM_JOB_ID\n');
        
        % close log file
        fclose(fid);
        
    otherwise % internal run
        
        cd(tDir)
        for f_run = 1:numT
           CaROISegSer(name, serverid, f_run); 
           fprintf('\n\n **************************************************** \n\n')
           cd(tDir)
        end
        
end

end
