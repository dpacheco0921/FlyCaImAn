% re-structure functions
% add ref and float indeces
% update all code that depends on this
% seems like the jitter was scaled to pixel units? why? do it again.
% delete variables not used
% make extra functions not dependent on obj.
% who uses: (replace)
% get_point_normal --> to_get_point_normal
% plotintersect_3D --> to_plot_point_3D
% get_avg_norm --> to_get_avg_norm
% get_xyz_jitter --> to_get_xyz_jitter

% make plot inside to_get_avg_norm compatible
% test if it runs
% update plotorthoplanes_2D plotorthoplanes_3D

%% _
cd('V:\Diego\LabScripts')
mgrep('traceObj')

edit  V:\Diego\LabScripts\paper.code\AuditoryPaper\figure_1_panel_c.m

%%

edit V:\Diego\LabScripts\paper.code\AuditoryPaper\figure_1_panel_c.m
obj_p1 = traceobj; % p1s and pc1s (load and smooth)

obj_pc2l = traceobj; % pc2l (load and smooth)

obj_pc2m = traceobj; % pc2m (load and smooth)

 
V:\Diego\LabScripts\paper.code\AuditoryPaper\figure_1_panel_c_2.m
% 2.1) reformat traces using traceobj

tobj = traceobj;

obj_p1_ibnwb = traceobj;

obj_p1 = traceobj;

 
V:\Diego\LabScripts\paper.code\AuditoryPaper\figure_1_supp_1.m
obj_p1 = traceobj; % p1s and pc1s (load and smooth)

obj_pc2l = traceobj; % pc2l (load and smooth)

obj_pc2m = traceobj; % pc2m (load and smooth)

obj_ta1 = traceobj; % ta1 (load and smooth)

obj_ta2 = traceobj; % ta2 (load and smooth)

obj_pmp7 = traceobj; % pmp7 (load and smooth)

 
V:\Diego\LabScripts\paper.code\AuditoryPaper\figure_1_supp_3.m
obj_p1_all = traceobj;

 
V:\Diego\LabScripts\preprocessCaIm\analysis_suite\st_plot_stim_pref_triad.m
obj_trace = traceobj;

 
V:\Diego\LabScripts\preprocessCaIm\tree_postpros\scriptsinfo_traceObj.m
%% 3) using cmtk & traceobj to reformat traces in matlab

obj = traceobj;

%% 4) using traceobj to vizualize and perform basic operations with traces in matlab

obj = traceobj; % make object

 
V:\Diego\LabScripts\preprocessCaIm\tree_postpros\to_collect_atlas_info.m
save([strrep(which('traceobj'), 'traceobj.m', ''), 'atlas_meta.mat'], 'atlas')

 
V:\Diego\LabScripts\preprocessCaIm\tree_postpros\to_dist_to_refpoimt.m
% notes: see traceobj

 
V:\Diego\LabScripts\preprocessCaIm\tree_postpros\to_get_plane2segline_intersect.m
% notes: see traceobj

 
V:\Diego\LabScripts\preprocessCaIm\tree_postpros\to_get_planeintersects.m
% notes: see traceobj

 
V:\Diego\LabScripts\preprocessCaIm\tree_postpros\to_pick_xyz_object.m
%   obj: traceobj obj

 
V:\Diego\LabScripts\preprocessCaIm\tree_postpros\to_plot_point_3D.m
%   obj: traceobj obj

 
V:\Diego\LabScripts\preprocessCaIm\tree_postpros\todo.m
mgrep('traceobj')

 
V:\Diego\LabScripts\preprocessCaIm\tree_postpros\traceObj.m
classdef traceobj < handle

        function obj = traceobj

            if exist([strrep(which('traceobj'), 'traceobj.m', ''), 'atlas_meta.mat'], 'file')

                load([strrep(which('traceobj'), 'traceobj.m', ''), 'atlas_meta.mat'])

 
V:\Diego\LabScripts\scriptsperData\atlasrelated\brain_nsybIVAi2IBNWB_accuracy.m
% reformat traces using traceobj

tobj = traceobj;

% 2.1) reformat traces using traceobj

tobj = traceobj;

obj_p1_nsybivai = traceobj;

 
V:\Diego\LabScripts\scriptsperData\atlasrelated\neuron_NP2631FruFLP.m
% reformat traces using traceobj

obj = traceobj;

 
V:\Diego\LabScripts\scriptsperData\atlasrelated\neuron_all.m
%% 3) reformat traces using cmtk and traceobj

 
V:\Diego\LabScripts\scriptsperData\atlasrelated\neuron_dsx.m
% reformat traces using traceobj

obj = traceobj;

 
V:\Diego\LabScripts\scriptsperData\atlasrelated\neuron_dsx_neuronjitter.m
obj_p1 = traceobj;

    eval([objname_{i}, ' = traceobj;'])

 
V:\Diego\LabScripts\scriptsperData\atlasrelated\neuron_pMP7_neuronjitter.m
obj_pmp7 = traceobj;