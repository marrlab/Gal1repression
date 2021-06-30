%% main script - run preprocessing, parameter inference and statistical analysis of Gal1 repression project/manuscript

%add path
addpath(genpath(pwd))

%% Preprocessing

%extract single cell information from segmentation files for WT and elp6
getPreprocessing;

%compute non-dividing cells for repression periods for WT and elp6
getNonDividing

%% Parameter estimation - parameter inference of non-repressor and repressor models for WT and elp6 and repressions 1 and 2

%getParameters(model,rep,strain,plt,server,pre_sol)

%model non-repressor, repression 1, strain WT, plots, no server, no
%preloading of parameters
getParameters(1,1,1,1,0,0)

%model repressor, repression 1, strain WT, plots, no server, no
%preloading of parameters
getParameters(2,1,1,1,0,0)

%model non-repressor, repression 2, strain WT, plots, no server, no
%preloading of parameters
getParameters(1,2,1,1,0,0)

%model repressor, repression 2, strain WT, plots, no server, no
%preloading of parameters
getParameters(2,2,1,1,0,0)

%model non-repressor, repression 1, strain elp6, plots, no server, no
%preloading of parameters
getParameters(1,1,2,1,0,0)

%model repressor, repression 1, strain elp6, plots, no server, no
%preloading of parameters
getParameters(2,1,2,1,0,0)

%model non-repressor, repression 2, strain elp6, plots, no server, no
%preloading of parameters
getParameters(1,2,2,1,0,0)

%model repressor, repression 2, strain elp6, plots, no server, no
%preloading of parameters
getParameters(2,2,2,1,0,0)

%% Statistical Analysis - comparison of GFP0 between repressor and non-repressor model

%compGFP0(strain1,rep1)

%WT, repression 1
compGFP0(1,1)

%WT, repression 2
compGFP0(1,2)

%elp6, repression 1
compGFP0(2,1)

%elp6, repression 2
compGFP0(2,2)

%% Statistical Analysis - comparison of estimated parameters

%compParameters(rep1,strain1,rep2,strain2,paired)

%repression 1 vs 2 for WT, paired
compParameters(1,1,2,1,1);

%repression 1 vs 2 for elp6, paired
compParameters(1,2,2,2,1);

%WT vs elp6 for repression 1, not paired
compParameters(1,1,1,2,0);

%WT vs elp6 for repression 2, not paired
compParameters(2,1,2,2,0);

%% Statistical Analysis - comparison of GFP0 and tdelay

%comp_GFP0_tdelay(rep1,strain1,rep2,strain2)

%repression 1 vs 2 for WT
comp_GFP0_tdelay(1,1,2,1)

%repression 1 vs 2 for elp6
comp_GFP0_tdelay(1,2,2,2)

%WT vs elp6 for repression 1
comp_GFP0_tdelay(1,1,1,2)

%WT vs elp6 for repression 2
comp_GFP0_tdelay(2,1,2,2)
