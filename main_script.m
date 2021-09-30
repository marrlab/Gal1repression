%% main script - run preprocessing, parameter inference and statistical analysis of Gal1 repression project/manuscript

%add path
addpath(genpath(pwd))

%% Preprocessing

%extract single cell information from segmentation files for WT and elp6
%!!here: manually add Phylocell path  in script getPreprocessing.m line 6!! 
getPreprocessing;

%compute non-dividing cells for repression periods for WT and elp6
getNonDividing

%% Parameter estimation - parameter inference of non-repressor and repressor models for WT and elp6 and repressions 1 and 2

%getParameters(model,rep,strain,plt,server,pre_sol)

%model = 1 or 2 / model == 1 - non-repressor model / model == 2 - repressor model
%rep = 1 or 2 / rep == 1 - repression 1 / rep == 2 - repression 2
%strain = 1 or 2 / strain == 1 - WT / strain == 2 - elp6
%plt = 0 or 1 / plt == 0 - no plots / plt == 1 - plots (data and fits)
%server = 0 or 1 / server == 0 - not on server / server == 1 - on server
%pre_sol = 0 or 1 / pre_sol == 0 - do not load precomputed parameter estimates / pre_sol == 1 - load precomputed parameter estimates

%!!add path to PESTO manually to script getParameters.m line 15 and potentially server paths lines 19 and 21!!

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

%strain1 = 1 - WT / = 2 - elp6
%rep1 = 1 - repression 1 / = 2 - repression 2

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

%rep1 = 1 - repression 1 / = 2 - repression 2
%strain1 = 1 - WT / = 2 - elp6
%rep2 = 1 - repression 1 / = 2 - repression 2
%strain2 = 1 - WT / = 2 - elp6
%paired = 0 - false / = 1 - true

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

%rep1 = 1 - repression 1 / = 2 - repression 2
%strain1 = 1 - WT / = 2 - elp6
%rep2 = 1 - repression 1 / = 2 - repression 2
%strain2 = 1 - WT / = 2 - elp6

%repression 1 vs 2 for WT
comp_GFP0_tdelay(1,1,2,1)

%repression 1 vs 2 for elp6
comp_GFP0_tdelay(1,2,2,2)

%WT vs elp6 for repression 1
comp_GFP0_tdelay(1,1,1,2)

%WT vs elp6 for repression 2
comp_GFP0_tdelay(2,1,2,2)
