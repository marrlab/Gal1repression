%% main script - run preprocessing, parameter inference and statistical analysis of Gal1 repression project/manuscript

%add path
addpath(genpath(pwd))

%% Preprocessing

%reformat the output files of Cell ACDC (.csv) to match the output from
%PhyloCell 
%only required to run if one is working on the raw data outputs from Cell ACDC 
%however, we here provide the already reformatted outputs 
% getReformatting;

%extract single cell information from segmentation files for WT and elp6
getPreprocessing;

%compute non-dividing cells for repression periods for WT and elp6
% getNonDividing(iexp), where iexp = 1 - main experiment / iexp = 2 - replicate experiment
%main experiment
getNonDividing(1)
%replicate experiment
getNonDividing(2)

%% Parameter estimation - parameter inference of non-repressor and repressor models for WT and elp6 and repressions 1 and 2

%getParameters(model,rep,strain,plt,server,pre_sol,iexp)

%model = 1 or 2 / model == 1 - non-repressor model / model == 2 - repressor model
%rep = 1 or 2 / rep == 1 - repression 1 / rep == 2 - repression 2
%strain = 1 or 2 / strain == 1 - WT / strain == 2 - elp6
%plt = 0 or 1 / plt == 0 - no plots / plt == 1 - plots (data and fits)
%server = 0 or 1 / server == 0 - not on server / server == 1 - on server
%pre_sol = 0 or 1 / pre_sol == 0 - do not load precomputed parameter estimates / pre_sol == 1 - load precomputed parameter estimates
%iexp = 1 - main experiment / iexp = 2 - replicate experiment

%main experiment
%model non-repressor, repression 1, strain WT, plots, no server, no
%preloading of parameters, main experiment
getParameters(1,1,1,1,0,0,1)

%model repressor, repression 1, strain WT, plots, no server, no
%preloading of parameters, main experiment
getParameters(2,1,1,1,0,0,1)

%model non-repressor, repression 2, strain WT, plots, no server, no
%preloading of parameters, main experiment
getParameters(1,2,1,1,0,0,1)

%model repressor, repression 2, strain WT, plots, no server, no
%preloading of parameters, main experiment
getParameters(2,2,1,1,0,0,1)

%model non-repressor, repression 1, strain elp6, plots, no server, no
%preloading of parameters, main experiment
getParameters(1,1,2,1,0,0,1)

%model repressor, repression 1, strain elp6, plots, no server, no
%preloading of parameters, main experiment
getParameters(2,1,2,1,0,0,1)

%model non-repressor, repression 2, strain elp6, plots, no server, no
%preloading of parameters, main experiment
getParameters(1,2,2,1,0,0,1)

%model repressor, repression 2, strain elp6, plots, no server, no
%preloading of parameters, main experiment
getParameters(2,2,2,1,0,0,1)

%replicate experiment
%model non-repressor, repression 1, strain WT, plots, no server, no
%preloading of parameters, replicate experiment
getParameters(1,1,1,1,0,0,2)

%model repressor, repression 1, strain WT, plots, no server, no
%preloading of parameters, replicate experiment
getParameters(2,1,1,1,0,0,2)

%model non-repressor, repression 2, strain WT, plots, no server, no
%preloading of parameters, replicate experiment
getParameters(1,2,1,1,0,0,2)

%model repressor, repression 2, strain WT, plots, no server, no
%preloading of parameters, replicate experiment
getParameters(2,2,1,1,0,0,2)

%model non-repressor, repression 1, strain elp6, plots, no server, no
%preloading of parameters, replicate experiment
getParameters(1,1,2,1,0,0,2)

%model repressor, repression 1, strain elp6, plots, no server, no
%preloading of parameters, replicate experiment
getParameters(2,1,2,1,0,0,2)

%model non-repressor, repression 2, strain elp6, plots, no server, no
%preloading of parameters, replicate experiment
getParameters(1,2,2,1,0,0,2)

%model repressor, repression 2, strain elp6, plots, no server, no
%preloading of parameters, replicate experiment
getParameters(2,2,2,1,0,0,2)

%% Statistical Analysis - comparison of GFP0 between repressor and non-repressor model

%compGFP0(strain1,rep1,iexp)

%strain1 = 1 - WT / = 2 - elp6
%rep1 = 1 - repression 1 / = 2 - repression 2
%iexp = 1 - main experiment / iexp = 2 - replicate experiment

%main experiment
%WT, repression 1, main experiment
compGFP0(1,1,1)

%WT, repression 2, main experiment
compGFP0(1,2,1)

%elp6, repression 1, main experiment
compGFP0(2,1,1)

%elp6, repression 2, main experiment
compGFP0(2,2,1)

%replicate experiment
%WT, repression 1, replicate experiment
compGFP0(1,1,2)

%WT, repression 2, replicate experiment
compGFP0(1,2,2)

%elp6, repression 1, replicate experiment
compGFP0(2,1,2)

%elp6, repression 2, replicate experiment
compGFP0(2,2,2)

%% Statistical Analysis - comparison of estimated parameters

%compParameters(rep1,strain1,rep2,strain2,paired,iexp)

%rep1 = 1 - repression 1 / = 2 - repression 2
%strain1 = 1 - WT / = 2 - elp6
%rep2 = 1 - repression 1 / = 2 - repression 2
%strain2 = 1 - WT / = 2 - elp6
%paired = 0 - false / = 1 - true
%iexp = 1 - main experiment / iexp = 2 - replicate experiment

%main experiment
%repression 1 vs 2 for WT, paired, main experiment
compParameters(1,1,2,1,1,1);

%repression 1 vs 2 for elp6, paired, main experiment
compParameters(1,2,2,2,1,1);

%WT vs elp6 for repression 1, not paired, main experiment
compParameters(1,1,1,2,0,1);

%WT vs elp6 for repression 2, not paired, main experiment
compParameters(2,1,2,2,0,1);

%replicate experiment
%repression 1 vs 2 for WT, paired, replicate experiment
compParameters(1,1,2,1,1,2);

%repression 1 vs 2 for elp6, paired, replicate experiment
compParameters(1,2,2,2,1,2);

%% Statistical Analysis - comparison of GFP0 and tdelay

%comp_GFP0_tdelay(rep1,strain1,rep2,strain2,iexp)

%rep1 = 1 - repression 1 / = 2 - repression 2
%strain1 = 1 - WT / = 2 - elp6
%rep2 = 1 - repression 1 / = 2 - repression 2
%strain2 = 1 - WT / = 2 - elp6
%iexp = 1 - main experiment / iexp = 2 - replicate experiment

%main experiment
%repression 1 vs 2 for WT, main experiment
comp_GFP0_tdelay(1,1,2,1,1)

%replicate experiment
%repression 1 vs 2 for WT, replicate experiment
comp_GFP0_tdelay(1,1,2,1,2)

%% Check differetn degradation rates (simulation)

simDeg
