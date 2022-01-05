function getParameters(model,rep,strain,plt,server,pre_sol,iexp)

%model = 1 or 2 / model == 1 - non-repressor model / model == 2 - repressor model
%rep = 1 or 2 / rep == 1 - repression 1 / rep == 2 - repression 2
%strain = 1 or 2 / strain == 1 - WT / strain == 2 - elp6
%plt = 0 or 1 / plt == 0 - no plots / plt == 1 - plots (data and fits)
%server = 0 or 1 / server == 0 - not on server / server == 1 - on server
%pre_sol = 0 or 1 / pre_sol == 0 - do not load precomputed parameter estimates / pre_sol == 1 - load precomputed parameter estimates
%iexp = 1 or 2 / iexp == 1 - main experimetn / iexp == 2 - replicate experiment 

clearvars -except model rep strain plt server pre_sol iexp
clc;

%add path for AMICI and PESTO
if server == 0 %adapt paths
    addpath(genpath('/Users/lea.schuh/Documents/PhD/ICB/Xenopus/PESTO-master'))
    addpath(genpath(pwd))
    if iexp == 1
        load('NonDividing1')
    else
        load('NonDividing2')
    end
else %adapt paths
    addpath(genpath('/home/icb/lea.schuh/PESTO-master'))
    addpath(genpath(pwd))
    if iexp == 1
        load('/home/icb/lea.schuh/Gal1/NonDividing1')
    else
        load('/home/icb/lea.schuh/Gal1/NonDividing2')
    end
end

if ~exist('./Results', 'dir')
    mkdir('./Results')
    addpath(genpath('./Results'))
end

for imodel = model
    for irep = rep
        if plt == 1
            figure
        end
        for istrain = strain
            
            clear scR DA data
            
            if irep == 1
                %normalize data by 1e7
                data = NonDividing{istrain}.r1./1e7;
            else
                data = NonDividing{istrain}.r2./1e7;
            end
            
            %define model parameters and upper and lower boundaries for
            %parameter estimation (in log10)
            %P0 - initital total GFP
            %t1 - tdelay - time point at which production rate is
            %turned off
            %rprod - production rate
            %rdeg - degradation rate
            %noise - width of Gaussian modeling the noise
            if imodel == 1
                parameters.name = {'P0' 'rprod' 'rdeg' 'noise'};
                parameters.min = [-10, -10, -10, -10];
                parameters.max = [1, 1, 1, 1];
                ind = 1:4;
            else
                parameters.name = {'P0' 't1' 'rprod' 'rdeg' 'noise'};
                parameters.min = [-10, -2, -10, -10, -10];
                parameters.max = [1, log10((size(data,2)-1)*3/60), 1, 1, 1];
                ind = 1:5;
            end
            
            %specifiy number of estimated parameters
            parameters.number = length(parameters.name);
            
            %for every total GFP trace do
            for icell = 1:size(data,1)
                
                %display which total GFP trace is currently fitted
                disp(sprintf('%d of %d for strain %d for repression %d and model %d', icell,size(data,1),istrain,irep,imodel))
                
                %create dat structure for optimization
                DA(1).y = data(icell,:)';
                DA(1).t = (1-1)*3/60:3/60:(size(data,2)-1)*3/60;
                
                if pre_sol == 0 %if do not load pre-computed parameters
                    
                    %set optimization options:
                    optionsPesto = PestoOptions();
                    
                    %do not show estimation window
                    optionsPesto.mode = 'silent';
                    
                    %no gradient used for optimization
                    optionsPesto.localOptimizerOptions.GradObj = 'off';
                    
                    %number of starts
                    optionsPesto.n_starts = 20;
                    
                    %perform optimization
                    scR(icell).sol = getMultiStarts(parameters,@(xi)logLikelihood_WTElp6(xi,DA,imodel),optionsPesto);
                    
                    %calculate Bayesian Information Criterion (BIC) for
                    %specific total GFP trace and model
                    scR(icell).sol.BIC = log(length(DA(1).y))*parameters.number-2*scR(icell).sol.MS.logPost(1);
                    
                    %if optimization did not converge - increase the
                    %starts to 50
                    if abs(scR(icell).sol.MS.logPost(1)-scR(icell).sol.MS.logPost(5)) > 0.1
                        optionsPesto = PestoOptions();
                        optionsPesto.mode = 'silent';
                        optionsPesto.localOptimizerOptions.GradObj = 'off';
                        optionsPesto.n_starts = 50;
                        scR(icell).sol = getMultiStarts(parameters,@(xi)logLikelihood_WTElp6(xi,DA,imodel),optionsPesto);
                        scR(icell).sol.BIC = log(length(DA(1).y))*parameters.number-2*scR(icell).sol.MS.logPost(1);
                    end
                    
                    %if optimization did not converge - increase the
                    %starts to 100
                    if abs(scR(icell).sol.MS.logPost(1)-scR(icell).sol.MS.logPost(5)) > 0.1
                        optionsPesto = PestoOptions();
                        optionsPesto.mode = 'silent';
                        optionsPesto.localOptimizerOptions.GradObj = 'off';
                        optionsPesto.n_starts = 100;
                        scR(icell).sol = getMultiStarts(parameters,@(xi)logLikelihood_WTElp6(xi,DA,imodel),optionsPesto);
                        scR(icell).sol.BIC = log(length(DA(1).y))*parameters.number-2*scR(icell).sol.MS.logPost(1);
                    end
                    
                    %if optimization did not converge - increase the
                    %starts to 200
                    if abs(scR(icell).sol.MS.logPost(1)-scR(icell).sol.MS.logPost(5)) > 0.1
                        optionsPesto = PestoOptions();
                        optionsPesto.mode = 'silent';
                        optionsPesto.localOptimizerOptions.GradObj = 'off';
                        optionsPesto.n_starts = 200;
                        scR(icell).sol = getMultiStarts(parameters,@(xi)logLikelihood_WTElp6(xi,DA,imodel),optionsPesto);
                        scR(icell).sol.BIC = log(length(DA(1).y))*parameters.number-2*scR(icell).sol.MS.logPost(1);
                    end
                    
                else %if load pre-computed parameters
                    if iexp == 1
                        load(sprintf('./Results/scR1_strain%d_rep%d_model%d',istrain,irep,imodel))
                    else
                        load(sprintf('./Results/scR2_strain%d_rep%d_model%d',istrain,irep,imodel))
                    end
                end
                
                if plt == 1
                    
                    %transform the parameters back from log10 space into normal space
                    par = 10.^(scR(icell).sol.MS.par(:,1));
                    
                    if imodel == 1
                        
                        %define the specific parameters for model
                        indA = ind(1:4);
                        P01 = par(indA(1));
                        b1 = par(indA(2));
                        c1 = par(indA(3));
                        sigmayA = par(indA(4))*ones(length(DA.y),1);
                        
                    else
                        
                        %define the specific parameters for model
                        indA = ind(1:5);
                        P01 = par(indA(1));
                        t_rep1 = par(indA(2));
                        b1 = par(indA(3));
                        c1 = par(indA(4));
                        sigmayA = par(indA(5))*ones(length(DA.y),1);
                    end
                    
                    %simulate with given parameter set
                    count = 1;
                    for t = DA(1).t
                        if imodel == 1
                            f1(count) = b1/(c1)+(P01-b1/(c1))*exp(-c1*t);
                        else
                            if t<t_rep1
                                f1(count) = b1/(c1)+(P01-b1/(c1))*exp(-c1*t);
                            else
                                P0_init = b1/(c1)+(P01-b1/(c1))*exp(-c1*t_rep1);
                                f1(count) = P0_init*exp(-c1*(t-t_rep1));
                            end
                        end
                        count = count+1;
                    end
                    f1 = f1';
                    
                    plot(DA(1).t,DA(1).y,'-','Color','k');
                    hold on
                    plot(DA(1).t,f1,':','Color','k');
                    hold on
                    
                    box off
                    set(gca,'linewidth',1.02)
                    set(gca,'FontSize',11)
                    set(gca,'FontName','Arial')
                    xticks([0:1:2])
                    xlabel('time(h)')
                    xlim([0,2])
                    ylabel('GFP intensity')
                end
            end
            
            %save the estimated parameters
            if pre_sol == 0
                if iexp == 1
                    save(sprintf('./Results/scR1_strain%d_rep%d_model%d',istrain,irep,imodel),'scR')
                else
                    save(sprintf('./Results/scR2_strain%d_rep%d_model%d',istrain,irep,imodel),'scR')
                end
            end
        end
    end
end

%count how many of the cells did not converge even after having used 200
%starts (otherwise increase starts and run again?)
count_unconverged = 0;
for icell = 1:length(scR)
    if abs(scR(icell).sol.MS.logPost(1)-scR(icell).sol.MS.logPost(5)) > 0.1
        count_unconverged = count_unconverged+1;
    end
end

disp(sprintf('%d of %d cells did not converge',count_unconverged,length(scR)))

end