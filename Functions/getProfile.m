
%model = 1 or 2 / model == 1 - non-repressor model / model == 2 - repressor model
%rep = 1 or 2 / rep == 1 - repression 1 / rep == 2 - repression 2
%strain = 1 or 2 / strain == 1 - WT / strain == 2 - elp6
%plt = 0 or 1 / plt == 0 - no plots / plt == 1 - plots (data and fits)
%server = 0 or 1 / server == 0 - not on server / server == 1 - on server
%pre_sol = 0 or 1 / pre_sol == 0 - do not load precomputed parameter estimates / pre_sol == 1 - load precomputed parameter estimates

model = 2;
rep = 1;
strain = 1;
plt = false;
server = false;
pre_sol = false;

clearvars -except model rep strain plt server pre_sol
clc;

%add path for AMICI and PESTO
if server == 0 %adapt paths
    addpath(genpath('/Users/lea.schuh/Documents/PhD/ICB/Xenopus/PESTO-master'))
    addpath(genpath('/Users/lea.schuh/Documents/PhD/ICB/Xenopus/AMICI-master'))
    addpath(genpath(pwd))
    load('./Data/NonDividing')
else %adapt paths
    addpath(genpath('/home/icb/lea.schuh/AMICI-master'))
    addpath(genpath('/home/icb/lea.schuh/PESTO-master'))
    addpath(genpath(pwd))
    load('/home/icb/lea.schuh/Gal1/NonDividing')
end

for imodel = model
    for irep = rep
        for istrain = strain
            
            clear scR DA data
            
            if irep == 1
                %normalize data by 1e7
                data = NonDividing{istrain}.I5r1./1e7;
            else
                data = NonDividing{istrain}.I5r2./1e7;
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
            
            if irep == 1
                Icell = 71;
            else
                Icell = 102;
            end
            
            for icell = Icell
                
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
                    optionsPesto.n_starts = 50;
                    
                    %perform optimization
                    scR(icell).sol = getMultiStarts(parameters,@(xi)logLikelihood_WTElp6(xi,DA,imodel),optionsPesto);
                    
                    %calculate Bayesian Information Criterion (BIC) for
                    %specific total GFP trace and model
                    scR(icell).sol.BIC = log(length(DA(1).y))*parameters.number-2*scR(icell).sol.MS.logPost(1);

                else %if load pre-computed parameters
                    load(sprintf('./Data/scR_strain%d_rep%d_model%d',istrain,irep,imodel))  
                end
                
                parameters = getParameterProfiles(scR(icell).sol,@(xi)logLikelihood_WTElp6(xi,DA,imodel),optionsPesto);
                
                figure
                for i = 1:length(parameters.name)-1
                    subplot(1,length(parameters.name)-1,i)
                    plot(10.^(parameters.P(i).par(i,:)),parameters.P(i).logPost,'-','Color','k')
                    hold on
                    plot(10.^(scR(icell).sol.MS.par(i,1)),scR(icell).sol.MS.logPost(1),'*','Color','k')
                    
                    set(gca,'FontSize',10)
                    box off
                    set(gca,'linewidth',1.02)
                    set(gca,'FontSize',11)
                    set(gca,'FontName','Arial')
                    if irep == 1
                        ylim([28.5,31.5])
                    else
                        ylim([110,113])
                    end
                end
                
                %save figure
                set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 13 4])
                if irep == 1
                    print('-dpdf','./Figures/Fig3Cnewleft','-painters')
                else
                    print('-dpdf','./Figures/Fig3Cnewright','-painters')
                end
                
                if plt == 1
                    
                    figure
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
            
        end
    end
end

