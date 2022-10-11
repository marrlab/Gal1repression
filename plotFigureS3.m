clearvars;
clc;

%load data of computed non-dividing cells
load('simData')
data = simData/1e2;

%figure('visible','off');

%istrain = 1 - WT / = 2 - elp6
%irep = 1 - repression 1 / = 2 - repression 2
istrain = 1;

figure

for  icell = 1:4
    
    subplot(1,4,icell)
    
    %load estimated parameter sets
    load('scR1_model1_sim');
    scR1_1 = scR;
    load('scR1_model2_sim');
    scR1_2 = scR;
    
    for imodel = 1:2
        
        %if the total GFP trace is better explained by the repressor
        %model do
        if imodel == 1
            plot((1-1)*3/60:3/60:(40-1)*3/60,data(icell,:),':','Color','k')
            hold on
        end
        
        par = 10.^(scR1_2(icell).sol.MS.par(:,1));
        indA = 1:5;
        P01 = par(indA(1));
        t_rep1 = par(indA(2));
        b1 = par(indA(3));
        c1 = par(indA(4));
        sigmayA = par(indA(5))*ones(40,1);
        
        %WT simulation
        count = 1;
        for t = (1-1)*3/60:3/60:(40-1)*3/60
            if t<t_rep1
                f1(count) = b1/(c1)+(P01-b1/(c1))*exp(-c1*t);
            else
                P0_init = b1/(c1)+(P01-b1/(c1))*exp(-c1*t_rep1);
                f1(count) = P0_init*exp(-c1*(t-t_rep1));
            end
            count = count+1;
        end
        f1 = f1';
        plot((1-1)*3/60:3/60:(40-1)*3/60,f1,'-','Color','k');
        
        par = 10.^(scR1_1(icell).sol.MS.par(:,1));
        indA = 1:4;
        P01 = par(indA(1));
        b1 = par(indA(2));
        c1 = par(indA(3));
        sigmayA = par(indA(4))*ones(40,1);
        
        %WT simulation
        count = 1;
        for t = (1-1)*3/60:3/60:(40-1)*3/60
            f1(count) = b1/(c1)+(P01-b1/(c1))*exp(-c1*t);
            count = count+1;
        end
        f1 = f1';
        
        plot((1-1)*3/60:3/60:(40-1)*3/60,f1,'-','Color','r');
        
    end
    
    ylabel('total GFP (a.u)')
    xlabel('repression time (h)')
    xticks([0,1,2])
    box off
    set(gca,'linewidth',1.02)
    set(gca,'FontSize',11)
    set(gca,'FontName','Arial')
    xlim([0,2])
    ylim([0,1.8])
    
end

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 26 4])
print('-dpdf','./Figures/FigS3','-painters')
