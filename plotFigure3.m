%add path
addpath(genpath(pwd))

%% Figure 3B top - example fits

clearvars;
clc;

%load data of computed non-dividing cells
load('NonDividing')
% figure('visible','off');
figure

%istrain = 1 - WT / = 2 - elp6
%irep = 1 - repression 1 / = 2 - repression 2

istrain = 1;

for  irep = 2
    
    %define color according to strain and repression
    if irep == 1
        c = [175,198,233]./255;
    else
        c = [33,68,120]./255;
    end
    
    %load estimated parameter sets
    load(sprintf('scR_strain%d_rep%d_model1',istrain,irep));
    scR1_1 = scR;
    load(sprintf('scR_strain%d_rep%d_model2',istrain,irep));
    scR1_2 = scR;
    
    %for each total GFP trace determine whether the repressor model is
    %required to fit the data (according to the BIC)
    for i = 1:size(scR1_1,2)
        BIC1_1(i) = scR1_1(i).sol.BIC;
    end
    for i = 1:size(scR1_2,2)
        BIC1_2(i) = scR1_2(i).sol.BIC;
    end
    ind1_2 = find(BIC1_2-BIC1_1<-10); %model 2 best
    ind1_1 = find(BIC1_2-BIC1_1>=-10);%model 1 best
    
%     if irep == 1
%         ind_rand = randsample(1:length(NonDividing{istrain}.I5r1),1);
%     else
%         ind_rand = randsample(1:length(NonDividing{istrain}.I5r2),1);
%     end
    
    if irep == 1
        ind_rand = 71;
    else
        ind_rand = 102;
    end
    
    display(sprintf('BIC non-repressor model: %d', BIC1_1(ind_rand)))
    display(sprintf('BIC repressor model: %d', BIC1_2(ind_rand)))
    
    %for each of the 10 randomly chosen cells, plot the total gFP trace
    %and the fit
    for icell = 1
        
        if ismember(ind_rand(icell),ind1_2)
            display('Better fitted by repressor model')
        else
            display('Better fitted by non-repressor model')
        end
        
        %if the total GFP trace is better explained by the repressor
        %model do
        
        if irep == 1
            plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r1(ind_rand(icell),1:40)./1e7,'-','Color',c)
            hold on
        else
            plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r2(ind_rand(icell),1:40)./1e7,'-','Color',c)
            hold on
        end
        
        par = 10.^(scR1_2(ind_rand(icell)).sol.MS.par(:,1));
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
        
        %if the total GFP trace is better explained by the non-repressor
        %model do
        
        if irep == 1
            plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r1(ind_rand(icell),1:40)./1e7,'-','Color',c)
            hold on
        else
            plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r2(ind_rand(icell),1:40)./1e7,'-','Color',c)
            hold on
        end
        
        par = 10.^(scR1_1(ind_rand(icell)).sol.MS.par(:,1));
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
        
        plot((1-1)*3/60:3/60:(40-1)*3/60,f1,'-','Color',[150,150,150]./255);
        
    end
    
    ylabel('GFP intensity')
    xlabel('time (h)')
    xticks([0,1,2])
    box off
    set(gca,'linewidth',1.02)
    set(gca,'FontSize',11)
    set(gca,'FontName','Arial')
    xlim([0,2])
    if irep == 1
        ylim([0,4])
    else
        ylim([0,0.5])
    end
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.5 5])
    if irep == 1
        print('-dpdf','./Figures/Fig3Btopleft','-painters')
    else
        print('-dpdf','./Figures/Fig3Btopright','-painters')
    end
    
end
%% Figure 3B bottom - profile likelihoods

%for irep = 1 and irep = 2
getProfile(1)

getProfile(2)

%% Figure 3C-D - fits and model selection of repressions 1 and 2 and WT and elp6

clearvars;
clc;

%load data of computed non-dividing cells 
load('NonDividing')
figure('visible','off');

%istrain = 1 - WT / = 2 - elp6
%irep = 1 - repression 1 / = 2 - repression 2

for istrain = 2
    for  irep = 2
        
        %define color according to strain and repression
        if istrain == 1
            if irep == 1
                c = [175,198,233]./255;
            else
                c = [33,68,120]./255;
            end
        else
            if irep == 1
                c = [205,135,222]./255;
            else
                c = [67,31,117]./255;
            end
        end
        
        %load estimated parameter sets 
        load(sprintf('scR_strain%d_rep%d_model1',istrain,irep));
        scR1_1 = scR;
        load(sprintf('scR_strain%d_rep%d_model2',istrain,irep));
        scR1_2 = scR;
        
        %for each total GFP trace determine whether the repressor model is
        %required to fit the data (according to the BIC)
        for i = 1:size(scR1_1,2)
            BIC1_1(i) = scR1_1(i).sol.BIC;
        end
        for i = 1:size(scR1_2,2)
            BIC1_2(i) = scR1_2(i).sol.BIC;
        end
        ind1_2 = find(BIC1_2-BIC1_1<-10); %model 2 best
        ind1_1 = find(BIC1_2-BIC1_1>=-10);%model 1 best
        
        %determine 10 randomly sampled cells which will be plotted 
        if irep == 1
            ind_rand = randsample(1:length(NonDividing{istrain}.I5r1),10);
        else
            ind_rand = randsample(1:length(NonDividing{istrain}.I5r2),10);
        end
        
        %for each of the 10 randomly chosen cells, plot the total gFP trace
        %and the fit 
        for icell = 1:10
            
            %if the total GFP trace is better explained by the repressor
            %model do
            if ismember(ind_rand(icell),ind1_2)
                
                if irep == 1
                    plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r1(ind_rand(icell),1:40)./1e7,'-','Color',c)
                    hold on
                else
                    plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r2(ind_rand(icell),1:40)./1e7,'-','Color',c)
                    hold on
                end
                
                par = 10.^(scR1_2(ind_rand(icell)).sol.MS.par(:,1));
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
            
            %if the total GFP trace is better explained by the non-repressor
            %model do    
            else
                
                if irep == 1
                    plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r1(ind_rand(icell),1:40)./1e7,'-','Color',c)
                    hold on
                else
                    plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r2(ind_rand(icell),1:40)./1e7,'-','Color',c)
                    hold on
                end
                
                par = 10.^(scR1_1(ind_rand(icell)).sol.MS.par(:,1));
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
                
                plot((1-1)*3/60:3/60:(40-1)*3/60,f1,'-','Color',[200,200,200]./255);
                
            end
        end
        
        ylabel('GFP intensity')
        xlabel('time (h)')
        xticks([0,1,2])
        box off
        set(gca,'linewidth',1.02)
        set(gca,'FontSize',11)
        set(gca,'FontName','Arial')
        xlim([0,2])
        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.5 5])
        if istrain == 1
            if irep == 1
                print('-dpdf','./Figures/Fig3Cleft','-painters')
            else
                print('-dpdf','./Figures/Fig3Cright','-painters')
            end
        else
            if irep == 1
                print('-dpdf','./Figures/Fig3Dleft','-painters')
            else
                print('-dpdf','./Figures/Fig3Dright','-painters')
            end
        end
        
    end
end

%% Figure 3E-F GFP0 vs selected model 

clearvars;
clc;

%istrain = 1 - WT / = 2 - elp6
%irep = 1 - repression 1 / = 2 - repression 2

%load total GFP traces of computed non-dividing cells
load('NonDividing')
figure('visible','off');

for istrain = 2
    for irep = 2
        
        %load estimated parameters 
        load(sprintf('scR_strain%d_rep%d_model1',istrain,irep));
        scR1_1 = scR;
        load(sprintf('scR_strain%d_rep%d_model2',istrain,irep));
        scR1_2 = scR;
        
        %determine for each total GFP trace whether the repressor model is
        %required to fit the data (according to the BIC)
        for i = 1:size(scR1_1,2)
            BIC1_1(i) = scR1_1(i).sol.BIC;
        end
        for i = 1:size(scR1_2,2)
            BIC1_2(i) = scR1_2(i).sol.BIC;
        end
        ind1_2 = find(BIC1_2-BIC1_1<-10); %model 2 best
        ind1_1 = find(BIC1_2-BIC1_1>=-10);%model 1 best

        %extract all estimated parameter sets per total GFP trace
        for icell = 1:length(scR1_1)
            clear par
            par = 10.^(scR1_1(icell).sol.MS.par(:,1));
            Par1(icell,:) = par';
        end
        
        for icell = 1:length(scR1_2)
            clear par
            par = 10.^(scR1_2(icell).sol.MS.par(:,1));
            Par2(icell,:) = par';
        end
        
        for ipar = 1

            %only consider total GFP traces which reuiqre the repressor
            %model 
            P1 = Par1(ind1_1,ipar)';
            P2 = Par1(ind1_2,ipar)';
            
            %show the number of cells better fitted by the non-repressor
            %model (1) and repressor model (2) 
            length(P1)
            length(P2)
            
            %show fractions of cells better fitted by the non-repressor
            %model (1) and repressor model (2) 
            length(P1)/(length(P1)+length(P2))
            length(P2)/(length(P1)+length(P2))
            
            index = [ones(length(P1),1);2*ones(length(P2),1)];
            
            %Mood's median test
            [p_mediantest,tab,chi2] = mediantest(P1,P2);
            Pval_mediantest(ipar) = p_mediantest;

            %plot with jitter 
            a = -0.2;
            b = 0.2;
            r1 = (b-a).*rand(length(P1),1) + a;
            plot(1+r1,P1,'.','Color',[200,200,200]./255,'Markersize',10)
            hold on
            
            r2 = (b-a).*rand(length(P2),1) + a;
            plot(2+r2,P2,'.','Color',[0,0,0]./255,'Markersize',10)
            hold on

            line([0.6,1.4],[median(P1),median(P1)],'Color','k','Linewidth',2)
            hold on
            line([1.6,2.4],[median(P2),median(P2)],'Color','k','Linewidth',2)
            
            hold on
            xlim([0,3])
            set(gca,'FontSize',10)
            ylim([0 inf])
            box off
            set(gca,'linewidth',1.02)
            set(gca,'FontSize',11)
            set(gca,'FontName','Arial')
            ylabel('GFP intensity')
            xlabel('model')
        end
        
        %save figure 
        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.5 5])
        if istrain == 1
            if irep == 1
                print('-dpdf','./Figures/Fig3Eleft','-painters')
            else
                print('-dpdf','./Figures/Fig3Eright','-painters')
            end
        else
            if irep == 1
                print('-dpdf','./Figures/Fig3Fleft','-painters')
            else
                print('-dpdf','./Figures/Fig3Fright','-painters')
            end
        end
    end
end

