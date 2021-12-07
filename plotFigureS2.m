%add path
addpath(genpath(pwd))

%% Figure S2A - full data set elp6 - replicate
clearvars;
clc;

% figure('visible','off');
figure
for istrain = 2
    
    clearvars -except istrain
    
    if istrain == 1
        
        clear Data
        %WT
        DataExpI5{1} = sprintf('ExpI5_pos1_Y208_WT');
        DataExpI5{2} = sprintf('ExpI5_pos2_Y208_WT');
        DataExpI5{3} = sprintf('ExpI5_pos4_Y208_WT');
        DataExpI5{4} = sprintf('ExpI5_pos5_Y208_WT');
        DataExpI5{6} = sprintf('ExpI5_pos8_Y208_WT');
        DataExpI5{7} = sprintf('ExpI5_pos9_Y208_WT');
        DataExpI5{8} = sprintf('ExpI5_pos10_Y208_WT');
        DataExpI5{9} = sprintf('ExpI5_pos11_Y208_WT');
        DataExpI5{10} = sprintf('ExpI5_pos12_Y208_WT');
        DataExpI5{11} = sprintf('ExpI5_pos13_Y208_WT');
        DataExpI5{12} = sprintf('ExpI5_pos14_Y208_WT');
        DataExpI5{13} = sprintf('ExpI5_pos15_Y208_WT');
        DataExpI5{14} = sprintf('ExpI5_pos16_Y208_WT');
        
    else
        
        clear Data
        %Elp6
%         DataExpI5{15} = sprintf('ExpI5_pos19_Y1474_Elp6');
%         DataExpI5{16} = sprintf('ExpI5_pos20_Y1474_Elp6');
%         DataExpI5{17} = sprintf('ExpI5_pos22_Y1474_Elp6');
%         DataExpI5{18} = sprintf('ExpI5_pos23_Y1474_Elp6');
%         DataExpI5{19} = sprintf('ExpI5_pos24_Y1474_Elp6');
%         DataExpI5{20} = sprintf('ExpI5_pos25_Y1474_Elp6');
%         DataExpI5{21} = sprintf('ExpI5_pos27_Y1474_Elp6');
%         DataExpI5{22} = sprintf('ExpI5_pos30_Y1474_Elp6');
%         DataExpI5{23} = sprintf('ExpI5_pos31_Y1474_Elp6');
%         DataExpI5{24} = sprintf('ExpI5_pos32_Y1474_Elp6');
        
        DataExpL{15} = sprintf('ExpL1_pos9_Y1474_Elp6');
        DataExpL{16} = sprintf('ExpL1_pos10_Y1474_Elp6');
        DataExpL{17} = sprintf('ExpL1_pos11_Y1474_Elp6');
        DataExpL{18} = sprintf('ExpL1_pos12_Y1474_Elp6');
        DataExpL{19} = sprintf('ExpL1_pos13_Y1474_Elp6');
        DataExpL{20} = sprintf('ExpL1_pos14_Y1474_Elp6');
        DataExpL{21} = sprintf('ExpL1_pos15_Y1474_Elp6');
        DataExpL{22} = sprintf('ExpL1_pos16_Y1474_Elp6');
        
        DataExpL{23} = sprintf('ExpL2_pos14_Y1474_Elp6');
        DataExpL{24} = sprintf('ExpL2_pos15_Y1474_Elp6');
        DataExpL{25} = sprintf('ExpL2_pos16_Y1474_Elp6');
        DataExpL{26} = sprintf('ExpL2_pos17_Y1474_Elp6');
        DataExpL{27} = sprintf('ExpL2_pos18_Y1474_Elp6');
        DataExpL{28} = sprintf('ExpL2_pos19_Y1474_Elp6');
        DataExpL{29} = sprintf('ExpL2_pos20_Y1474_Elp6');
        DataExpL{30} = sprintf('ExpL2_pos21_Y1474_Elp6');
        DataExpL{31} = sprintf('ExpL2_pos22_Y1474_Elp6');
        DataExpL{32} = sprintf('ExpL2_pos23_Y1474_Elp6');
        DataExpL{33} = sprintf('ExpL2_pos24_Y1474_Elp6');
        DataExpL{34} = sprintf('ExpL2_pos25_Y1474_Elp6');
        DataExpL{35} = sprintf('ExpL2_pos26_Y1474_Elp6');
        DataExpL{36} = sprintf('ExpL2_pos27_Y1474_Elp6');
        DataExpL{37} = sprintf('ExpL2_pos28_Y1474_Elp6');
    end
    
%     Data = DataExpI5;
    Data = DataExpL;
    
    Data = Data(~cellfun('isempty',Data));
    count = 0;
    
    for i = 1:length(Data)
        
        clearvars -except Data i istrain
        
        %load data of position
        loadData = sprintf('S%s',Data{i});
        load(loadData);
        
        %for every cell in position do
        for iS = 1:length(S)
            
            %color definition
            c = [221,175,233;160,160,160;203,133,221;160,160,160;66,30,115]./255;
            
            %repression 0
            ind81 = find(S{iS}.DF:S{iS}.LF==81);
            if isempty(ind81) == 0
                plot((81-length(S{iS}.GFPabs(1:ind81)))*3/60:3/60:(81-1)*3/60,S{iS}.GFPabs(1:ind81)./10e6,'-','Color',c(1,:))
                hold on
            else
                ind81 = 1;
            end
            
            %induction 1
            ind141 = find(S{iS}.DF:S{iS}.LF==141);
            if isempty(ind141) == 0
                plot((141-length(S{iS}.GFPabs(ind81:ind141)))*3/60:3/60:(141-1)*3/60,S{iS}.GFPabs(ind81:ind141)./10e6,'-','Color',c(2,:))
                hold on
            else
                ind141 = 1;
            end
            
            %repression 1
            ind221 = find(S{iS}.DF:S{iS}.LF==221);
            if isempty(ind221) == 0
                plot((221-length(S{iS}.GFPabs(ind141:ind221)))*3/60:3/60:(221-1)*3/60,S{iS}.GFPabs(ind141:ind221)./10e6,'-','Color',c(3,:))
                hold on
            else
                ind221 = 1;
            end
            
            %induction 2
            ind281 = find(S{iS}.DF:S{iS}.LF==281);
            if isempty(ind281) == 0
                plot((281-length(S{iS}.GFPabs(ind221:ind281)))*3/60:3/60:(281-1)*3/60,S{iS}.GFPabs(ind221:ind281)./10e6,'-','Color',c(4,:))
                hold on
            else
                ind281 = 1;
            end
            
            %repression 2
            plot((S{iS}.LF-length(S{iS}.GFPabs(ind281:end)))*3/60:3/60:(S{iS}.LF-1)*3/60,S{iS}.GFPabs(ind281:end)./10e6,'-','Color',c(5,:))
            hold on
        end
        
    end
end
ylabel('GFP intensity')
yticks([0,2,4])
xlabel('time (h)')
xticks([0,4,7,11,14,16])
box off
set(gca,'linewidth',1.02)
set(gca,'FontSize',11)
set(gca,'FontName','Arial')
ylim([0,4])
xlim([0,16])
set(gcf, 'DefaultFigureRenderer', 'painters');
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.5 5])
print('-dpdf','./Figures/FigS2A','-painters')

%% Figure S2B - non-dividing total GFP traces for repressions r1 and r2 for elp6 - replicate

clearvars;
clc;

load('NonDividing2')
% figure('visible','off');

%istrain = 1 - WT / = 2 - elp6
%irep = 1 - repression 1 / = 2 - repression 2

for istrain = 2
    for irep = 1:2
        
        figure
        
        %define color according to strain and repression
        if irep == 1
            c = [203,133,221]./255;
        else
            c = [66,30,115]./255;
        end
        
        %plot the non-dividing cells of specified strain and repression 
        %also plot mean total GFP trace and maximal mean total GFP value
        if irep == 1
            plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r1(:,1:40)./10e6,'-','Color',c)
            hold on
            plot((1-1)*3/60:3/60:(40-1)*3/60,mean(NonDividing{istrain}.I5r1(:,1:40))./10e6,'-','Color','k','Linewidth',2)
            hold on
            indmax = find(mean(NonDividing{istrain}.I5r1(:,1:40))==max(mean(NonDividing{istrain}.I5r1(:,1:40))));
            time = (1-1)*3/60:3/60:(40-1)*3/60;
            plot(time(indmax),max(mean(NonDividing{istrain}.I5r1(:,1:40)))./10e6,'.','Color','k','Markersize',12)
            
            %bootstrap the time to maximal mean total GFP 
            for isample = 1:100000
                S = datasample(1:size(NonDividing{istrain}.I5r1,1),round(size(NonDividing{istrain}.I5r1,1)));
                indmax = find(mean(NonDividing{istrain}.I5r1(S,1:40))==max(mean(NonDividing{istrain}.I5r1(S,1:40))));
                T(isample) = time(indmax(1));
            end
            
            %get mean and std of time to maximal mean total GFP 
            display(sprintf('Mean of elp6 time maximal mean total GFP repression r1 is %d', mean(T)))
            display(sprintf('Standard deviation of elp6 time maximal mean total GFP repression r1 is %d', std(T)))
            display(sprintf('Number of elp6 cells repression r1 is %d', size(NonDividing{istrain}.I5r1,1)))
            
            line([mean(T),mean(T)],[max(mean(NonDividing{istrain}.I5r1(:,1:40)))./10e6-0.4,max(mean(NonDividing{istrain}.I5r1(:,1:40)))./10e6+0.4],'Color','k','Linewidth',1)
            hold on
            line([mean(T)-std(T),mean(T)+std(T)],[max(mean(NonDividing{istrain}.I5r1(:,1:40)))./10e6,max(mean(NonDividing{istrain}.I5r1(:,1:40)))./10e6],'Color','k','Linewidth',1)
            hold on
            line([mean(T)+std(T),mean(T)+std(T)],[max(mean(NonDividing{istrain}.I5r1(:,1:40)))./10e6-0.2,max(mean(NonDividing{istrain}.I5r1(:,1:40)))./10e6+0.2],'Color','k','Linewidth',1)
            hold on
            line([mean(T)-std(T),mean(T)-std(T)],[max(mean(NonDividing{istrain}.I5r1(:,1:40)))./10e6-0.2,max(mean(NonDividing{istrain}.I5r1(:,1:40)))./10e6+0.2],'Color','k','Linewidth',1)
    
        else
            plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r2(:,1:40)./10e6,'-','Color',c)
            hold on
            plot((1-1)*3/60:3/60:(40-1)*3/60,mean(NonDividing{istrain}.I5r2(:,1:40))./10e6,'-','Color','w','Linewidth',2)
            hold on
            indmax = find(mean(NonDividing{istrain}.I5r2(:,1:40))==max(mean(NonDividing{istrain}.I5r2(:,1:40))));
            time = (1-1)*3/60:3/60:(40-1)*3/60;
            plot(time(indmax),max(mean(NonDividing{istrain}.I5r2(:,1:40)))./10e6,'.','Color','w','Markersize',12)
            
            %bootstrap the time to maximal mean total GFP 
            for isample = 1:100000
                S = datasample(1:size(NonDividing{istrain}.I5r2,1),round(size(NonDividing{istrain}.I5r2,1)));
                indmax = find(mean(NonDividing{istrain}.I5r2(S,1:40))==max(mean(NonDividing{istrain}.I5r2(S,1:40))));
                T(isample) = time(indmax);
            end
            
            display(sprintf('Mean of elp6 time maximal mean total GFP repression r2 is %d', mean(T)))
            display(sprintf('Standard deviation of elp6 time maximal mean total GFP repression r2 is %d', std(T)))
            display(sprintf('Number of elp6 cells repression r2 is %d', size(NonDividing{istrain}.I5r2,1)))
           
            line([mean(T),mean(T)],[max(mean(NonDividing{istrain}.I5r2(:,1:40)))./10e6-0.4,max(mean(NonDividing{istrain}.I5r2(:,1:40)))./10e6+0.4],'Color','w','Linewidth',1)
            hold on
            line([mean(T)-std(T),mean(T)+std(T)],[max(mean(NonDividing{istrain}.I5r2(:,1:40)))./10e6,max(mean(NonDividing{istrain}.I5r2(:,1:40)))./10e6],'Color','w','Linewidth',1)
            hold on
            line([mean(T)+std(T),mean(T)+std(T)],[max(mean(NonDividing{istrain}.I5r2(:,1:40)))./10e6-0.2,max(mean(NonDividing{istrain}.I5r2(:,1:40)))./10e6+0.2],'Color','w','Linewidth',1)
            hold on
            line([mean(T)-std(T),mean(T)-std(T)],[max(mean(NonDividing{istrain}.I5r2(:,1:40)))./10e6-0.2,max(mean(NonDividing{istrain}.I5r2(:,1:40)))./10e6+0.2],'Color','w','Linewidth',1)     

        end
        
        ylabel('GFP intensity')
        if irep == 1
            yticks([0,1,2])
            ylim([0,2.5])
        else
            yticks([0,2,4])
            ylim([0,4])
        end
        xlabel('time (h)')
        xticks([0,2,4])
        box off
        set(gca,'linewidth',1.02)
        set(gca,'FontSize',11)
        set(gca,'FontName','Arial')
        xlim([0,2])
        set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7 5])
        if irep == 1
            print('-dpdf','./Figures/FigS2Bleft','-painters')
        else
            print('-dpdf','./Figures/FigS2Bright','-painters')
        end
    end
end


%% Figure S2C - random example fits and model selection of repressions r1 and r2 elp6 - replicate

%istrain = 1 - WT / = 2 - elp6
%irep = 1 - repression 1 / = 2 - repression 2

for istrain = 2
    for  irep = 1:2
        
        clearvars -except istrain irep
        clc;
        
        %load data of computed non-dividing cells
        load('NonDividing2')
%         figure('visible','off');
        figure
        
        %define color according to strain and repression
        if irep == 1
            c = [203,133,221]./255;
        else
            c = [66,30,115]./255;
        end
        
        %load estimated parameter sets
        load(sprintf('scR2_strain%d_rep%d_model1_truncated',istrain,irep));
        scR1_1 = scR;
        load(sprintf('scR2_strain%d_rep%d_model2_truncated',istrain,irep));
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
                    plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r1(ind_rand(icell),1:40)./1e7,':','Color',c)
                    hold on
                else
                    plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r2(ind_rand(icell),1:40)./1e7,':','Color',c)
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
                    plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r1(ind_rand(icell),1:40)./1e7,':','Color',c)
                    hold on
                else
                    plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r2(ind_rand(icell),1:40)./1e7,':','Color',c)
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
                
                plot((1-1)*3/60:3/60:(40-1)*3/60,f1,'-','Color','r');
                
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
        if irep == 1
            print('-dpdf','./Figures/FigS2Cleft','-painters')
        else
            print('-dpdf','./Figures/FigS2Cright','-painters')
        end
    end
end

%% Figure S2D - GFP0 vs selected model for repressions r1 and r2 elp6 - replicate

%istrain = 1 - WT / = 2 - elp6
%irep = 1 - repression 1 / = 2 - repression 2
clc;

for istrain = 2
    for irep = 1:2
        
        clearvars -except istrain irep
        %         clc;
        
        %load total GFP traces of computed non-dividing cells
%         figure('visible','off');
        load('NonDividing2')
        figure
        
        %load estimated parameters
        load(sprintf('scR2_strain%d_rep%d_model1_truncated',istrain,irep));
        scR1_1 = scR;
        load(sprintf('scR2_strain%d_rep%d_model2_truncated',istrain,irep));
        scR1_2 = scR;
%         load(sprintf('scR_strain%d_rep%d_model1',istrain,irep));
%         scR1_1 = scR;
%         load(sprintf('scR_strain%d_rep%d_model2',istrain,irep));
%         scR1_2 = scR;
        
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
            display(sprintf('%d elp6 cells of repression r%d are better fitted by a non-repressor model',length(P1),irep))
            display(sprintf('%d elp6 cells of repression r%d are better fitted by a repressor model',length(P2),irep))
            
            %show fractions of cells better fitted by the non-repressor
            %model (1) and repressor model (2)
            display(sprintf('%d percent of elp6 cells repression r%d are better fitted by a non-repressor model',length(P1)/(length(P1)+length(P2))*100,irep))
            display(sprintf('%d percent of elp6 cells repression r%d are better fitted by a repressor model',length(P2)/(length(P1)+length(P2))*100,irep))
            
            index = [ones(length(P1),1);2*ones(length(P2),1)];
            
            %Mood's median test
            [p_mediantest,tab,chi2] = mediantest(P1,P2);
            Pval_mediantest(ipar) = p_mediantest;
            
            Pval_mediantest
            
            %plot with jitter
            a = -0.2;
            b = 0.2;
            r1 = (b-a).*rand(length(P1),1) + a;
            plot(1+r1,P1,'.','Color','r','Markersize',10)
            hold on
            
            r2 = (b-a).*rand(length(P2),1) + a;
            plot(2+r2,P2,'.','Color','k','Markersize',10)
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
        if irep == 1
            print('-dpdf','./Figures/FigS2Dleft','-painters')
        else
            print('-dpdf','./Figures/FigS2Dright','-painters')
        end
    end
end

%% Figure S2E - time to maximal mean total GFP for elp6 repressor cells of repressions r1 and r2 - replicate

clearvars;
clc;

% load('NonDividing')
load('NonDividing2')
count = 1;

for strain1 = 2
    for rep1 = 1:2
        
        clear ind1_2 BIC1_1 BIC1_2
        
        %istrain = 1 - WT / = 2 - elp6
        %irep = 1 - repression 1 / = 2 - repression 2
        
        %load all estimated parameter sets for both models and repressions
%         load(sprintf('scR_strain%d_rep%d_model%d',strain1,rep1,1))
%         scR1_1 = scR;
%         load(sprintf('scR_strain%d_rep%d_model%d',strain1,rep1,2))
%         scR1_2 = scR;
        load(sprintf('scR2_strain%d_rep%d_model%d_truncated',strain1,rep1,1))
        scR1_1 = scR;
        load(sprintf('scR2_strain%d_rep%d_model%d_truncated',strain1,rep1,2))
        scR1_2 = scR;
        
        %extract BIC values for all single-cell trajectories for data set 1
        for i = 1:size(scR1_1,2)
            BIC1_1(i) = scR1_1(i).sol.BIC;
        end
        for i = 1:size(scR1_2,2)
            BIC1_2(i) = scR1_2(i).sol.BIC;
        end
        
        %decide whether single-cell requires repressor model or not for data set 1
        ind1_2 = find(BIC1_2-BIC1_1<-10); %model 2 best
        ind1_1 = find(BIC1_2-BIC1_1>=-10);%model 1 best
        
        istrain = strain1;
        irep = rep1;
        
        %plot the non-dividing cells of specified strain and repression
        %also plot mean total GFP trace and maximal mean total GFP value
%         time = (1-1)*3/60:3/60:(40-1)*3/60;
        time = (1-1)*3/60:3/60:(81-1)*3/60;
        if irep == 1
            for isample = 1:100000
                S = datasample(1:length(ind1_2),length(ind1_2));
                indmax = find(mean(NonDividing{istrain}.I5r1(ind1_2(S),1:40))==max(mean(NonDividing{istrain}.I5r1(ind1_2(S),1:40))));
                T(count,isample) = time(indmax);
            end
            
            display(sprintf('Mean of WT time maximal mean total GFP repression r1 is %d for repressor cells', mean(T(count,:))))
            display(sprintf('Standard deviation of WT time maximal mean total GFP repression r1 is %d for repressor cells', std(T(count,:))))
            display(sprintf('Number of WT repressor cells in repression r1 is %d', length(ind1_2)))
            
        else
            for isample = 1:100000
                S = datasample(1:length(ind1_2),length(ind1_2));
                indmax = find(mean(NonDividing{istrain}.I5r2(ind1_2(S),1:40))==max(mean(NonDividing{istrain}.I5r2(ind1_2(S),1:40))));
                T(count,isample) = time(indmax);
            end
            
            display(sprintf('Mean of WT time maximal mean total GFP repression r2 is %d for repressor cells', mean(T(count,:))))
            display(sprintf('Standard deviation of WT time maximal mean total GFP repression r2 is %d for repressor cells', std(T(count,:))))
            display(sprintf('Number of WT repressor cells in repression r2 is %d', length(ind1_2)))
            
        end
        
        count = count+1;
        
    end
end

% figure('visible','off');
figure
for icount = 1:size(T,1)
    
    if icount == 1
        c = [203,133,221]./255;
    else
        c = [66,30,115]./255;
    end
    
    line([mean(T(icount,:)),mean(T(icount,:))],[icount-0.4,icount+0.4],'Color',c,'Linewidth',1)
    hold on
    line([mean(T(icount,:))-std(T(icount,:)),mean(T(icount,:))+std(T(icount,:))],[icount,icount],'Color',c,'Linewidth',1)
    hold on
    line([mean(T(icount,:))+std(T(icount,:)),mean(T(icount,:))+std(T(icount,:))],[icount-0.2,icount+0.2],'Color',c,'Linewidth',1)
    hold on
    line([mean(T(icount,:))-std(T(icount,:)),mean(T(icount,:))-std(T(icount,:))],[icount-0.2,icount+0.2],'Color',c,'Linewidth',1)
    
end

% xlim([0,5])
set(gca,'FontSize',10)
% xlim([0, 2])
xlim([0, 2])
box off
set(gca,'linewidth',1.02)
set(gca,'FontSize',11)
set(gca,'FontName','Arial')

%save figure
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.5 5])
print('-dpdf','./Figures/FigS2E','-painters')

%% Figure S2F - comparison of estimated parameters elp6 - replicate

%strain = 1 - WT / = 2 - elp6
%rep = 1 - repression 1 / = 2 - repression 2

strain1 = 2;
rep1 = 1;
strain2 = 2;
rep2 = 2;
paired = true;

clearvars -except rep1 strain1 rep2 strain2 paired
clc;

%load all estimated parameter sets for both models and repressions
load(sprintf('scR2_strain%d_rep%d_model%d_truncated',strain1,rep1,1))
scR1_1 = scR;
load(sprintf('scR2_strain%d_rep%d_model%d_truncated',strain1,rep1,2))
scR1_2 = scR;

load(sprintf('scR2_strain%d_rep%d_model%d_truncated',strain2,rep2,1))
scR2_1 = scR;
load(sprintf('scR2_strain%d_rep%d_model%d_truncated',strain2,rep2,2))
scR2_2 = scR;

%extract BIC values for all single-cell trajectories for data set 1
for i = 1:size(scR1_1,2)
    BIC1_1(i) = scR1_1(i).sol.BIC;
end
for i = 1:size(scR1_2,2)
    BIC1_2(i) = scR1_2(i).sol.BIC;
end

%decide whether single-cell requires repressor model or not for data set 1
ind1_2 = find(BIC1_2-BIC1_1<-10); %model 2 best
ind1_1 = find(BIC1_2-BIC1_1>=-10);%model 1 best

%extract BIC values for all single-cell trajectories for data set 2
for i = 1:size(scR2_1,2)
    BIC2_1(i) = scR2_1(i).sol.BIC;
end
for i = 1:size(scR2_2,2)
    BIC2_2(i) = scR2_2(i).sol.BIC;
end

%decide whether single-cell requires repressor model or not for data set 2
ind2_2 = find(BIC2_2-BIC2_1<-10);
ind2_1 = find(BIC2_2-BIC2_1>=-10);

%get the data sets for both experiments
load('NonDividing2')

%extract data from non-dividing cells structure
if rep1 == 1
    data1 = NonDividing{strain1}.I5r1;
    momID1 = NonDividing{strain1}.I5momIDr1;
    mompos1 = NonDividing{strain1}.I5momposr1;
else
    data1 = NonDividing{strain1}.I5r2;
    momID1 = NonDividing{strain1}.I5momIDr2;
    mompos1 = NonDividing{strain1}.I5momposr2;
end

if rep2 == 1
    data2 = NonDividing{strain2}.I5r1;
    momID2 = NonDividing{strain2}.I5momIDr1;
    mompos2 = NonDividing{strain2}.I5momposr1;
else
    data2 = NonDividing{strain2}.I5r2;
    momID2 = NonDividing{strain2}.I5momIDr2;
    mompos2 = NonDividing{strain2}.I5momposr2;
end

%extract the estimated parameter sets per cell and data set
for icell = 1:length(scR1_2)
    clear par
    par = 10.^(scR1_2(icell).sol.MS.par(:,1));
    Par1(icell,:) = par';
end

for icell = 1:length(scR2_2)
    clear par
    par = 10.^(scR2_2(icell).sol.MS.par(:,1));
    Par2(icell,:) = par';
end

%reduce parameter sets to the cells requiring the repressor model
Par1 = Par1(ind1_2,:);
Par2 = Par2(ind2_2,:);

%if comparison between cells in repressions 1 and 2
if paired == 1
    %extract the cell IDs and positions of cells requiring repressor model in
    %repressions 1 and 2
    momInfo1 = [momID1(ind1_2)',mompos1(ind1_2)'];
    momInfo2 = [momID2(ind2_2)',mompos2(ind2_2)'];
    [~,index_momInfo1,index_momInfo2] = intersect(momInfo1,momInfo2,'rows');
    Par1 = Par1(index_momInfo1,:);
    Par2 = Par2(index_momInfo2,:);
end

%determine the color according to strain and repression
if strain1 == 1
    if rep1 == 1
        c1 = [175,198,233]./255;
    else
        c1 = [33,68,120]./255;
    end
else
    if rep1 == 1
        c1 = [205,135,222]./255;
    else
        c1 = [67,31,117]./255;
    end
end

if strain2 == 1
    if rep2 == 1
        c2 = [175,198,233]./255;
    else
        c2 = [33,68,120]./255;
    end
else
    if rep2 == 1
        c2 = [205,135,222]./255;
    else
        c2 = [67,31,117]./255;
    end
end

% figure('visible','off');
figure

%for each of the estimated parameters GFP0, delay, rprod and rdeg do
for ipar = 1:4
    
    clear P1_1 index h_kstest2 p_kstest2 h_ranksum p_ranksum
    
    subplot(1,4,ipar)
    
    %get correct estimnated parameters
    P1 = Par1(:,ipar)';
    P2 = Par2(:,ipar)';
    index = [ones(length(P1),1);2*ones(length(P2),1)];
    
    %plot with jitter
    a = -0.2;
    b = 0.2;
    r1 = (b-a).*rand(length(P1),1) + a;
    
    if paired
        for i = 1:size(P1,2)
            line([(1+r1(i)),(2+r1(i))],[P1(i),P2(i)],'Color',[200,200,200]./255);
            hold on
        end
        r2 = r1;
    else
        r2 = (b-a).*rand(length(P1),1) + a;
    end
    
    sum(P1>P2)/length(P1)*100
    
    plot(1+r1,P1,'.','Color',c1,'Markersize',10)
    hold on
    plot(2+r2,P2,'.','Color',c2,'Markersize',10)
    hold on
    line([0.6,1.4],[median(P1),median(P1)],'Color','k','Linewidth',2)
    hold on
    line([1.6,2.4],[median(P2),median(P2)],'Color','k','Linewidth',2)
    hold on
    
    xlim([0,3])
    set(gca,'FontSize',10,'XTickLabelRotation',90)
    ylim([0 inf])
    box off
    set(gca,'linewidth',1.02)
    set(gca,'FontSize',11)
    set(gca,'FontName','Arial')
    
end

%save figureNonD
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 13 4])
print('-dpdf','./Figures/FigS2F','-painters')

%% Figure S2G - GFP0 vs tdelay for elp6 - replicate

%strain = 1 - WT / = 2 - elp6
%rep = 1 - repression 1 / = 2 - repression 2

strain1 = 2;
rep1 = 1;
strain2 = 2;
rep2 = 2;

clearvars -except rep1 strain1 rep2 strain2 paired
clc;

%load all estimated parameter sets for both models and repressions
load(sprintf('scR2_strain%d_rep%d_model%d_truncated',strain1,rep1,1))
scR1_1 = scR;
load(sprintf('scR2_strain%d_rep%d_model%d_truncated',strain1,rep1,2))
scR1_2 = scR;

load(sprintf('scR2_strain%d_rep%d_model%d_truncated',strain2,rep2,1))
scR2_1 = scR;
load(sprintf('scR2_strain%d_rep%d_model%d_truncated',strain2,rep2,2))
scR2_2 = scR;

%extract BIC values for all single-cell trajectories for data set 1
for i = 1:size(scR1_1,2)
    BIC1_1(i) = scR1_1(i).sol.BIC;
end
for i = 1:size(scR1_2,2)
    BIC1_2(i) = scR1_2(i).sol.BIC;
end

%decide whether single-cell requires repressor model or not for data set 1
ind1_2 = find(BIC1_2-BIC1_1<-10); %model 2 best
ind1_1 = find(BIC1_2-BIC1_1>=-10);%model 1 best

%extract BIC values for all single-cell trajectories for data set 2
for i = 1:size(scR2_1,2)
    BIC2_1(i) = scR2_1(i).sol.BIC;
end
for i = 1:size(scR2_2,2)
    BIC2_2(i) = scR2_2(i).sol.BIC;
end

%decide whether single-cell requires repressor model or not for data set 2
ind2_2 = find(BIC2_2-BIC2_1<-10);
ind2_1 = find(BIC2_2-BIC2_1>=-10);

%get the data sets for both experiments
load('NonDividing2')

%extract relevant information from non-dividing cell structure
if rep1 == 1
    data1 = NonDividing{strain1}.I5r1;
    momID1 = NonDividing{strain1}.I5momIDr1;
    mompos1 = NonDividing{strain1}.I5momposr1;
else
    data1 = NonDividing{strain1}.I5r2;
    momID1 = NonDividing{strain1}.I5momIDr2;
    mompos1 = NonDividing{strain1}.I5momposr2;
end

if rep2 == 1
    data2 = NonDividing{strain2}.I5r1;
    momID2 = NonDividing{strain2}.I5momIDr1;
    mompos2 = NonDividing{strain2}.I5momposr1;
else
    data2 = NonDividing{strain2}.I5r2;
    momID2 = NonDividing{strain2}.I5momIDr2;
    mompos2 = NonDividing{strain2}.I5momposr2;
end

%extract the estimated parameter sets per cell and data set
for icell = 1:length(scR1_2)
    clear par
    par = 10.^(scR1_2(icell).sol.MS.par(:,1));
    Par1(icell,:) = par';
end

for icell = 1:length(scR2_2)
    clear par
    par = 10.^(scR2_2(icell).sol.MS.par(:,1));
    Par2(icell,:) = par';
end

%reduce parameter sets to the cells requiring the repressor model
Par1 = Par1(ind1_2,:);
Par2 = Par2(ind2_2,:);

%determine color according ot strain and repression
if strain1 == 1
    if rep1 == 1
        c1 = [175,198,233]./255;
    else
        c1 = [33,68,120]./255;
    end
else
    if rep1 == 1
        c1 = [205,135,222]./255;
    else
        c1 = [67,31,117]./255;
    end
end

if strain2 == 1
    if rep2 == 1
        c2 = [175,198,233]./255;
    else
        c2 = [33,68,120]./255;
    end
else
    if rep2 == 1
        c2 = [205,135,222]./255;
    else
        c2 = [67,31,117]./255;
    end
end

% figure('visible','off');
figure

%plot GFP0 vs tdelay for both data sets
plot(Par2(:,1),Par2(:,2),'.','Color',c2,'Markersize',10)
hold on
plot(Par1(:,1),Par1(:,2),'.','Color',c1,'Markersize',10)
hold on

%get linear regression fit
con1 = polyfit(Par1(:,1),Par1(:,2),1);
% Display evaluated equation y = m*x + b
disp(['Repression 1: Equation is y = ' num2str(con1(1)) '*x + ' num2str(con1(2))]);
% Evaluate fit equation using polyval
y1_est = polyval(con1,Par1(:,1));
% Add trend line to plot
hold on
plot(Par1(:,1),y1_est,'-','Color',c1,'LineWidth',2)
hold on

con2 = polyfit(Par2(:,1),Par2(:,2),1);
% Display evaluated equation y = m*x + b
disp(['Repression 2: Equation is y = ' num2str(con2(1)) '*x + ' num2str(con2(2))]);
% Evaluate fit equation using polyval
y2_est = polyval(con2,Par2(:,1));
% Add trend line to plot
hold on
plot(Par2(:,1),y2_est,'-','Color',c2,'LineWidth',2)

ylim([0 inf])
box off
set(gca,'linewidth',1.02)
set(gca,'FontSize',11)
set(gca,'FontName','Arial')
xlabel('GFP0')
ylabel('tdelay')
if strain1 == 1
    xlim([0,8])
    xticks([0,4,8])
    yticks([0,1])
end

%save figure
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5 5])
print('-dpdf','./Figures/FigS2G','-painters')