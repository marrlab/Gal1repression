%add path
addpath(genpath(pwd))

%% Figure 2A - full data set WT
clearvars;
clc;

figure('visible','off');
for istrain = 1
    
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
        DataExpI5{15} = sprintf('ExpI5_pos19_Y1474_Elp6');
        DataExpI5{16} = sprintf('ExpI5_pos20_Y1474_Elp6');
        DataExpI5{17} = sprintf('ExpI5_pos22_Y1474_Elp6');
        DataExpI5{18} = sprintf('ExpI5_pos23_Y1474_Elp6');
        DataExpI5{19} = sprintf('ExpI5_pos24_Y1474_Elp6');
        DataExpI5{20} = sprintf('ExpI5_pos25_Y1474_Elp6');
        DataExpI5{21} = sprintf('ExpI5_pos27_Y1474_Elp6');
        DataExpI5{22} = sprintf('ExpI5_pos30_Y1474_Elp6');
        DataExpI5{23} = sprintf('ExpI5_pos31_Y1474_Elp6');
        DataExpI5{24} = sprintf('ExpI5_pos32_Y1474_Elp6');
        
    end
    
    Data = DataExpI5;
    
    Data = Data(~cellfun('isempty',Data));
    count = 0;
    
    for i = 1:length(Data)
        
        clearvars -except Data i istrain
        
        %load different positions 
        loadData = sprintf('./Data/S%s',Data{i});
        load(loadData);
        
        %for every cell in that position do
        for iS = 1:length(S)
            
            %color definition
            if istrain == 1
                c = [0,0,0;160,160,160;175,198,233;160,160,160;33,68,120]./255;
            else
                c = [0,0,0;160,160,160;205,135,222;160,160,160;67,31,117]./255;
            end
            
            %repression 0
            ind81 = find(S{iS}.DF:S{iS}.LF==81);
            if isempty(ind81) == 0
                plot((81-length(S{iS}.GFPabs(1:ind81)))*3/60:3/60:(81-1)*3/60,S{iS}.GFPabs(1:ind81),'-','Color',c(1,:))
                hold on
            else
                ind81 = 1;
            end
            
            %induction 1
            ind141 = find(S{iS}.DF:S{iS}.LF==141);
            if isempty(ind141) == 0
                plot((141-length(S{iS}.GFPabs(ind81:ind141)))*3/60:3/60:(141-1)*3/60,S{iS}.GFPabs(ind81:ind141),'-','Color',c(2,:))
                hold on
            else
                ind141 = 1;
            end
            
            %repression 1
            ind221 = find(S{iS}.DF:S{iS}.LF==221);
            if isempty(ind221) == 0
                plot((221-length(S{iS}.GFPabs(ind141:ind221)))*3/60:3/60:(221-1)*3/60,S{iS}.GFPabs(ind141:ind221),'-','Color',c(3,:))
                hold on
            else
                ind221 = 1;
            end
            
            %induction 2
            ind281 = find(S{iS}.DF:S{iS}.LF==281);
            if isempty(ind281) == 0
                plot((281-length(S{iS}.GFPabs(ind221:ind281)))*3/60:3/60:(281-1)*3/60,S{iS}.GFPabs(ind221:ind281),'-','Color',c(4,:))
                hold on
            else
                ind281 = 1;
            end
            
            %repression 2
            plot((S{iS}.LF-length(S{iS}.GFPabs(ind281:end)))*3/60:3/60:(S{iS}.LF-1)*3/60,S{iS}.GFPabs(ind281:end),'-','Color',c(5,:))
            hold on
        end
        
    end
end
ylabel('GFP intensity')
yticks([0,5e7,10e7,15e7])
xlabel('time (h)')
xticks([0,4,7,11,14,16])
box off
set(gca,'linewidth',1.02)
set(gca,'FontSize',11)
set(gca,'FontName','Arial')
ylim([0,15e7])
xlim([0,16])
set(gcf, 'DefaultFigureRenderer', 'painters');
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 9 5])
print('-dpdf','./Figures/Fig2A','-painters')

%% Figure 2B - full data set elp6
clearvars;
clc;

figure('visible','off');
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
        DataExpI5{15} = sprintf('ExpI5_pos19_Y1474_Elp6');
        DataExpI5{16} = sprintf('ExpI5_pos20_Y1474_Elp6');
        DataExpI5{17} = sprintf('ExpI5_pos22_Y1474_Elp6');
        DataExpI5{18} = sprintf('ExpI5_pos23_Y1474_Elp6');
        DataExpI5{19} = sprintf('ExpI5_pos24_Y1474_Elp6');
        DataExpI5{20} = sprintf('ExpI5_pos25_Y1474_Elp6');
        DataExpI5{21} = sprintf('ExpI5_pos27_Y1474_Elp6');
        DataExpI5{22} = sprintf('ExpI5_pos30_Y1474_Elp6');
        DataExpI5{23} = sprintf('ExpI5_pos31_Y1474_Elp6');
        DataExpI5{24} = sprintf('ExpI5_pos32_Y1474_Elp6');
        
    end
    
    Data = DataExpI5;
    
    Data = Data(~cellfun('isempty',Data));
    count = 0;
    
    for i = 1:length(Data)
        
        clearvars -except Data i istrain
        
        %load data of position
        loadData = sprintf('./Data/S%s',Data{i});
        load(loadData);
        
        %for every cell in position do
        for iS = 1:length(S)
            
            %color definition
            if istrain == 1
                c = [0,0,0;160,160,160;175,198,233;160,160,160;33,68,120]./255;
            else
                c = [0,0,0;160,160,160;205,135,222;160,160,160;67,31,117]./255;
            end
            
            %repression 0
            ind81 = find(S{iS}.DF:S{iS}.LF==81);
            if isempty(ind81) == 0
                plot((81-length(S{iS}.GFPabs(1:ind81)))*3/60:3/60:(81-1)*3/60,S{iS}.GFPabs(1:ind81),'-','Color',c(1,:))
                hold on
            else
                ind81 = 1;
            end
            
            %induction 1
            ind141 = find(S{iS}.DF:S{iS}.LF==141);
            if isempty(ind141) == 0
                plot((141-length(S{iS}.GFPabs(ind81:ind141)))*3/60:3/60:(141-1)*3/60,S{iS}.GFPabs(ind81:ind141),'-','Color',c(2,:))
                hold on
            else
                ind141 = 1;
            end
            
            %repression 1
            ind221 = find(S{iS}.DF:S{iS}.LF==221);
            if isempty(ind221) == 0
                plot((221-length(S{iS}.GFPabs(ind141:ind221)))*3/60:3/60:(221-1)*3/60,S{iS}.GFPabs(ind141:ind221),'-','Color',c(3,:))
                hold on
            else
                ind221 = 1;
            end
            
            %induction 2
            ind281 = find(S{iS}.DF:S{iS}.LF==281);
            if isempty(ind281) == 0
                plot((281-length(S{iS}.GFPabs(ind221:ind281)))*3/60:3/60:(281-1)*3/60,S{iS}.GFPabs(ind221:ind281),'-','Color',c(4,:))
                hold on
            else
                ind281 = 1;
            end
            
            %repression 2
            plot((S{iS}.LF-length(S{iS}.GFPabs(ind281:end)))*3/60:3/60:(S{iS}.LF-1)*3/60,S{iS}.GFPabs(ind281:end),'-','Color',c(5,:))
            hold on
        end
        
    end
end
ylabel('GFP intensity')
yticks([0,5e7,10e7,15e7])
xlabel('time (h)')
xticks([0,4,7,11,14,16])
box off
set(gca,'linewidth',1.02)
set(gca,'FontSize',11)
set(gca,'FontName','Arial')
ylim([0,15e7])
xlim([0,16])
set(gcf, 'DefaultFigureRenderer', 'painters');
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 9 5])
print('-dpdf','./Figures/Fig2B','-painters')

%% Figure 1D-E - non-dividing total GFP traces for repressions 1 and 2 for WT and elp6

clearvars;
clc;

load('NonDividing')
figure('visible','off');

%istrain = 1 - WT / = 2 - elp6
%irep = 1 - repression 1 / = 2 - repression 2

for istrain = 2
    for irep = 2
        
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
        
        %plot the non-dividing cells of specified strain and repression 
        %also plot mean total GFP trace and maximal mean total GFP value
        if irep == 1
            plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r1(:,1:40),'-','Color',c)
            hold on
            plot((1-1)*3/60:3/60:(40-1)*3/60,mean(NonDividing{istrain}.I5r1(:,1:40)),'-','Color','k')
            hold on
            indmax = find(mean(NonDividing{istrain}.I5r1(:,1:40))==max(mean(NonDividing{istrain}.I5r1(:,1:40))));
            time = (1-1)*3/60:3/60:(40-1)*3/60;
            time(indmax)
            plot(time(indmax),max(mean(NonDividing{istrain}.I5r1(:,1:40))),'.','Color','k','Markersize',12)
        else
            plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.I5r2(:,1:40),'-','Color',c)
            hold on
            plot((1-1)*3/60:3/60:(40-1)*3/60,mean(NonDividing{istrain}.I5r2(:,1:40)),'-','Color','k')
            hold on
            indmax = find(mean(NonDividing{istrain}.I5r2(:,1:40))==max(mean(NonDividing{istrain}.I5r2(:,1:40))));
            time = (1-1)*3/60:3/60:(40-1)*3/60;
            time(indmax)
            plot(time(indmax),max(mean(NonDividing{istrain}.I5r2(:,1:40))),'.','Color','k','Markersize',12)
        end
        
        ylabel('GFP intensity')
        if irep == 1
            yticks([0,2e7,4e7])
            ylim([0,5e7])
        else
            yticks([0,5e7,10e7,15e7])
            ylim([0,15e7])
        end
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
                print('-dpdf','./Figures/Fig2Dleft','-painters')
            else
                print('-dpdf','./Figures/Fig2Dright','-painters')
            end
        else
            if irep == 1
                print('-dpdf','./Figures/Fig2Eleft','-painters')
            else
                print('-dpdf','./Figures/Fig2Eright','-painters')
            end
        end
    end
end