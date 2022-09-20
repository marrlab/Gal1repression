%add path
addpath(genpath(pwd))

%% Figure 1B - full data set WT
clearvars;
clc;

figure('visible','off');
clearvars

%WT
DataExpI5{1} = sprintf('ExpI5_pos1_Y208_WT');
DataExpI5{2} = sprintf('ExpI5_pos2_Y208_WT');
DataExpI5{3} = sprintf('ExpI5_pos4_Y208_WT');
DataExpI5{4} = sprintf('ExpI5_pos5_Y208_WT');
DataExpI5{5} = sprintf('ExpI5_pos8_Y208_WT');
DataExpI5{6} = sprintf('ExpI5_pos9_Y208_WT');
DataExpI5{7} = sprintf('ExpI5_pos10_Y208_WT');
DataExpI5{8} = sprintf('ExpI5_pos11_Y208_WT');
DataExpI5{9} = sprintf('ExpI5_pos12_Y208_WT');
DataExpI5{10} = sprintf('ExpI5_pos13_Y208_WT');
DataExpI5{11} = sprintf('ExpI5_pos14_Y208_WT');
DataExpI5{12} = sprintf('ExpI5_pos15_Y208_WT');
DataExpI5{13} = sprintf('ExpI5_pos16_Y208_WT');

Data = DataExpI5;

Data = Data(~cellfun('isempty',Data));
count = 0;

for i = 1:length(Data)
    
    clearvars -except Data i istrain
    
    %load different positions
    loadData = sprintf('S%s',Data{i});
    load(loadData);
    
    %for every cell in that position do
    for iS = 1:length(S)
        
        %color definition
        c = [175,198,233;160,160,160;117,157,233;160,160,160;33,68,120]./255;
        
        %repression 0
        ind81 = find(S{iS}.DF:S{iS}.LF==81);
        if isempty(ind81) == 0
            plot((81-length(S{iS}.GFPabs(1:ind81)))*3/60:3/60:(81-1)*3/60,...
                S{iS}.GFPabs(1:ind81)./10e6,'-','Color',c(1,:))
            hold on
        else
            ind81 = 1;
        end
        
        %induction 1
        ind141 = find(S{iS}.DF:S{iS}.LF==141);
        if isempty(ind141) == 0
            plot((141-length(S{iS}.GFPabs(ind81:ind141)))*3/60:3/60:(141-1)*3/60,...
                S{iS}.GFPabs(ind81:ind141)./10e6,'-','Color',c(2,:))
            hold on
        else
            ind141 = 1;
        end
        
        %repression 1
        ind221 = find(S{iS}.DF:S{iS}.LF==221);
        if isempty(ind221) == 0
            plot((221-length(S{iS}.GFPabs(ind141:ind221)))*3/60:3/60:(221-1)*3/60,...
                S{iS}.GFPabs(ind141:ind221)./10e6,'-','Color',c(3,:))
            hold on
        else
            ind221 = 1;
        end
        
        %induction 2
        ind281 = find(S{iS}.DF:S{iS}.LF==281);
        if isempty(ind281) == 0
            plot((281-length(S{iS}.GFPabs(ind221:ind281)))*3/60:3/60:(281-1)*3/60,...
                S{iS}.GFPabs(ind221:ind281)./10e6,'-','Color',c(4,:))
            hold on
        else
            ind281 = 1;
        end
        
        %repression 2
        plot((320-length(S{iS}.GFPabs(ind281:end)))*3/60:3/60:(320-1)*3/60,...
            S{iS}.GFPabs(ind281:end)./10e6,'-','Color',c(5,:))
        hold on
    end
    
end

ylabel('total Gal1-GFP fluorescence (a.u.)')
yticks([0,5,10,15])
xlabel('experimental time (h)')
xticks([0,4,7,11,14,16])
box off
set(gca,'linewidth',1.02)
set(gca,'FontSize',11)
set(gca,'FontName','Arial')
ylim([0,15])
xlim([0,16])
set(gcf, 'DefaultFigureRenderer', 'painters');
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.5 5])
print('-dpdf','./Figures/Fig1B','-painters')

%% Figure 1C - normalized total GFP for repressions r1 and r2

clearvars;
clc;

figure('visible','off');
clearvars

i_allcell_r1 = 1;
i_allcell_r2 = 1;

%WT
DataExpI5{1} = sprintf('ExpI5_pos1_Y208_WT');
DataExpI5{2} = sprintf('ExpI5_pos2_Y208_WT');
DataExpI5{3} = sprintf('ExpI5_pos4_Y208_WT');
DataExpI5{4} = sprintf('ExpI5_pos5_Y208_WT');
DataExpI5{5} = sprintf('ExpI5_pos8_Y208_WT');
DataExpI5{6} = sprintf('ExpI5_pos9_Y208_WT');
DataExpI5{7} = sprintf('ExpI5_pos10_Y208_WT');
DataExpI5{8} = sprintf('ExpI5_pos11_Y208_WT');
DataExpI5{9} = sprintf('ExpI5_pos12_Y208_WT');
DataExpI5{10} = sprintf('ExpI5_pos13_Y208_WT');
DataExpI5{11} = sprintf('ExpI5_pos14_Y208_WT');
DataExpI5{12} = sprintf('ExpI5_pos15_Y208_WT');
DataExpI5{13} = sprintf('ExpI5_pos16_Y208_WT');

Data = DataExpI5;

Data = Data(~cellfun('isempty',Data));
count = 0;

for i = 1:length(Data)
    
    clearvars -except Data i istrain i_allcell_r1 i_allcell_r2 GFPmean_r1 GFPmean_r2
    
    %load different positions
    loadData = sprintf('S%s',Data{i});
    load(loadData);
    
    %for every cell in that position do
    for iS = 1:length(S)
        
        %color definition
        c = [175,198,233;160,160,160;117,157,233;160,160,160;33,68,120]./255;
        
        %repression 0
        ind81 = find(S{iS}.DF:S{iS}.LF==81);
        if isempty(ind81) == 0
        else
            ind81 = 1;
        end
        
        %induction 1
        ind141 = find(S{iS}.DF:S{iS}.LF==141);
        if isempty(ind141) == 0
        else
            ind141 = 1;
        end
        
        %repression 1
        ind221 = find(S{iS}.DF:S{iS}.LF==221);
        if isempty(ind221) == 0
            GFPmean_r1(i_allcell_r1,221-length(S{iS}.GFPabs(ind141:ind221))+1-140:221-140) = S{iS}.GFPabs(ind141:ind221)./10e6;
            i_allcell_r1 = i_allcell_r1+1;
            hold on
        else
            ind221 = 1;
        end
        
        %induction 2
        ind281 = find(S{iS}.DF:S{iS}.LF==281);
        if isempty(ind281) == 0
            hold on
        else
            ind281 = 1;
        end
        
        %repression 2
        GFPmean_r2(i_allcell_r2,320-length(S{iS}.GFPabs(ind281:end))+1-280:320-280) = S{iS}.GFPabs(ind281:end)./10e6;
            i_allcell_r2 = i_allcell_r2+1;
        hold on
    end
    
end

plot(0*3/60:3/60:39*3/60,mean(GFPmean_r1(:,1:40))/mean(GFPmean_r1(:,1)),'-','Color',c(3,:),'Linewidth',1)
hold on
plot(0*3/60:3/60:(319-280)*3/60,mean(GFPmean_r2)/mean(GFPmean_r2(:,1)),'-','Color',c(5,:),'Linewidth',1)

ylabel('normalized Gal1-GFP fluorescence (a.u.)')
yticks([0,0.5,1,1.5])
xlabel('experimental time (h)')
xticks([0,1,2])
box off
set(gca,'linewidth',1.02)
set(gca,'FontSize',11)
set(gca,'FontName','Arial')
ylim([0,2])
xlim([0,2])
set(gcf, 'DefaultFigureRenderer', 'painters');
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.5 5])
print('-dpdf','./Figures/Fig1C','-painters')
