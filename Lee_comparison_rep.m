clearvars;
clc;

figure
clearvars

i_allcell_r1 = 1;
i_allcell_r2 = 1;

%WT
DataExp{1} = sprintf('ExpL1_pos1_Y208_WT');
DataExp{2} = sprintf('ExpL1_pos5_Y208_WT');
DataExp{3} = sprintf('ExpL1_pos6_Y208_WT');
DataExp{4} = sprintf('ExpL1_pos7_Y208_WT');
DataExp{5} = sprintf('ExpL1_pos8_Y208_WT');
DataExp{6} = sprintf('ExpL2_pos1_Y208_WT');
DataExp{7} = sprintf('ExpL2_pos4_Y208_WT');
DataExp{8} = sprintf('ExpL2_pos6_Y208_WT');
DataExp{9} = sprintf('ExpL2_pos9_Y208_WT');
DataExp{10} = sprintf('ExpL2_pos10_Y208_WT');
DataExp{11} = sprintf('ExpL2_pos11_Y208_WT');
DataExp{12} = sprintf('ExpL2_pos12_Y208_WT');
DataExp{13} = sprintf('ExpI7_pos4_Y208_WT');

Data = DataExp;

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
print('-dpdf','./Figures/FigS1B_new','-painters')