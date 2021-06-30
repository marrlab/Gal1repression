clearvars;
clc;

figure('visible','off');

for istrain = 1:2
    
    clearvars -except iexp istrain divtime_r1all divtime_r2all
    
    %load data
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
    
    divtime_r1all{istrain} = [];
    divtime_r2all{istrain} = [];
    
    %go over all positions
    for i = 1:length(Data)
        
        clearvars -except Data i divtime_r1all divtime_r2all istrain iexp
        
        %load data for position
        loadData = sprintf('/Users/lea.schuh/Documents/GitHub/Gal1repression/Repression/S%s',Data{i});
        load(loadData);
        
        %reorder cells according to detection frame and get summary
        %matrix
        clear cell_info
        for iS = 1:length(S)
            cell_info(iS,:) = [S{iS}.N, S{iS}.M, S{iS}.DF];
        end
        [val,indval] = sort(cell_info(:,3));
        cell_info = cell_info(indval,:);
        
        %reorder cells
        for j = 1:length(S)
            Snew{j} = S{indval(j)};
        end
        S = Snew;
        
        %for every cell do
        for iS = 1:length(S)
            
            %only take cells into consideration if DF < 180 (cell in
            %first two hours of repression 1)
            if S{iS}.DF <= 180
                
                %if cell detected before repression 1
                if S{iS}.DF < 141
                    
                    %cell gets division time of 0
                    divtime_r1(iS) = 0;
                    
                    %if cell detected during first two hours of repression 1
                else
                    
                    %get the mother cell ID
                    ind = find(cell_info(:,1)==cell_info(iS,2));
                    
                    %if the mother cell is recorded do
                    if isempty(ind) == 0
                        
                        %if mother cell detected before repression 1 do
                        if cell_info(ind,3) < 141
                            
                            %division time of this cell is given by the
                            %the detection frame to the cell and the
                            %start frame of repression 1 (as during
                            %induction yeast cells do not proliferate)
                            divtime_r1(iS) = cell_info(iS,3)-140;
                            
                            %if mother cell is detected during repression 1
                            %do
                        else
                            
                            %division time form cell given by detection
                            %frame daughter - detection frame mother
                            divtime_r1(iS) = cell_info(iS,3)-cell_info(ind,3);
                        end
                    else
                        
                        %do not consider cell if mother cell is not
                        %recorded
                        divtime_r1(iS) = 0;
                    end
                end
            end
            
            %if cell born before repression 2
            if S{iS}.DF < 281
                
                %set vision time to 0 (do not count)
                divtime_r2(iS) = 0;
                
                %if cell born during repression 2
            else
                
                %get the mother cell ID
                ind = find(cell_info(:,1)==cell_info(iS,2));
                
                %if the mother cell is recorded do
                if isempty(ind) == 0
                    
                    %if mother cell detected before repression 2 do
                    if cell_info(ind,3) < 281
                        
                        %division time of this cell is given by the
                        %the detection frame to the cell and the
                        %start frame of repression 2 (as during
                        %induction yeast cells do not proliferate)
                        divtime_r2(iS) = cell_info(iS,3)-280;
                        
                        %if mother cell is detected during repression 2
                        %do
                    else
                        
                        %division time form cell given by detection
                        %frame daughter - detection frame mother
                        divtime_r2(iS) = cell_info(iS,3)-cell_info(ind,3);
                        
                    end
                else
                    
                    %do not consider cell if mother cell is not
                    %recorded
                    divtime_r2(iS) = 0;
                end
            end
            
        end
        
        %record all positions together
        divtime_r1all{istrain} = [divtime_r1all{istrain},divtime_r1];
        divtime_r2all{istrain} = [divtime_r2all{istrain},divtime_r2];
    end
end

%     boxplot(Pall,index,'positions',index,'labels',labels(unique(index)));
divtime_r1all{1} = divtime_r1all{1}(divtime_r1all{1}~=0);
divtime_r1all{2} = divtime_r1all{2}(divtime_r1all{2}~=0);
divtime_r2all{1} = divtime_r2all{1}(divtime_r2all{1}~=0);
divtime_r2all{2} = divtime_r2all{2}(divtime_r2all{2}~=0);

%plot with jitter
a = -0.2;
b = 0.2;
r = (b-a).*rand(length(divtime_r1all{1}),1) + a;
plot(1+r,divtime_r1all{1}*3/60,'.','Color',[175,198,233]./255,'Markersize',10)
hold on
line([0.6,1.4],[median(divtime_r1all{1}*3/60),median(divtime_r1all{1}*3/60)],'Color','k','Linewidth',2)

a = -0.2;
b = 0.2;
r = (b-a).*rand(length(divtime_r1all{2}),1) + a;
plot(2+r,divtime_r1all{2}*3/60,'.','Color',[205,135,222]./255,'Markersize',10)
hold on
line([1.6,2.4],[median(divtime_r1all{2}*3/60),median(divtime_r1all{2}*3/60)],'Color','k','Linewidth',2)

a = -0.2;
b = 0.2;
r = (b-a).*rand(length(divtime_r2all{1}),1) + a;
plot(3+r,divtime_r2all{1}*3/60,'.','Color',[33,68,120]./255,'Markersize',10)
hold on
line([2.6,3.4],[median(divtime_r2all{1}*3/60),median(divtime_r2all{1}*3/60)],'Color','k','Linewidth',2)

a = -0.2;
b = 0.2;
r = (b-a).*rand(length(divtime_r2all{2}),1) + a;
plot(4+r,divtime_r2all{2}*3/60,'.','Color',[67,31,117]./255,'Markersize',10)
hold on
line([3.6,4.4],[median(divtime_r2all{2}*3/60),median(divtime_r2all{2}*3/60)],'Color','k','Linewidth',2)

ylim([0 2])
box off
set(gca,'linewidth',1.02)
set(gca,'FontSize',11)
set(gca,'FontName','Arial')
xlim([0,4.5])

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5 5])
print('-dpdf','./Figures/SuppFig1')

kruskalwallis([divtime_r1all{1}*3/60,divtime_r1all{2}*3/60,divtime_r2all{1}*3/60,divtime_r2all{2}*3/60],...
    [ones(1,length(divtime_r1all{1})),2*ones(1,length(divtime_r1all{2})),3*ones(1,length(divtime_r2all{1})),4*ones(1,length(divtime_r2all{2}))])

