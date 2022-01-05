function getNonDividing(iexp)

clearvars -except iexp
clc;

n_total = 0;

for istrain = 1:2
    
    clearvars -except istrain NonDividing fig_count iexp n_total
    
    %get correct data sets (Wt vs elp6)
    if istrain == 1
        clear Data
        
        if iexp == 1
            DataExp{1} = sprintf('ExpI5_pos1_Y208_WT');
            DataExp{2} = sprintf('ExpI5_pos2_Y208_WT');
            DataExp{3} = sprintf('ExpI5_pos4_Y208_WT');
            DataExp{4} = sprintf('ExpI5_pos5_Y208_WT');
            DataExp{5} = sprintf('ExpI5_pos8_Y208_WT');
            DataExp{6} = sprintf('ExpI5_pos9_Y208_WT');
            DataExp{7} = sprintf('ExpI5_pos10_Y208_WT');
            DataExp{8} = sprintf('ExpI5_pos11_Y208_WT');
            DataExp{9} = sprintf('ExpI5_pos12_Y208_WT');
            DataExp{10} = sprintf('ExpI5_pos13_Y208_WT');
            DataExp{11} = sprintf('ExpI5_pos14_Y208_WT');
            DataExp{12} = sprintf('ExpI5_pos15_Y208_WT');
            DataExp{13} = sprintf('ExpI5_pos16_Y208_WT');
        else
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
        end
        
    else
        
        clear Data
        if iexp == 1
            DataExp{1} = sprintf('ExpI5_pos19_Y1474_Elp6');
            DataExp{2} = sprintf('ExpI5_pos20_Y1474_Elp6');
            DataExp{3} = sprintf('ExpI5_pos22_Y1474_Elp6');
            DataExp{4} = sprintf('ExpI5_pos23_Y1474_Elp6');
            DataExp{5} = sprintf('ExpI5_pos24_Y1474_Elp6');
            DataExp{6} = sprintf('ExpI5_pos25_Y1474_Elp6');
            DataExp{7} = sprintf('ExpI5_pos27_Y1474_Elp6');
            DataExp{8} = sprintf('ExpI5_pos30_Y1474_Elp6');
            DataExp{9} = sprintf('ExpI5_pos31_Y1474_Elp6');
            DataExp{10} = sprintf('ExpI5_pos32_Y1474_Elp6');
        else
            DataExp{1} = sprintf('ExpL1_pos9_Y1474_Elp6');
            DataExp{2} = sprintf('ExpL1_pos10_Y1474_Elp6');
            DataExp{3} = sprintf('ExpL1_pos11_Y1474_Elp6');
            DataExp{4} = sprintf('ExpL1_pos12_Y1474_Elp6');
            DataExp{5} = sprintf('ExpL1_pos13_Y1474_Elp6');
            DataExp{6} = sprintf('ExpL1_pos14_Y1474_Elp6');
            DataExp{7} = sprintf('ExpL1_pos15_Y1474_Elp6');
            DataExp{8} = sprintf('ExpL1_pos16_Y1474_Elp6');
            DataExp{9} = sprintf('ExpL2_pos14_Y1474_Elp6');
            DataExp{10} = sprintf('ExpL2_pos15_Y1474_Elp6');
            DataExp{11} = sprintf('ExpL2_pos16_Y1474_Elp6');
            DataExp{12} = sprintf('ExpL2_pos17_Y1474_Elp6');
            DataExp{13} = sprintf('ExpL2_pos18_Y1474_Elp6');
            DataExp{14} = sprintf('ExpL2_pos19_Y1474_Elp6');
            DataExp{15} = sprintf('ExpL2_pos20_Y1474_Elp6');
            DataExp{16} = sprintf('ExpL2_pos21_Y1474_Elp6');
            DataExp{17} = sprintf('ExpL2_pos22_Y1474_Elp6');
            DataExp{18} = sprintf('ExpL2_pos23_Y1474_Elp6');
            DataExp{19} = sprintf('ExpL2_pos24_Y1474_Elp6');
            DataExp{20} = sprintf('ExpL2_pos25_Y1474_Elp6');
            DataExp{21} = sprintf('ExpL2_pos26_Y1474_Elp6');
            DataExp{22} = sprintf('ExpL2_pos27_Y1474_Elp6');
            DataExp{23} = sprintf('ExpL2_pos28_Y1474_Elp6');
        end
    end
    
    Data = DataExp;
    
    Data = Data(~cellfun('isempty',Data));
    
    %initialization of structures to save information on non-dividing cells
    momGFPr1 = zeros(500,40); %total GFP repression 1 for cells of specific position
    momGFPr1all = zeros(1,40); %total GFP repression 1 for all non-dividing cells of all positions
    momGFPr2 = zeros(500,40); %total GFP repression 2 for cells of specific position
    momGFPr2all = zeros(1,40); %total GFP repression 2 for all non-dividing cells of all positions
    mom_IDr1all = []; %cell IDs of mother cells for repression 1 of all positions
    mom_IDr2all = []; %cell IDs of mother cells for repression 2 of all positions
    mom_posr1all = []; %cell positions of mother cells for repression 1 of all positions
    mom_posr2all= []; %cell positions of mother cells for repression 2 of all positions
    daughters_IDr1all = []; %cell IDs of daughter cells for repressoin 1 of all positions
    mother_count_r1all = []; %number of daughters per mother cell in repression 1 for all positions
    daughters_IDr2all = []; %cell IDs of daughter cells for repression 2 of all positions
    mother_count_r2all = []; %number of daughters per mother cell in repression 2 for all positions
    
    n_cells_r1 = 0;
    n_cells_r2 = 0;
    
    %for every position do
    for i = 1:length(Data)
        
        clearvars -except Data istrain i momGFPr1all momGFPr2all...
            NonDividing n_cells_r1 n_cells_r2 mom_IDr1all mom_IDr2all...
            mom_posr1all mom_posr2all daughters_IDr1all mother_count_r1all...
            daughters_IDr2all mother_count_r2all iexp n_total
        
        loadData = sprintf('S%s',Data{i});
        load(loadData);
        
        count_mother_r1 = 1;
        count_mother_r2 = 1;
        count_daughter_r1 = 1;
        count_daughter_r2 = 1;
        
        daughters_IDr1 = [];
        daughters_IDr2 = [];
        
        clear cell_info
        
        %count all cells in WT and elp6
        n_total = n_total+length(S);
        
        %create summary matrix of all cells in that position
        for iS = 1:length(S)
            cell_info(iS,:) = [S{iS}.N, S{iS}.M, S{iS}.DF];
        end
        
        %sort according to the detection frame such that we look at the
        %cells in a timely ordered fashion
        [val,indval] = sort(cell_info(:,3));
        cell_info = cell_info(indval,:);
        
        %also reorder the cells
        for j = 1:length(S)
            Snew{j} = S{indval(j)};
        end
        S = Snew;
        
        %for every cell of position do
        for iS = 1:length(S)
            
            %only take cells into consideration if DF < 180  (= cell
            %present in first two hours of repression 1 = present between frames 141 and 180)
            if S{iS}.DF <= 180
                
                %if cell born before repression 1 (= frame 141) = mother cell
                if S{iS}.DF < 141
                    
                    %update cell count of mother cells in repression 1
                    n_cells_r1 = n_cells_r1+1;
                    
                    %find the frame number of cell for frame 140 (last
                    %frame before start of repression 1)
                    f141 = find(S{iS}.DF:S{iS}.LF == 140);
                    
                    %save total GFP of mother cell for frames 141:180
                    %(first two hours of repression 1)
                    momGFPr1(count_mother_r1,:) = S{iS}.GFPabs(f141+1:f141+40);
                    
                    %save mother ID and position and initialize daughter
                    %count to 0 for this mother cell
                    mom_IDr1(count_mother_r1) = cell_info(iS,1);
                    mom_posr1(count_mother_r1) = i;
                    mother_count_r1(count_mother_r1) = 0;
                    
                    %update counted number of mother cells for repression 1
                    count_mother_r1 = count_mother_r1+1;
                    
                else %if not a mother cell at start of repression 1 (cell detected after frame 140)
                    
                    %find the cell ID of mother cell of that cell
                    imother = find(cell_info(iS,2)==cell_info(:,1));
                    
                    %find the original predecessor cell of that daughter cell present at start of
                    %repression 1
                    if isempty(imother) == 0
                        
                        %while the predecessor was not present at frame
                        %141, continue looking for the mother of that
                        %predecessor cell until predecessor cell at start
                        %of repression 1 was found
                        while cell_info(imother,3)>=141
                            imother = find(cell_info(imother,2)==cell_info(:,1));
                        end
                        
                        %if this predecessor can be found do
                        if isempty(imother) == 0
                            
                            %identify the corresponding mother cell in
                            %momGFPr1 to add up daughter total GFP to
                            %correct mother total GFP
                            ind = length(find(cell_info(1:imother,3)<141));
                            
                            %add daughter GFP to mother GFP for the
                            %overlapping time frames between mother and
                            %daughter for the first two hours of repression
                            %1
                            momGFPr1(ind,S{iS}.DF-140:180-140) = ...
                                momGFPr1(ind,S{iS}.DF-140:180-140)+S{iS}.GFPabs(1:180-S{iS}.DF+1);
                            
                            %save daughter ID and update the count of
                            %daughter cells for this particular mother cell
                            daughters_IDr1(count_daughter_r1) = cell_info(iS,1);
                            mother_count_r1(ind) = mother_count_r1(ind)+1;
                            
                            %update counted number of daughter cells
                            count_daughter_r1 = count_daughter_r1+1;
                        end
                    end
                end
            end
            
            %if cell born before repression 2 (= frame 281) = mother cell
            if S{iS}.DF < 281
                
                %update cell count of mother cells in repression 1
                n_cells_r2 = n_cells_r2+1;
                
                %find the frame number of cell for frame 280 (last
                %frame before start of repression 2)
                f281 = find(S{iS}.DF:S{iS}.LF == 280);
                
                %save total GFP of mother cell for frames 281:last frame
                %(first two hours of repression 2)
                momGFPr2(count_mother_r2,:) = S{iS}.GFPabs(f281+1:f281+40);
                
                %save mother ID and position and initialize daughter
                %count to 0 for this mother cell
                mom_IDr2(count_mother_r2) = cell_info(iS,1);
                mom_posr2(count_mother_r2) = i;
                mother_count_r2(count_mother_r2) = 0;
                
                %update counted number of mother cells for repression 2
                count_mother_r2 = count_mother_r2+1;
                
            else %if not a mother cell at start of repression 2 (cell detected after frame 280)
                
                %find the cell ID of mother cell of that cell
                imother = find(cell_info(iS,2)==cell_info(:,1));
                
                %find the original predecessor cell of that daughter cell present at start of
                %repression 2
                if isempty(imother) == 0
                    
                    %while the predecessor was not present at frame
                    %281, continue looking for the mother of that
                    %predecessor cell until predecessor cell at start
                    %of repression 2 was found
                    while cell_info(imother,3)>=281
                        imother = find(cell_info(imother,2)==cell_info(:,1));
                    end
                    
                    %if this predecessor can be found do
                    if isempty(imother) == 0
                        
                        %identify the corresponding mother cell in
                        %momGFPr2 to add up daughter total GFP to
                        %correct mother total GFP
                        ind = length(find(cell_info(1:imother,3)<281));
                        
                        %add daughter GFP to mother GFP for the
                        %overlapping time frames between mother and
                        %daughter for the first two hours of repression 2
                        momGFPr2(ind,S{iS}.DF-280:320-280) = ...
                            momGFPr2(ind,S{iS}.DF-280:320-280)+S{iS}.GFPabs(1:320-S{iS}.DF+1);
                        
                        %save daughter ID and update the count of
                        %daughter cells for this particular mother cell
                        daughters_IDr2(count_daughter_r2) = cell_info(iS,1);
                        mother_count_r2(ind) = mother_count_r2(ind)+1;
                        
                        %update counted number of daughter cells
                        count_daughter_r2 = count_daughter_r2+1;
                        
                    end
                end
            end
            
        end
        
        %save non-dividng cell information of that position with information from other
        %positions (for parameter definitions - see above)
        momGFPr1all = [momGFPr1all; momGFPr1];
        momGFPr2all = [momGFPr2all; momGFPr2];
        mom_IDr1all = [mom_IDr1all, mom_IDr1];
        mom_IDr2all = [mom_IDr2all, mom_IDr2];
        mom_posr1all = [mom_posr1all, mom_posr1];
        mom_posr2all = [mom_posr2all, mom_posr2];
        daughters_IDr1all = [daughters_IDr1all,daughters_IDr1];
        mother_count_r1all = [mother_count_r1all,mother_count_r1];
        daughters_IDr2all = [daughters_IDr2all,daughters_IDr2];
        mother_count_r2all = [mother_count_r2all,mother_count_r2];
    end
    
    %remove all empty MATLAB cells of structure containing total GFP of non-dividing cell
    %(specified too many MATLAB cells at start)
    momGFPr1all(~any(momGFPr1all,2),:) = [];
    momGFPr2all(~any(momGFPr2all,2),:) = [];
    
    %save everything in MATLAB structure (for experiment I5)
    NonDividing{istrain}.r1 = momGFPr1all;
    NonDividing{istrain}.r2 = momGFPr2all;
    NonDividing{istrain}.nr1 =  n_cells_r1;
    NonDividing{istrain}.nr2 =  n_cells_r2;
    NonDividing{istrain}.momIDr1 =  mom_IDr1all;
    NonDividing{istrain}.momIDr2 =  mom_IDr2all;
    NonDividing{istrain}.momposr1 =  mom_posr1all;
    NonDividing{istrain}.momposr2 =  mom_posr2all;
    NonDividing{istrain}.momcountr1 = mother_count_r1all;
    NonDividing{istrain}.momcountr2 = mother_count_r2all;
    NonDividing{istrain}.daughtersIDr1 = daughters_IDr1all;
    NonDividing{istrain}.daughtersIDr2 = daughters_IDr2all;
    
    n_total
end

if iexp == 1
    save('./Melanie/Data/NonDividing1','NonDividing')
    %     xlswrite('./Data/NonDividing1_WT_rep1.xlsx',NonDividing{1}.I5r1,1,'A1');
    %     xlswrite('./Data/NonDividing1_WT_rep2.xlsx',NonDividing{1}.I5r2,1,'A1');
    %     xlswrite('./Data/NonDividing1_elp6_rep1.xlsx',NonDividing{2}.I5r1,1,'A1');
    %     xlswrite('./Data/NonDividing1_elp6_rep2.xlsx',NonDividing{2}.I5r2,1,'A1');
else
    save('./Melanie/Data/NonDividing2','NonDividing')
    %     xlswrite('./Data/NonDividing2_WT_rep1.xlsx',NonDividing{1}.I5r1,1,'A1');
    %     xlswrite('./Data/NonDividing2_WT_rep2.xlsx',NonDividing{1}.I5r2,1,'A1');
    %     xlswrite('./Data/NonDividing2_elp6_rep1.xlsx',NonDividing{2}.I5r1,1,'A1');
    %     xlswrite('./Data/NonDividing2_elp6_rep2.xlsx',NonDividing{2}.I5r2,1,'A1');
end

