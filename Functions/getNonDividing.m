clearvars;
clc;

for istrain = 1:2
    
    clearvars -except istrain NonDividing fig_count
    
    %get correct data sets (Wt vs elp6)
    if istrain == 1
        clear Data
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
            daughters_IDr2all mother_count_r2all
        
        loadData = sprintf('./Data/S%s',Data{i});
        load(loadData);
        
        count_mother_r1 = 1;
        count_mother_r2 = 1;
        count_daughter_r1 = 1;
        count_daughter_r2 = 1;
        
        clear cell_info
        
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
                momGFPr2(count_mother_r2,:) = S{iS}.GFPabs(f281+1:end);
                
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
                        momGFPr2(ind,S{iS}.DF-280:S{iS}.LF-280) = ...
                            momGFPr2(ind,S{iS}.DF-280:S{iS}.LF-280)+S{iS}.GFPabs;
                        
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
    NonDividing{istrain}.I5r1 = momGFPr1all;
    NonDividing{istrain}.I5r2 = momGFPr2all;
    NonDividing{istrain}.I5nr1 =  n_cells_r1;
    NonDividing{istrain}.I5nr2 =  n_cells_r2;
    NonDividing{istrain}.I5momIDr1 =  mom_IDr1all;
    NonDividing{istrain}.I5momIDr2 =  mom_IDr2all;
    NonDividing{istrain}.I5momposr1 =  mom_posr1all;
    NonDividing{istrain}.I5momposr2 =  mom_posr2all;
    NonDividing{istrain}.I5momcountr1 = mother_count_r1all;
    NonDividing{istrain}.I5momcountr2 = mother_count_r2all;
    NonDividing{istrain}.I5daughtersIDr1 = daughters_IDr1all;
    NonDividing{istrain}.I5daughtersIDr2 = daughters_IDr2all;
end

save('./Data/NonDividing','NonDividing')
