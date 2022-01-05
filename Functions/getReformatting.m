for iexp = 1:3
    
    clear P
    
    if iexp == 1
        P = [1,5:16]; %L1
    elseif iexp == 2
        P = [1,4,6,9:12,14:28]; %L2
    else
        P = 4; % I7
    end
    
    for ipos = 1:length(P)
        
        clearvars -except ipos iexp P
        clc
        
        if iexp == 1 %L1
            if P(ipos) < 10
                datafile = sprintf('/Volumes/GAL1MEMORY/20210914_1-16_WT_17-32_elp6/TIFFs/Position_%d/Images/exp001_pos%d_s0%d_acdc_output.csv',P(ipos),P(ipos)+8,P(ipos));
                T = readtable(datafile);
            else
                datafile = sprintf('/Volumes/GAL1MEMORY/20210914_1-16_WT_17-32_elp6/TIFFs/Position_%d/Images/exp001_pos%d_s%d_acdc_output.csv',P(ipos),P(ipos)+8,P(ipos));
                T = readtable(datafile);
            end
            
        elseif iexp == 2 %L2
            add = [0,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4];
            if P(ipos) < 10
                datafile = sprintf('/Volumes/GAL1MEMORY/20211130_1-16_WT_17-32_elp6/TIFFs/Position_%d/Images/exp001_pos%d_s0%d_acdc_output.csv',P(ipos),P(ipos)+add(ipos),P(ipos));
                T = readtable(datafile);
            else
                datafile = sprintf('/Volumes/GAL1MEMORY/20211130_1-16_WT_17-32_elp6/TIFFs/Position_%d/Images/exp001_pos%d_s%d_acdc_output.csv',P(ipos),P(ipos)+add(ipos),P(ipos));
                T = readtable(datafile);
            end
            
        else %I7
            datafile = sprintf('/Volumes/GAL1MEMORY/20191220/TIFFs/Position_%d/Images/WT_1-5_ELP2_6-10_ELP3_11-15_ELP4_16-20_ELP6_21-25_R22_26-30_pos%d_s0%d_acdc_output.csv',P(ipos),P(ipos),P(ipos));
            T = readtable(datafile);
        end
        
        count1 = 1;
        
        for i = 1:max(T.Cell_ID)
            clear rows
            rows = find(T.Cell_ID==i);
            if ~isempty(rows)
                segmentation.tcells1(count1).N = T.Cell_ID(rows(1));
                segmentation.tcells1(count1).detectionFrame = T.frame_i(rows(1))+1;
                segmentation.tcells1(count1).lastFrame = T.frame_i(rows(end))+1;
                segmentation.tcells1(count1).mother = T.relative_ID(rows(1));
                
                count2 = 1;
                for j = 1:length(rows)
                    
                    segmentation.tcells1(count1).Obj(count2).image = T.frame_i(rows(j))+1;
                    segmentation.tcells1(count1).Obj(count2).fluoMean(2) = T.GFP_mean(rows(j));
                    segmentation.tcells1(count1).Obj(count2).area = T.cell_area_pxl(rows(j));
                    count2 = count2+1;
                    
                end
                
                count1 = count1+1;
            end
        end
        
        if iexp == 1 %L1
            if P(ipos) < 9
                saveS = sprintf('./Melanie/Segmentation/ExpL1_pos%d_Y208_WT',P(ipos));
                save(saveS,'segmentation')
            else
                saveS = sprintf('./Melanie/Segmentation/ExpL1_pos%d_Y1474_Elp6',P(ipos));
                save(saveS,'segmentation')
            end   
        elseif iexp == 2 %L2
            if P(ipos) < 14
                saveS = sprintf('./Melanie/Segmentation/ExpL2_pos%d_Y208_WT',P(ipos));
                save(saveS,'segmentation')
            else
                saveS = sprintf('./Melanie/Segmentation/ExpL2_pos%d_Y1474_Elp6',P(ipos));
                save(saveS,'segmentation')
            end 
        else %I7
            saveS = sprintf('./Melanie/Segmentation/ExpI7_pos%d_Y208_WT',P(ipos));
            save(saveS,'segmentation')
        end
        
    end
end