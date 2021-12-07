for ipos = 27
    
    clearvars -except ipos
    clc
    
%     if ipos < 10
%         datafile = sprintf('/Volumes/GAL1MEMORY/20210914_1-16_WT_17-32_elp6/TIFFs/Position_%d/Images/exp001_pos%d_s0%d_acdc_output.csv',ipos,ipos+8,ipos);
%         T = readtable(datafile);
%     else
%         datafile = sprintf('/Volumes/GAL1MEMORY/20210914_1-16_WT_17-32_elp6/TIFFs/Position_%d/Images/exp001_pos%d_s%d_acdc_output.csv',ipos,ipos+8,ipos);
%         T = readtable(datafile);
%     end

    if ipos < 10
        datafile = sprintf('/Volumes/GAL1MEMORY/20211130_1-16_WT_17-32_elp6/TIFFs/Position_%d/Images/exp001_pos%d_s0%d_acdc_output.csv',ipos,ipos+3,ipos);
        T = readtable(datafile);
    else
        datafile = sprintf('/Volumes/GAL1MEMORY/20211130_1-16_WT_17-32_elp6/TIFFs/Position_%d/Images/exp001_pos%d_s%d_acdc_output.csv',ipos,ipos+4,ipos);
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
    
%     if ipos < 9
%        saveS = sprintf('./Data/ExpL1_pos%d_Y208_WT',ipos);
%         save(saveS,'segmentation') 
%     else
%         saveS = sprintf('./Data/ExpL1_pos%d_Y1474_Elp6',ipos);
%         save(saveS,'segmentation')
%     end
    
    if ipos < 14
       saveS = sprintf('./Data/ExpL2_pos%d_Y208_WT',ipos);
        save(saveS,'segmentation') 
    else
        saveS = sprintf('./Data/ExpL2_pos%d_Y1474_Elp6',ipos);
        save(saveS,'segmentation')
    end
end