addpath(genpath(pwd))

%add Phylocell path
%if no Phylocell installed - skip this step and continue using the
%SExpI5_posX_Y_Z files in the data folder
addpath(genpath('/Users/lea.schuh/Downloads/phyloCell-dev-gcharvin-2-2.2/phyloCell'))
clearvars;
clc;

%define the data set - define correctly segmented positions for WT and elp6
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

Data = DataExpI5;

Data = Data(~cellfun('isempty',Data));

for i = 1:length(Data)
    
    clearvars -except Data i
    clc;
    
    %load specific single cell data information
    loadData = sprintf('./Segmentation/%s',Data{i});
    load(loadData);
    
    close all force
    
    count = 1;
    %go thorugh each cell and extract relevant single cell information
    %required for further analysis
    for icell = 1:length(segmentation.tcells1)
        %reject MATLAb cells without any single cell information
        if segmentation.tcells1(icell).N ~= 0
            %extract and safe only the important information
            DF(count) = segmentation.tcells1(icell).detectionFrame;
            LF(count) = segmentation.tcells1(icell).lastFrame;
            BF(count) = segmentation.tcells1(icell).birthFrame;
            M(count) = segmentation.tcells1(icell).mother;
            N(count) = segmentation.tcells1(icell).N;
            
            for itime = 1:length(segmentation.tcells1(icell).Obj)
                %sort according to image as some are not sequential
                I{count}(itime) = segmentation.tcells1(icell).Obj(itime).image;
                usGFPrel{count}(itime) = segmentation.tcells1(icell).Obj(itime).fluoMean(2); %GFP channel
                usArea{count}(itime) = segmentation.tcells1(icell).Obj(itime).area;
            end
            
            [valI,indI] = sort(I{count},'ascend');
            GFPrel{count} = usGFPrel{count}(indI);
            %get total GFP per cell per time
            GFPabs{count} = usGFPrel{count}(indI) .* usArea{count}(indI);
            Area{count} = usArea{count}(indI);
            
            %collect and save information
            S1{count}.DF = DF(count);
            S1{count}.LF = LF(count);
            S1{count}.M = M(count);
            S1{count}.N = N(count);
            S1{count}.GFPabs = GFPabs{count};
            S1{count}.Area = Area{count};
            count = count + 1;
        end
    end
    
    %make summary matrix out of all cells
    for iS1 = 1:length(S1)
        SAll1(iS1,1) = S1{iS1}.DF;
        SAll1(iS1,2) = S1{iS1}.LF;
        SAll1(iS1,4) = S1{iS1}.M;
        SAll1(iS1,5) = S1{iS1}.N;
    end
    
    count = 1;
    for iS1 = 1:length(S1)
        %if cell not tracked til last frame - discard
        if S1{iS1}.LF == mode(SAll1(:,2))
            %assure that for every frame we have the total GFP value
            if length(S1{iS1}.DF:S1{iS1}.LF) == length(S1{iS1}.GFPabs)
                %discard cells with DF before DF of mother cell
                if SAll1(iS1,4) > 0 %if not initial cell
                    Mcell = find(SAll1(:,5)==SAll1(iS1,4));
                    if S1{iS1}.DF > SAll1(Mcell,1)
                        S2{count} = S1{iS1};
                        count = count + 1;
                    end
                else
                    S2{count} = S1{iS1};
                    count = count + 1;
                end
            end
        end
    end
    
    %rename to S
    S = S2;
    
    saveS = sprintf('./Data/S%s',Data{i});
    save(saveS,'S');
    
end
