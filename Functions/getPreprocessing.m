addpath(genpath(pwd))

%add Phylocell path
%if no Phylocell installed - skip this step and continue using the
%SExpI5_posX_Y_Z files in the data folder
addpath(genpath('/Users/lea.schuh/Downloads/phyloCell-dev-gcharvin-2-2.2/phyloCell'))
clearvars;
clc;

if ~exist('./Data', 'dir')
    mkdir('./Data')
    addpath(genpath('./Data'))
end

% %define the data set - define correctly segmented positions for WT and elp6
%main experiment
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

DataExp{14} = sprintf('ExpI5_pos19_Y1474_Elp6');
DataExp{15} = sprintf('ExpI5_pos20_Y1474_Elp6');
DataExp{16} = sprintf('ExpI5_pos22_Y1474_Elp6');
DataExp{17} = sprintf('ExpI5_pos23_Y1474_Elp6');
DataExp{18} = sprintf('ExpI5_pos24_Y1474_Elp6');
DataExp{19} = sprintf('ExpI5_pos25_Y1474_Elp6');
DataExp{20} = sprintf('ExpI5_pos27_Y1474_Elp6');
DataExp{21} = sprintf('ExpI5_pos30_Y1474_Elp6');
DataExp{22} = sprintf('ExpI5_pos31_Y1474_Elp6');
DataExp{23} = sprintf('ExpI5_pos32_Y1474_Elp6');

%replicate experiment 
DataExp{24} = sprintf('ExpL1_pos1_Y208_WT');
DataExp{25} = sprintf('ExpL1_pos5_Y208_WT');
DataExp{26} = sprintf('ExpL1_pos6_Y208_WT');
DataExp{27} = sprintf('ExpL1_pos7_Y208_WT');
DataExp{28} = sprintf('ExpL1_pos8_Y208_WT');
DataExp{29} = sprintf('ExpL2_pos1_Y208_WT');
DataExp{30} = sprintf('ExpL2_pos4_Y208_WT');
DataExp{31} = sprintf('ExpL2_pos6_Y208_WT');
DataExp{32} = sprintf('ExpL2_pos9_Y208_WT');
DataExp{33} = sprintf('ExpL2_pos10_Y208_WT');
DataExp{34} = sprintf('ExpL2_pos11_Y208_WT');
DataExp{35} = sprintf('ExpL2_pos12_Y208_WT');
DataExp{36} = sprintf('ExpI7_pos4_Y208_WT');

DataExp{37} = sprintf('ExpL1_pos9_Y1474_Elp6');
DataExp{38} = sprintf('ExpL1_pos10_Y1474_Elp6');
DataExp{39} = sprintf('ExpL1_pos11_Y1474_Elp6');
DataExp{40} = sprintf('ExpL1_pos12_Y1474_Elp6');
DataExp{41} = sprintf('ExpL1_pos13_Y1474_Elp6');
DataExp{42} = sprintf('ExpL1_pos14_Y1474_Elp6');
DataExp{43} = sprintf('ExpL1_pos15_Y1474_Elp6');
DataExp{44} = sprintf('ExpL1_pos16_Y1474_Elp6');
DataExp{45} = sprintf('ExpL2_pos14_Y1474_Elp6');
DataExp{46} = sprintf('ExpL2_pos15_Y1474_Elp6');
DataExp{47} = sprintf('ExpL2_pos16_Y1474_Elp6');
DataExp{48} = sprintf('ExpL2_pos17_Y1474_Elp6');
DataExp{49} = sprintf('ExpL2_pos18_Y1474_Elp6');
DataExp{50} = sprintf('ExpL2_pos19_Y1474_Elp6');
DataExp{51} = sprintf('ExpL2_pos20_Y1474_Elp6');
DataExp{52} = sprintf('ExpL2_pos21_Y1474_Elp6');
DataExp{53} = sprintf('ExpL2_pos22_Y1474_Elp6');
DataExp{54} = sprintf('ExpL2_pos23_Y1474_Elp6');
DataExp{55} = sprintf('ExpL2_pos24_Y1474_Elp6');
DataExp{56} = sprintf('ExpL2_pos25_Y1474_Elp6');
DataExp{57} = sprintf('ExpL2_pos26_Y1474_Elp6');
DataExp{58} = sprintf('ExpL2_pos27_Y1474_Elp6');
DataExp{59} = sprintf('ExpL2_pos28_Y1474_Elp6');

Data = DataExp;

Data = Data(~cellfun('isempty',Data));

for i = 1:length(Data)
    
    clearvars -except Data i
    clc;
    
    %load specific single cell data information
    loadData = sprintf('%s',Data{i});
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
%             BF(count) = segmentation.tcells1(icell).birthFrame;
            M(count) = segmentation.tcells1(icell).mother;
            N(count) = segmentation.tcells1(icell).N;
            
            for itime = 1:min(length(segmentation.tcells1(icell).Obj),length(segmentation.tcells1(icell).detectionFrame:320))
                %sort according to image as some are not sequential
                I{count}(itime) = segmentation.tcells1(icell).Obj(itime).image;
                usGFPrel{count}(itime) = segmentation.tcells1(icell).Obj(itime).fluoMean(2); %GFP channel
                usArea{count}(itime) = segmentation.tcells1(icell).Obj(itime).area;
            end
            
            %if cell detected after frame 320
            if min(length(segmentation.tcells1(icell).Obj),length(segmentation.tcells1(icell).detectionFrame:320))==0
                I{count}(itime) = NaN;
                usGFPrel{count}(itime) = NaN;
                usArea{count}(itime) = NaN;
            end
            
            if ~isempty(I{count})
                [valI,indI] = sort(I{count},'ascend');
                GFPrel{count} = usGFPrel{count}(indI);
                %get total GFP per cell per time
                GFPabs{count} = usGFPrel{count}(indI) .* usArea{count}(indI);
                Area{count} = usArea{count}(indI);
            else
                GFPrel{count} = NaN;
                %get total GFP per cell per time
                GFPabs{count} = NaN;
                Area{count} = NaN;
            end
            
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
        if S1{iS1}.LF >= 320
            %assure that for every frame we have the total GFP value
            if length(S1{iS1}.DF:320) == length(S1{iS1}.GFPabs)
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
    
    saveS = sprintf('./Melanie/Data/S%s',Data{i});
    save(saveS,'S');
    
end
