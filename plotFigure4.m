%% Figure 4 - comparison of estimated parameters WT

%strain = 1 - WT / = 2 - elp6
%rep = 1 - repression 1 / = 2 - repression 2

strain1 = 1;
rep1 = 1;
strain2 = 1;
rep2 = 2;
paired = true;

clearvars -except rep1 strain1 rep2 strain2 paired
clc;

%load all estimated parameter sets for both models and repressions
load(sprintf('scR1_strain%d_rep%d_model%d',strain1,rep1,1))
scR1_1 = scR;
load(sprintf('scR1_strain%d_rep%d_model%d',strain1,rep1,2))
scR1_2 = scR;

load(sprintf('scR1_strain%d_rep%d_model%d',strain2,rep2,1))
scR2_1 = scR;
load(sprintf('scR1_strain%d_rep%d_model%d',strain2,rep2,2))
scR2_2 = scR;

%extract BIC values for all single-cell trajectories for data set 1
for i = 1:size(scR1_1,2)
    BIC1_1(i) = scR1_1(i).sol.BIC;
end
for i = 1:size(scR1_2,2)
    BIC1_2(i) = scR1_2(i).sol.BIC;
end

%decide whether single-cell requires repressor model or not for data set 1
ind1_2 = find(BIC1_2-BIC1_1<-10); %model 2 best
ind1_1 = find(BIC1_2-BIC1_1>=-10);%model 1 best

%extract BIC values for all single-cell trajectories for data set 2
for i = 1:size(scR2_1,2)
    BIC2_1(i) = scR2_1(i).sol.BIC;
end
for i = 1:size(scR2_2,2)
    BIC2_2(i) = scR2_2(i).sol.BIC;
end

%decide whether single-cell requires repressor model or not for data set 2
ind2_2 = find(BIC2_2-BIC2_1<-10);
ind2_1 = find(BIC2_2-BIC2_1>=-10);

%get the data sets for both experiments
load('NonDividing1')

%extract data from non-dividing cells structure
if rep1 == 1
    data1 = NonDividing{strain1}.r1;
    momID1 = NonDividing{strain1}.momIDr1;
    mompos1 = NonDividing{strain1}.momposr1;
else
    data1 = NonDividing{strain1}.r2;
    momID1 = NonDividing{strain1}.momIDr2;
    mompos1 = NonDividing{strain1}.momposr2;
end

if rep2 == 1
    data2 = NonDividing{strain2}.r1;
    momID2 = NonDividing{strain2}.momIDr1;
    mompos2 = NonDividing{strain2}.momposr1;
else
    data2 = NonDividing{strain2}.r2;
    momID2 = NonDividing{strain2}.momIDr2;
    mompos2 = NonDividing{strain2}.momposr2;
end

%extract the estimated parameter sets per cell and data set
for icell = 1:length(scR1_2)
    clear par
    par = 10.^(scR1_2(icell).sol.MS.par(:,1));
    Par1(icell,:) = par';
end

for icell = 1:length(scR2_2)
    clear par
    par = 10.^(scR2_2(icell).sol.MS.par(:,1));
    Par2(icell,:) = par';
end

%reduce parameter sets to the cells requiring the repressor model
Par1 = Par1(ind1_2,:);
Par2 = Par2(ind2_2,:);

%if comparison between cells in repressions 1 and 2
if paired == 1
    %extract the cell IDs and positions of cells requiring repressor model in
    %repressions 1 and 2
    momInfo1 = [momID1(ind1_2)',mompos1(ind1_2)'];
    momInfo2 = [momID2(ind2_2)',mompos2(ind2_2)'];
    [~,index_momInfo1,index_momInfo2] = intersect(momInfo1,momInfo2,'rows');
    Par1 = Par1(index_momInfo1,:);
    Par2 = Par2(index_momInfo2,:);
end

%determine the color according to strain and repression
if rep1 == 1
    c1 = [117,157,233]./255;
else
    c1 = [33,68,120]./255;
end

if rep2 == 1
    c2 = [117,157,233]./255;
else
    c2 = [33,68,120]./255;
end

figure('visible','off');
%for each of the estimated parameters GFP0, delay, rprod and rdeg do
for ipar = 1:4
    
    clear P1_1 index h_kstest2 p_kstest2 h_ranksum p_ranksum
    
    subplot(1,4,ipar)
    
    %get correct estimnated parameters
    P1 = Par1(:,ipar)';
    P2 = Par2(:,ipar)';
    index = [ones(length(P1),1);2*ones(length(P2),1)];
    
    %plot with jitter
    a = -0.2;
    b = 0.2;
    r1 = (b-a).*rand(length(P1),1) + a;
    
    if paired
        for i = 1:size(P1,2)
            line([(1+r1(i)),(2+r1(i))],[P1(i),P2(i)],'Color',[200,200,200]./255);
            hold on
        end
        r2 = r1;
    else
        r2 = (b-a).*rand(length(P1),1) + a;
    end
    
    sum(P1>P2)/length(P1)*100
    
    plot(1+r1,P1,'.','Color',c1,'Markersize',10)
    hold on
    plot(2+r2,P2,'.','Color',c2,'Markersize',10)
    hold on
    line([0.6,1.4],[median(P1),median(P1)],'Color','k','Linewidth',2)
    hold on
    line([1.6,2.4],[median(P2),median(P2)],'Color','k','Linewidth',2)
    hold on
    
    xlim([0,3])
    set(gca,'FontSize',10,'XTickLabelRotation',90)
    ylim([0 inf])
    box off
    set(gca,'linewidth',1.02)
    set(gca,'FontSize',11)
    set(gca,'FontName','Arial')
    
end

%save figure
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 13 4])
print('-dpdf','./Figures/Fig4','-painters')
