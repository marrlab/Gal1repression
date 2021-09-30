function comp_GFP0_tdelay(rep1,strain1,rep2,strain2)

%rep1 = 1 - repression 1 / = 2 - repression 2
%strain1 = 1 - WT / = 2 - elp6
%rep2 = 1 - repression 1 / = 2 - repression 2
%strain2 = 1 - WT / = 2 - elp6

%not paired as we look at all cells for regression fit
paired = 0;

clearvars -except rep1 strain1 rep2 strain2 paired
clc;

%load all estimated parameter sets for both models and repressions
load(sprintf('scR_strain%d_rep%d_model%d',strain1,rep1,1))
scR1_1 = scR;
load(sprintf('scR_strain%d_rep%d_model%d',strain1,rep1,2))
scR1_2 = scR;

load(sprintf('scR_strain%d_rep%d_model%d',strain2,rep2,1))
scR2_1 = scR;
load(sprintf('scR_strain%d_rep%d_model%d',strain2,rep2,2))
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
load('NonDividing')

if rep1 == 1
    data1 = NonDividing{strain1}.I5r1;
    momID1 = NonDividing{strain1}.I5momIDr1;
    mompos1 = NonDividing{strain1}.I5momposr1;
    numdaughters1 = NonDividing{strain1}.I5momcountr1;
    IDdaughters1 = NonDividing{strain1}.I5daughtersIDr1;
    countmom1 = 1;
    for imom1 = 1:length(momID1)
        M1{imom1} = IDdaughters1(countmom1:(countmom1+numdaughters1(imom1)-1));
        countmom1 = countmom1+numdaughters1(imom1);
    end
else
    data1 = NonDividing{strain1}.I5r2;
    momID1 = NonDividing{strain1}.I5momIDr2;
    mompos1 = NonDividing{strain1}.I5momposr2;
    numdaughters1 = NonDividing{strain1}.I5momcountr2;
    IDdaughters1 = NonDividing{strain1}.I5daughtersIDr2;
    countmom1 = 1;
    for imom1 = 1:length(momID1)
        M1{imom1} = IDdaughters1(countmom1:(countmom1+numdaughters1(imom1)-1));
        countmom1 = countmom1+numdaughters1(imom1);
    end
end

if rep2 == 1
    data2 = NonDividing{strain2}.I5r1;
    momID2 = NonDividing{strain2}.I5momIDr1;
    mompos2 = NonDividing{strain2}.I5momposr1;
    numdaughters2 = NonDividing{strain1}.I5momcountr1;
    IDdaughters2 = NonDividing{strain1}.I5daughtersIDr1;
    countmom2 = 1;
    for imom2 = 1:length(momID2)
        M2{imom2} = IDdaughters2(countmom2:(countmom2+numdaughters2(imom2)-1));
        countmom2 = countmom2+numdaughters2(imom2);
    end
else
    data2 = NonDividing{strain2}.I5r2;
    momID2 = NonDividing{strain2}.I5momIDr2;
    mompos2 = NonDividing{strain2}.I5momposr2;
    numdaughters2 = NonDividing{strain1}.I5momcountr2;
    IDdaughters2 = NonDividing{strain1}.I5daughtersIDr2;
    countmom2 = 1;
    for imom2 = 1:length(momID2)
        M2{imom2} = IDdaughters2(countmom2:(countmom2+numdaughters2(imom2)-1));
        countmom2 = countmom2+numdaughters2(imom2);
    end
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

if rep1 < rep2
    %get linear regression fits for data set 1
    c1 = polyfit(Par1(:,1),Par1(:,2),1);
    % Display evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c1(1)) '*x + ' num2str(c1(2))])
    
    %get linear regression fits for data set 2
    c2 = polyfit(Par2(:,1),Par2(:,2),1);
    % Display evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c2(1)) '*x + ' num2str(c2(2))])
    
    %test whether repression slope of data set 1 == 0
    X=[Par1(:,1) ones(size(Par1(:,1),1),1)];
    y=Par1(:,2);
    [b,bint,r,rint,stats]=regress(y,X);
    stats1 = stats(3);
    
    %test whether repression slope of data set 2 == 0
    X=[Par2(:,1) ones(size(Par2(:,1),1),1)];
    y=Par2(:,2);
    [b,bint,r,rint,stats]=regress(y,X);
    stats2 = stats(3);
    
    sol_GFP0_tdelay.c1 = c1;
    sol_GFP0_tdelay.c2 = c2;
    sol_GFP0_tdelay.stats1 = stats1;
    sol_GFP0_tdelay.stats2 = stats2;
    
end

if rep1 == rep2
    
     %get linear regression fits for data set 1
    c1 = polyfit(Par1(:,1),Par1(:,2),1);
    % Display evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c1(1)) '*x + ' num2str(c1(2))])
    
    %get linear regression fits for data set 2
    c2 = polyfit(Par2(:,1),Par2(:,2),1);
    % Display evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c2(1)) '*x + ' num2str(c2(2))])
    
    C = table([Par2(:,2);Par1(:,2)],[Par2(:,1);Par1(:,1)],[zeros(length(Par2(:,1)),1);...
        ones(length(Par1(:,1)),1)],'VariableNames',{'t','GFP0','rep'});
    C.rep = categorical(C.rep);
    
    fit = fitlm(C,'t~GFP0*rep');
    sol_GFP0_tdelay.c1 = c1;
    sol_GFP0_tdelay.c2 = c2;
    sol_GFP0_tdelay.fit = fit;

end

save(sprintf('./Results/sol_GFP0_tdelay_%d_%d_%d_%d',rep1,strain1,rep2,strain2),'sol_GFP0_tdelay')

end