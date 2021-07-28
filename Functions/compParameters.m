function compParameters(rep1,strain1,rep2,strain2,paired)

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

%for each of the estimated parameters GFP0, delay, rprod and rdeg do
for ipar = 1:4
    
    clear P1_1 index h_kstest2 p_kstest2 h_ranksum p_ranksum
    
    %get correct estimnated parameters 
    P1 = Par1(:,ipar)';
    P2 = Par2(:,ipar)';
    index = [ones(length(P1),1);2*ones(length(P2),1)];
    
    if paired == 0
        [p_mediantest,tab,chi2] = mediantest(P1,P2);
        Pval(ipar) = p_mediantest;
    else
        [p_signtest,h_signtest,stats] = signtest(P1,P2);
        Pval(ipar) = p_signtest;
    end
    
    M(ipar,1) = median(P1);
    M(ipar,2) = median(P2);

end

sol_Par.Pval = Pval;
sol_Par.Median = M;


save(sprintf('./Results/sol_Par_%d_%d_%d_%d_%d',rep1,strain1,rep2,strain2,paired),'sol_Par')

end