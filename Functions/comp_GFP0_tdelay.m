function comp_GFP0_tdelay(rep1,strain1,rep2,strain2,iexp)

%rep1 = 1 - repression 1 / = 2 - repression 2
%strain1 = 1 - WT / = 2 - elp6
%rep2 = 1 - repression 1 / = 2 - repression 2
%strain2 = 1 - WT / = 2 - elp6
%iexp = 1 - main experiment / = 2 - replicate experiment

clearvars -except rep1 strain1 rep2 strain2 paired iexp
clc;

%load all estimated parameter sets for both models and repressions
if iexp == 1
    load(sprintf('scR1_strain%d_rep%d_model%d',strain1,rep1,1))
    scR1_1 = scR;
    load(sprintf('scR1_strain%d_rep%d_model%d',strain1,rep1,2))
    scR1_2 = scR;
    
    load(sprintf('scR1_strain%d_rep%d_model%d',strain2,rep2,1))
    scR2_1 = scR;
    load(sprintf('scR1_strain%d_rep%d_model%d',strain2,rep2,2))
    scR2_2 = scR;
else
    load(sprintf('scR2_strain%d_rep%d_model%d',strain1,rep1,1))
    scR1_1 = scR;
    load(sprintf('scR2_strain%d_rep%d_model%d',strain1,rep1,2))
    scR1_2 = scR;
    
    load(sprintf('scR2_strain%d_rep%d_model%d',strain2,rep2,1))
    scR2_1 = scR;
    load(sprintf('scR2_strain%d_rep%d_model%d',strain2,rep2,2))
    scR2_2 = scR;
end

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

%split into low (<3) and high (>3) inducing cells
Par1a = Par1(Par1(:,1)<3,:);
Par2a = Par2(Par2(:,1)<3,:);

Par1b = Par1(Par1(:,1)>3,:);
Par2b = Par2(Par2(:,1)>3,:);

%get linear regression fits for data set 1 < 3
c1 = polyfit(Par1a(:,1),Par1a(:,2),1);
% Display evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(c1(1)) '*x + ' num2str(c1(2))])

%get linear regression fits for data set 2 < 3
c2 = polyfit(Par2a(:,1),Par2a(:,2),1);
% Display evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(c2(1)) '*x + ' num2str(c2(2))])

%test whether repression slope of data set 1 == 0
X=[Par1a(:,1) ones(size(Par1a(:,1),1),1)];
y=Par1a(:,2);
[b,bint,r,rint,stats]=regress(y,X);
stats1 = stats(3);

%test whether repression slope of data set 2 == 0
X=[Par2a(:,1) ones(size(Par2a(:,1),1),1)];
y=Par2a(:,2);
[b,bint,r,rint,stats]=regress(y,X);
stats2 = stats(3);

sol_GFP0_tdelay.c1a = c1;
sol_GFP0_tdelay.c2a = c2;
sol_GFP0_tdelay.stats1a = stats1;
sol_GFP0_tdelay.stats2a = stats2;

%get linear regression fits for data set 1 > 3
c1 = polyfit(Par1b(:,1),Par1b(:,2),1);
% Display evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(c1(1)) '*x + ' num2str(c1(2))])

%get linear regression fits for data set 2 > 3
c2 = polyfit(Par2b(:,1),Par2b(:,2),1);
% Display evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(c2(1)) '*x + ' num2str(c2(2))])

%test whether repression slope of data set 1 == 0
X=[Par1b(:,1) ones(size(Par1b(:,1),1),1)];
y=Par1b(:,2);
[b,bint,r,rint,stats]=regress(y,X);
stats1 = stats(3);

%test whether repression slope of data set 2 == 0
X=[Par2b(:,1) ones(size(Par2b(:,1),1),1)];
y=Par2b(:,2);
[b,bint,r,rint,stats]=regress(y,X);
stats2 = stats(3);

sol_GFP0_tdelay.c1b = c1;
sol_GFP0_tdelay.c2b = c2;
sol_GFP0_tdelay.stats1b = stats1;
sol_GFP0_tdelay.stats2b = stats2;

if iexp == 1
    save(sprintf('./Results/sol1_GFP0_tdelay_%d_%d_%d_%d',rep1,strain1,rep2,strain2),'sol_GFP0_tdelay')
else
    save(sprintf('./Results/sol2_GFP0_tdelay_%d_%d_%d_%d',rep1,strain1,rep2,strain2),'sol_GFP0_tdelay')
end

end