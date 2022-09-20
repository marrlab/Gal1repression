clearvars;
clc;

strain1 = 1;
strain2 = 1;
rep2 = 2;
rep1 = 2;


load(sprintf('scR1_strain%d_rep%d_model%d',strain1,rep1,1))
scR1 = scR;
load(sprintf('scR1_strain%d_rep%d_model%d',strain1,rep1,2))
scR2 = scR;

load(sprintf('scR1_strain%d_rep%d_model%d_var0',strain2,rep2,1))
scR1_var = scR;
load(sprintf('scR1_strain%d_rep%d_model%d_var0',strain2,rep2,2))
scR2_var = scR;

%get BIC values per total GFP trace
for i = 1:size(scR1,2)
    BIC1(i) = scR1(i).sol.BIC;
end
for i = 1:size(scR2,2)
    BIC2(i) = scR2(i).sol.BIC;
end
%identify which model better fits each trace
ind2 = find(BIC2-BIC1<-10); %model 2 best
ind1 = find(BIC2-BIC1>=-10);%model 1 best

%get BIC values per total GFP trace
for i = 1:size(scR1_var,2)
    BIC1_var(i) = scR1_var(i).sol.BIC;
end
for i = 1:size(scR2_var,2)
    BIC2_var(i) = scR2_var(i).sol.BIC;
end
%identify which model better fits each trace
ind2_var = find(BIC2_var-BIC1_var<-10); %model 2 best
ind1_var = find(BIC2_var-BIC1_var>=-10);%model 1 best

%get estimated parameters
for icell = 1:length(scR2)
    clear par
    par = 10.^(scR2(icell).sol.MS.par(:,1));
    Par1(icell,:) = par';
end

for icell = 1:length(scR2_var)
    clear par
    par = 10.^(scR2_var(icell).sol.MS.par(:,1));
    Par2(icell,:) = par';
end

subset_ind = [];
for iind2_var = 1:length(ind2_var)
    if isempty(find(ind2_var(iind2_var)==ind2))
        subset(iind2_var) = 0;
    else
        subset(iind2_var) = 1;
        subset_ind = [subset_ind,ind2_var(iind2_var)];
    end
end

disp(sprintf('%d (of %d) of previously %d cells are again found to be repressor cells.', sum(subset),length(ind2_var), length(ind2)))

Par1_1 = Par1(subset_ind,:);
Par2_1 = Par2(subset_ind,:);

figure 

c1 = [0,0,0];
c2 = [0,0,0];

for ipar = 1:4
    
    subplot(1,4,ipar)
    
    P1 = Par1_1(:,ipar)';
    P2 = Par2_1(:,ipar)';
        
    %[p_signtest,h_signtest,stats] = signtest(P1,P2);
    [h_ttest,p_ttest,] = ttest(P1,P2);
    Pval(ipar) = p_ttest;
    
    index = [ones(length(P1),1);2*ones(length(P2),1)];
    
    %plot with jitter
    a = -0.2;
    b = 0.2;
    r1 = (b-a).*rand(length(P1),1) + a;
    
    for i = 1:size(P1,2)
        line([(1+r1(i)),(2+r1(i))],[P1(i),P2(i)],'Color',[200,200,200]./255);
        hold on
    end
    r2 = r1;
    
%     sum(P1>P2)/length(P1)*100
    
    plot(1+r1,P1,'.','Color',c1,'Markersize',10)
    hold on
    plot(2+r2,P2,'^','Color',c2,'Markersize',4)
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
% set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 13 4])
% print('-dpdf','./Figures/FigS2C','-painters')