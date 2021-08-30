
clearvars;
clc;

load('/Users/lea.schuh/Documents/PhD/IFE/Gal1repression/Data/NonDividing')

figure('visible','off');
count = 1;
S = 100000;

for istrain = 1:2
    for irep = 1:2
        
        for isample=1:S

            %get maximal point of all cells 
            if irep==1
                subset = datasample(1:size(NonDividing{istrain}.I5r1,1),round(0.25*size(NonDividing{istrain}.I5r1,1)));
                indmax = find(mean(NonDividing{istrain}.I5r1(subset,:))==max(mean(NonDividing{istrain}.I5r1(subset,:))));
            else
                subset = datasample(1:size(NonDividing{istrain}.I5r2,1),round(0.25*size(NonDividing{istrain}.I5r2,1)));
                indmax = find(mean(NonDividing{istrain}.I5r2(subset,:))==max(mean(NonDividing{istrain}.I5r2(subset,:))));
            end
            time = (1-1)*3/60:3/60:(40-1)*3/60;
            T(count,isample) = time(indmax);
        end

        a = -0.2;
        b = 0.2;
        r1 = (b-a).*rand(S,1) + a;

        line([count-0.4,count+0.4],[mean(T(count,:)),mean(T(count,:))],'Color','k','Linewidth',1)
        hold on
        line([count,count],[mean(T(count,:))-std(T(count,:)),mean(T(count,:))+std(T(count,:))],'Color','k','Linewidth',1)
        hold on
        line([count-0.2,count+0.2],[mean(T(count,:))+std(T(count,:)),mean(T(count,:))+std(T(count,:))],'Color','k','Linewidth',1)
        hold on
        line([count-0.2,count+0.2],[mean(T(count,:))-std(T(count,:)),mean(T(count,:))-std(T(count,:))],'Color','k','Linewidth',1)
        count = count+1;
    end
end

xlim([0,5])
set(gca,'FontSize',10)
ylim([0, 2])
box off
set(gca,'linewidth',1.02)
set(gca,'FontSize',11)
set(gca,'FontName','Arial')
xlabel('data set')
ylabel('maximal mean total GFP')

%save figure
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 5.5 5])
print('-dpdf','./Figures/Fig2G','-painters')