%add path
addpath(genpath(pwd))

%% Figure 2B (left and right) - non-dividing total GFP traces for repressions r1 and r2 for WT

clearvars;
clc;

load('NonDividing1')
figure('visible','off');

%istrain = 1 - WT / = 2 - elp6
%irep = 1 - repression 1 / = 2 - repression 2
istrain = 1;

for irep = 1:2
    
    %define color according to strain and repression
    
    if irep == 1
        c = [117,157,233]./255;
    else
        c = [33,68,120]./255;
    end
    
    %plot the non-dividing cells of specified strain and repression
    %also plot mean total GFP trace and maximal mean total GFP value
    if irep == 1
        
        plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.r1(:,1:40)./10e6,'-','Color',c)
        hold on
        plot((1-1)*3/60:3/60:(40-1)*3/60,mean(NonDividing{istrain}.r1(:,1:40))./10e6,'-','Color','k')
        hold on
        indmax = find(mean(NonDividing{istrain}.r1(:,1:40))==max(mean(NonDividing{istrain}.r1(:,1:40))));
        time = (1-1)*3/60:3/60:(40-1)*3/60;
        plot(time(indmax),max(mean(NonDividing{istrain}.r1(:,1:40)))./10e6,'.','Color','k','Markersize',12)
        hold on
        
        %bootstrap the time to maximal mean total GFP
        for isample = 1:100000
            S = datasample(1:size(NonDividing{istrain}.r1,1),round(size(NonDividing{istrain}.r1,1)));
            indmax = find(mean(NonDividing{istrain}.r1(S,1:40))==max(mean(NonDividing{istrain}.r1(S,1:40))));
            T(isample) = time(indmax);
        end
        
        %get mean and std of time to maximal mean total GFP
        display(sprintf('Mean of WT time maximal mean total GFP repression r1 is %d', mean(T)))
        display(sprintf('Standard deviation of WT time maximal mean total GFP repression r1 is %d', std(T)))
        display(sprintf('Number of WT cells repression r1 is %d', size(NonDividing{istrain}.r1,1)))
        
        line([mean(T),mean(T)],[max(mean(NonDividing{istrain}.r1(:,1:40)))./10e6-0.4,...
            max(mean(NonDividing{istrain}.r1(:,1:40)))./10e6+0.4],'Color','k','Linewidth',1)
        hold on
        line([mean(T)-std(T),mean(T)+std(T)],[max(mean(NonDividing{istrain}.r1(:,1:40)))./10e6,...
            max(mean(NonDividing{istrain}.r1(:,1:40)))./10e6],'Color','k','Linewidth',1)
        hold on
        line([mean(T)+std(T),mean(T)+std(T)],[max(mean(NonDividing{istrain}.r1(:,1:40)))./10e6-0.2,...
            max(mean(NonDividing{istrain}.r1(:,1:40)))./10e6+0.2],'Color','k','Linewidth',1)
        hold on
        line([mean(T)-std(T),mean(T)-std(T)],[max(mean(NonDividing{istrain}.r1(:,1:40)))./10e6-0.2,...
            max(mean(NonDividing{istrain}.r1(:,1:40)))./10e6+0.2],'Color','k','Linewidth',1)
        
    else
        
        plot((1-1)*3/60:3/60:(40-1)*3/60,NonDividing{istrain}.r2(:,1:40)./10e6,'-','Color',c)
        hold on
        plot((1-1)*3/60:3/60:(40-1)*3/60,mean(NonDividing{istrain}.r2(:,1:40))./10e6,'-','Color','w')
        hold on
        indmax = find(mean(NonDividing{istrain}.r2(:,1:40))==max(mean(NonDividing{istrain}.r2(:,1:40))));
        time = (1-1)*3/60:3/60:(40-1)*3/60;
        plot(time(indmax),max(mean(NonDividing{istrain}.r2(:,1:40)))./10e6,'.','Color','w','Markersize',12)
        hold on
        
        %bootstrap the time to maximal mean total GFP
        for isample = 1:100000
            S = datasample(1:size(NonDividing{istrain}.r2,1),round(size(NonDividing{istrain}.r2,1)));
            indmax = find(mean(NonDividing{istrain}.r2(S,1:40))==max(mean(NonDividing{istrain}.r2(S,1:40))));
            T(isample) = time(indmax);
        end
        
        %get mean and std of time to maximal mean total GFP
        display(sprintf('Mean of WT time maximal mean total GFP repression r2 is %d', mean(T)))
        display(sprintf('Standard deviation of WT time maximal mean total GFP repression r2 is %d', std(T)))
        display(sprintf('Number of WT cells in repression r2 is %d', size(NonDividing{istrain}.r2,1)))
        
        line([mean(T),mean(T)],[max(mean(NonDividing{istrain}.r2(:,1:40)))./10e6-0.4,...
            max(mean(NonDividing{istrain}.r2(:,1:40)))./10e6+0.4],'Color','w','Linewidth',1)
        hold on
        line([mean(T)-std(T),mean(T)+std(T)],[max(mean(NonDividing{istrain}.r2(:,1:40)))./10e6,...
            max(mean(NonDividing{istrain}.r2(:,1:40)))./10e6],'Color','w','Linewidth',1)
        hold on
        line([mean(T)+std(T),mean(T)+std(T)],[max(mean(NonDividing{istrain}.r2(:,1:40)))./10e6-0.2,...
            max(mean(NonDividing{istrain}.r2(:,1:40)))./10e6+0.2],'Color','w','Linewidth',1)
        hold on
        line([mean(T)-std(T),mean(T)-std(T)],[max(mean(NonDividing{istrain}.r2(:,1:40)))./10e6-0.2,...
            max(mean(NonDividing{istrain}.r2(:,1:40)))./10e6+0.2],'Color','w','Linewidth',1)    
        
    end
    
    ylabel('total GFP (a.u)')
    if irep == 1
        yticks([0,2,4])
        ylim([0,5])
    else
        yticks([0,5,10,15])
        ylim([0,15])
    end
    xlabel('repression time (h)')
    xticks([0,1,2])
    box off
    set(gca,'linewidth',1.02)
    set(gca,'FontSize',11)
    set(gca,'FontName','Arial')
    xlim([0,2])
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 7 5])
    if irep == 1
        print('-dpdf','./Figures/Fig2Bleft','-painters')
    else
        print('-dpdf','./Figures/Fig2Bright','-painters')
    end
end
