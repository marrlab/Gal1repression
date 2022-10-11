function sD = getSimulation(rprod, rdeg)

figure

timepoints1 = 0:3/60:0.5; 
timepoints2 = 0:3/60:1.5-3/60;

for irun = 1:1
    %rate constants
    p.rprod = rprod;
    p.rdeg = rdeg;
    
    %initial state
    tspan1 = [0, 0.5]; %hours
    x0 = 100;
    
    %specify reaction network
    pfun = @propensities_state;
    pfun2 = @propensities_state2;
    
    S = [1, -1];

    stoich_matrix = S';
    
    % simulation
    [t,x] = firstReactionMethod(stoich_matrix, pfun, tspan1, x0, p);
    
    for i = 1:length(timepoints1)
        ind = find(t<=timepoints1(i),1,'last');
        sD1{1}(irun,i) = x(ind,end);
    end
    
    stoich_matrix2 = -1;
    
    x02 = sD1{1}(end);
    
    p2.rdeg = rdeg;
    
    tspan2 = [0, 1.5]; %hours
    
    [t2,x2] = firstReactionMethod(stoich_matrix2, pfun2, tspan2, x02, p2);
    
    for i = 1:length(timepoints2)
        ind = find(t2<=timepoints2(i),1,'last');
        sD2{1}(irun,i) = x2(ind,end);
    end
    % plot
%         figure();
        sD1{1} = sD1{1}+normrnd(0,5,1,length(sD1{1}));
        plot(timepoints1,sD1{1}, ':', 'Color','k'); 
        hold on
        sD2{1} = sD2{1}+normrnd(0,5,1,length(sD2{1}));
        plot(timepoints2+0.5,[sD1{1}(end),sD2{1}(2:end)], ':', 'Color','k');
        set(gca,'XLim',[0,2]);
        set(gca,'YLim',[0,Inf]);
        xlabel('time (s)');
        ylabel('molecules');
        legend({'mRNA'});
        hold on
%         if irun == 100
% %             plot(timepoints, rprod/rdeg*(1-exp(-rdeg*timepoints))+x0*exp(-rdeg*timepoints), 'Color', 'k','Linewidth',2);
%             plot(timepoints, x0*exp(-rdeg*timepoints), 'Color', 'k','Linewidth',2);
%             hold on
%             plot(timepoints,mean(sD{1}),'Color','r','Linewidth',2)
%             
%             figure
%             plot(timepoints,var(sD{1}),'Color','r','Linewidth',2)
%             hold on
%             %plot(timepoints, 1/rdeg.*exp(-2*rdeg*timepoints).*(exp(rdeg*timepoints)-1)...
%             %    .*(rdeg*x0+rprod*exp(rdeg*timepoints)), 'Color', 'k','Linewidth',2);
%             plot(timepoints, x0.*exp(-2*rdeg*timepoints).*(exp(rdeg*timepoints)-1), 'Color', 'k','Linewidth',2);
%         end

    sD = [sD1{1},sD2{1}(2:end)];
    
end
end

function a = propensities_state(x, p)

a(1) = p.rprod;
a(2) = p.rdeg*x;

a = a';

end

function a = propensities_state2(x, p)

a(1) = p.rdeg*x;

a = a';

end
