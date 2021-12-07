function compGFP0(strain1,rep1, iexp)

clearvars -except strain1 rep1 iexp
clc;

%strain1 = 1 - WT / = 2 - elp6
%rep1 = 1 - repression 1 / = 2 - repression 2

if iexp == 1
    load('NonDividing')
else
    load('NonDividing2')
end

for istrain = strain1
    for irep = rep1
        
        %load estaimted parameter sets
        if iexp == 1
            load(sprintf('scR_strain%d_rep%d_model1',istrain,irep));
            scR1_1 = scR;
            load(sprintf('scR_strain%d_rep%d_model2',istrain,irep));
            scR1_2 = scR;
        else
            load(sprintf('scR2_strain%d_rep%d_model1_truncated',istrain,irep));
            scR1_1 = scR;
            load(sprintf('scR2_strain%d_rep%d_model2_truncated',istrain,irep));
            scR1_2 = scR;
        end
        
        %get BIC values per total GFP trace
        for i = 1:size(scR1_1,2)
            BIC1_1(i) = scR1_1(i).sol.BIC;
        end
        for i = 1:size(scR1_2,2)
            BIC1_2(i) = scR1_2(i).sol.BIC;
        end
        %identify which model better fits each trace
        ind1_2 = find(BIC1_2-BIC1_1<-10); %model 2 best
        ind1_1 = find(BIC1_2-BIC1_1>=-10);%model 1 best
        
        %get estimated parameters
        for icell = 1:length(scR1_1)
            clear par
            par = 10.^(scR1_1(icell).sol.MS.par(:,1));
            Par1(icell,:) = par';
        end
        
        for icell = 1:length(scR1_2)
            clear par
            par = 10.^(scR1_2(icell).sol.MS.par(:,1));
            Par2(icell,:) = par';
        end
        
        for ipar = 1
            
            %compare the GFP0 (par=1) estimated parameters between total GFP
            %traces better fitted by the non-repressor (1) and repressor (2) model  
            P1 = Par1(ind1_1,ipar)';
            P2 = Par1(ind1_2,ipar)';
            
            index = [ones(length(P1),1);2*ones(length(P2),1)];
            
            %Mood's median test 
            [p_mediantest,tab,chi2] = mediantest(P1,P2);
            Pval(ipar) = p_mediantest;
            Pval
        end
        
        %save p-values of GFP0 comparison
        sol_GFP0.Pval = Pval;
        
        if iexp == 1
            save(sprintf('./Results/sol_GFP0_%d_%d',istrain,irep),'sol_GFP0')
        else
            save(sprintf('./Results/sol2_GFP0_%d_%d',istrain,irep),'sol_GFP0')
        end
    end
end

